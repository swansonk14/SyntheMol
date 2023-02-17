"""Contains classes and functions for performing a tree search to generate molecules combinatorially."""
import itertools
import json
import math
import pickle
from collections import Counter
from datetime import datetime
from functools import cache, cached_property, partial
from pathlib import Path
from typing import Any, Callable, Literal, Optional

import numpy as np
import pandas as pd
import torch
from rdkit import Chem
from sklearn.metrics import pairwise_distances
from tap import Tap
from tqdm import trange

from chem_utils.molecular_fingerprints import compute_fingerprint, compute_fingerprints
from chemprop.utils import load_checkpoint
from reactions import Reaction, set_allowed_reaction_smiles
from real_reactions import REAL_REACTIONS
from synnet_reactions import SYNNET_REACTIONS
from train_model import sklearn_predict


class Args(Tap):
    # Required args
    fragment_path: Path
    """Path to CSV file containing molecular building blocks."""
    reaction_to_reagents_path: Path
    """Path to JSON file containing mapping from REAL reactions to allowed reagents."""
    search_type: Literal['random', 'greedy', 'mcts']
    """Type of search to perform."""
    save_dir: Path
    """Path to directory where the generated molecules will be saved."""

    # Model args
    model_path: Optional[Path] = None
    """Path to a directory of model checkpoints or to a specific PKL or PT file containing a trained model."""
    model_type: Optional[Literal['rf', 'mlp', 'chemprop']] = None
    """Type of model to train. 'rf' = random forest. 'mlp' = multilayer perceptron."""
    fingerprint_type: Optional[Literal['morgan', 'rdkit']] = None
    """Type of fingerprints to use as input features."""
    fragment_to_model_score_path: Optional[Path] = None
    """Path to JSON file containing a dictionary mapping fragments to model prediction scores."""

    # Data args
    smiles_column: str = 'smiles'
    """Name of the column containing SMILES."""
    fragment_id_column: str = 'Reagent_ID'
    """Name of the column in fragment_path that contains fragment IDs."""

    # Train similarity args
    train_similarity: bool = False
    """Whether to incorporate similarity to the train set as part of the MCTS score."""
    train_path: Optional[Path] = None
    """Path to CSV file containing SMILES for the molecules in the training set.
    Only used if train_similarity is True."""
    fragment_to_train_similarity_path: Optional[Path] = None
    """Path to JSON file containing a dictionary mapping fragments to train similarities
     (top 10 nearest neighbor Tanimoto). Only used if train_similarity is True."""
    train_hits_path: Optional[Path] = None
    """Path to CSV file containing SMILES for the active molecules in the training set.
    Only used if train_similarity is True."""
    fragment_to_train_hits_similarity_path: Optional[Path] = None
    """Path to JSON file containing a dictionary mapping fragments to train hits similarities
     (nearest neighbor Tversky). Only used if train_similarity is True."""

    # Search args
    max_reactions: int = 1
    """Maximum number of reactions that can be performed to expand fragments into molecules."""
    n_rollout: int = 3
    """The number of times to run the tree search."""
    c_puct: float = 10.0
    """The hyperparameter that encourages exploration."""
    num_expand_nodes: Optional[int] = None
    """The number of tree nodes to expand when extending the child nodes in the search tree."""
    synnet_rxn: bool = False
    """Whether to include SynNet reactions in addition to REAL reactions."""
    binarize_scoring: float = 0
    """If > 0, then molecule scores are binarized based on whether they are >= this threshold."""
    noise: bool = False
    """Whether to add uniformly random noise to molecule scores."""
    noise_limit: float = 0.2
    """If noise is True, this is the maximum noise added to molecule scores."""
    rng_seed: int = 0
    """Seed for random number generators."""
    fragment_diversity: bool = False
    """Whether to encourage the use of diverse fragments by modifying the score."""
    debug: bool = False
    """Whether to print out additional statements for debugging."""
    count_nodes: bool = False
    """Whether to count the number of nodes in the tree search.
    Note: The memory usage will grow linearly with the number of nodes searched since
    a set of all nodes must be stored in order to accurately count the number of nodes."""


class TreeNode:
    """A node in a tree search representing a step in the molecule construction process."""

    def __init__(self,
                 c_puct: float,
                 scoring_fn: Callable[[str], float],
                 node_id: Optional[int] = None,
                 fragments: Optional[tuple[str]] = None,
                 unique_reagents: Optional[set[int]] = None,
                 construction_log: Optional[tuple[dict[str, Any]]] = None,
                 rollout_num: Optional[int] = None) -> None:
        """Initializes the TreeNode object.

        :param c_puct: The hyperparameter that encourages exploration.
        :param scoring_fn: A function that takes as input a SMILES representing a molecule and returns a score.
        :param node_id: The ID of the node, which should correspond to the order in which the node was visited.
        :param fragments: A tuple of SMILES containing the fragments for the next reaction.
                         The first element is the currently constructed molecule while the remaining elements
                         are the fragments that are about to be added.
        :param unique_reagents: A set of reagent indices used in this node.
        :param construction_log: A tuple of dictionaries containing information about each reaction.
        :param rollout_num: The number of the rollout on which this node was created.
        """
        self.c_puct = c_puct
        self.scoring_fn = scoring_fn
        self.node_id = node_id
        self.fragments = fragments if fragments is not None else tuple()
        self.unique_reagents = unique_reagents if unique_reagents is not None else set()
        self.construction_log = construction_log if construction_log is not None else tuple()
        # TODO: maybe change sum (really mean) to max since we don't care about finding the best leaf node, just the best node along the way?
        self.W = 0.0  # The sum of the leaf node values for leaf nodes that descend from this node.
        self.N = 0  # The number of times this node has been expanded.
        self.rollout_num = rollout_num

    @classmethod
    def compute_score(cls, fragments: tuple[str], scoring_fn: Callable[[str], float]) -> float:
        """Computes the score of the fragments.

        :param fragments: A tuple of SMILES containing the fragments for the next reaction.
                         The first element is the currently constructed molecule while the remaining elements
                         are the fragments that are about to be added.
        :param scoring_fn: A function that takes as input a SMILES representing a molecule and returns a score.
        :return: The score of the fragments.
        """
        # TODO: change this!!! to weight the fragments differently
        return sum(scoring_fn(fragment) for fragment in fragments) / len(fragments) if len(fragments) > 0 else 0.0

    @cached_property
    def P(self) -> float:
        """The property score of this node. (Note: The value is cached, so it assumes the node is immutable.)"""
        return self.compute_score(fragments=self.fragments, scoring_fn=self.scoring_fn)

    def Q(self) -> float:
        """Value that encourages exploitation of nodes with high reward."""
        return self.W / self.N if self.N > 0 else 0.0

    def U(self, n: int) -> float:
        """Value that encourages exploration of nodes with few visits."""
        return self.c_puct * self.P * math.sqrt(1 + n) / (1 + self.N)

    @property
    def num_fragments(self) -> int:
        """Returns the number of fragments for the upcoming reaction.

        :return: The number of fragments for the upcoming reaction.
        """
        return len(self.fragments)

    @property
    def num_reactions(self) -> int:
        """Returns the number of reactions used so far to generate the current molecule.

        :return: The number of reactions used so far to generate the current molecule.
        """
        return len(self.construction_log)

    def __hash__(self) -> int:
        return hash(self.fragments)

    def __eq__(self, other: Any) -> bool:
        if not isinstance(other, TreeNode):
            return False

        return self.fragments == other.fragments


class TreeSearcher:
    """A class that runs a tree search to generate high scoring molecules."""

    def __init__(self,
                 search_type: Literal['random', 'greedy', 'mcts'],
                 fragment_smiles_to_id: dict[str, int],
                 max_reactions: int,
                 scoring_fn: Callable[[str], float],
                 n_rollout: int,
                 c_puct: float,
                 num_expand_nodes: Optional[int],
                 reactions: list[Reaction],
                 rng_seed: int,
                 fragment_diversity: bool,
                 debug: bool,
                 count_nodes: bool) -> None:
        """Creates the TreeSearcher object.

        :param search_type: Type of search to perform.
        :param fragment_smiles_to_id: A dictionary mapping fragment SMILES to their IDs.
        :param max_reactions: The maximum number of reactions to use to construct a molecule.
        :param scoring_fn: A function that takes as input a SMILES representing a molecule and returns a score.
        :param n_rollout: The number of times to run the tree search.
        :param c_puct: The hyperparameter that encourages exploration.
        :param num_expand_nodes: The number of tree nodes to expand when extending the child nodes in the search tree.
        :param reactions: A list of chemical reactions that can be used to combine molecular fragments.
        :param rng_seed: Seed for the random number generator.
        :param fragment_diversity: Whether to encourage the use of diverse fragments by modifying the score.
        :param debug: Whether to print out additional statements for debugging.
        :param count_nodes: Whether to count the number of nodes in the tree search.
                                Note: The memory usage will grow linearly with the number of nodes searched since
                                a set of all nodes must be stored in order to accurately count the number of nodes.
        """
        self.search_type = search_type
        self.fragment_smiles_to_id = fragment_smiles_to_id
        # TODO: maybe just allow all building blocks?
        self.all_fragments = list(dict.fromkeys(
            fragment
            for reaction in REAL_REACTIONS
            for reagent in reaction.reagents
            for fragment in reagent.allowed_smiles
        ))
        self.max_reactions = max_reactions
        self.scoring_fn = scoring_fn
        self.n_rollout = n_rollout
        self.c_puct = c_puct
        self.num_expand_nodes = num_expand_nodes
        self.reactions = reactions
        self.rng = np.random.default_rng(seed=rng_seed)
        self.fragment_diversity = fragment_diversity
        self.debug = debug

        self.rollout_num = 0
        self.TreeNodeClass = partial(TreeNode, c_puct=c_puct, scoring_fn=scoring_fn)
        self.root = self.TreeNodeClass(node_id=1, rollout_num=0)
        self.state_map: dict[TreeNode, TreeNode] = {self.root: self.root}
        self.reagent_counts = Counter()

        self.node_fragments_set = {self.root.fragments} if count_nodes else None

    def random_choice(self, array: list[Any], size: Optional[int] = None, replace: bool = True) -> Any:
        if size is None:
            return array[self.rng.integers(len(array))]

        return [
            array[i] for i in self.rng.choice(len(array), size=size, replace=replace)
        ]

    # TODO: documentation
    @classmethod
    def get_reagent_matches_per_mol(cls, reaction: Reaction, fragments: tuple[str]) -> list[list[int]]:
        return [
            [
                reagent_index
                for reagent_index, reagent in enumerate(reaction.reagents)
                if reagent.has_match(fragment)
            ]
            for fragment in fragments
        ]

    # TODO: documentation
    def get_next_fragments(self, fragments: tuple[str]) -> list[str]:
        # Initialize list of allowed fragments
        available_fragments = []

        # Loop through each reaction
        for reaction in self.reactions:
            # Get indices of the reagents in this reaction
            reagent_indices = set(range(reaction.num_reagents))

            # Skip reaction if there's no room to add more reagents
            if len(fragments) >= reaction.num_reagents:
                continue

            # For each mol, get a list of indices of reagents it matches
            reagent_matches_per_mol = self.get_reagent_matches_per_mol(reaction=reaction, fragments=fragments)

            # Loop through products of reagent indices that the mols match to
            # and for each product, if it matches to all separate reagents,
            # then include the missing reagents in the set of unfilled reagents
            for matched_reagent_indices in itertools.product(*reagent_matches_per_mol):
                matched_reagent_indices = set(matched_reagent_indices)

                if len(matched_reagent_indices) == len(fragments):
                    for index in sorted(reagent_indices - matched_reagent_indices):
                        available_fragments += reaction.reagents[index].allowed_smiles

        # Remove duplicates but maintain order for reproducibility and avoid sorting sets for speed
        # Note: requires Python 3.7+ for ordered dictionaries
        available_fragments = list(dict.fromkeys(available_fragments))

        return available_fragments

    def get_reactions_for_fragments(self, fragments: tuple[str]) -> list[tuple[Reaction, dict[str, int]]]:
        matching_reactions = []

        for reaction in self.reactions:
            if len(fragments) != reaction.num_reagents:
                continue

            # For each mol, get a list of indices of reagents it matches
            reagent_matches_per_mol = self.get_reagent_matches_per_mol(reaction=reaction, fragments=fragments)

            # Include every assignment of fragments to reagents that fills all the reagents
            for matched_reagent_indices in itertools.product(*reagent_matches_per_mol):
                if len(set(matched_reagent_indices)) == reaction.num_reagents:
                    fragment_to_reagent_index = dict(zip(fragments, matched_reagent_indices))
                    matching_reactions.append((reaction, fragment_to_reagent_index))

        return matching_reactions

    def run_all_reactions(self, node: TreeNode) -> list[TreeNode]:
        matching_reactions = self.get_reactions_for_fragments(fragments=node.fragments)

        product_nodes = []
        product_set = set()
        for reaction, fragment_to_reagent_index in matching_reactions:
            # Put fragments in the right order for the reaction
            fragments = sorted(node.fragments, key=lambda frag: fragment_to_reagent_index[frag])

            # Run reaction
            products = reaction.run_reactants(fragments)

            if len(products) == 0:
                raise ValueError('Reaction failed to produce products.')

            assert all(len(product) == 1 for product in products)

            # Convert product mols to SMILES (and remove Hs)
            products = [Chem.MolToSmiles(Chem.RemoveHs(product[0])) for product in products]

            # Filter out products that have already been created and deduplicate
            products = list(dict.fromkeys(product for product in products if product not in product_set))
            product_set |= set(products)

            # Create reaction log
            reaction_log = {
                'reaction_id': reaction.id,
                'reagent_ids': tuple(self.fragment_smiles_to_id.get(fragment, -1) for fragment in fragments),
            }

            product_nodes += [
                self.TreeNodeClass(
                    fragments=(product,),
                    unique_reagents=node.unique_reagents,
                    construction_log=node.construction_log + (reaction_log,),
                    rollout_num=self.rollout_num
                )
                for product in products
            ]

        return product_nodes

    def expand_node(self, node: TreeNode) -> list[TreeNode]:
        # Run all valid reactions on the current fragments to generate new molecules
        new_nodes = self.run_all_reactions(node=node)

        # Add all possible next fragments to the current fragment
        if node.num_fragments == 0:
            next_fragments = self.all_fragments
        else:
            next_fragments = self.get_next_fragments(fragments=node.fragments)

        # Limit the number of next fragments to expand
        if self.num_expand_nodes is not None and len(next_fragments) > self.num_expand_nodes:
            next_fragments = self.random_choice(next_fragments, size=self.num_expand_nodes, replace=False)

        # Convert next fragment tuples to TreeNodes
        new_nodes += [
            self.TreeNodeClass(
                fragments=node.fragments + (next_fragment,),
                unique_reagents=node.unique_reagents | {next_fragment},
                construction_log=node.construction_log,
                rollout_num=self.rollout_num
            )
            for next_fragment in next_fragments
        ]

        # Remove duplicates
        new_nodes = list(dict.fromkeys(new_nodes))

        return new_nodes

    def compute_mcts_score(self, node: TreeNode, total_visit_count: int) -> float:
        """Computes the MCTS score of a TreeNode."""
        mcts_score = node.Q() + node.U(n=total_visit_count)

        # Reduce MCTS score when node includes a common fragment
        if self.fragment_diversity and node.num_fragments > 0:
            max_reagent_count = max(self.reagent_counts[reagent_id] for reagent_id in node.unique_reagents)
            mcts_score /= np.exp((max_reagent_count - 1) / 100)  # -1 b/c every fragment appears once as its own node

        return mcts_score

    def rollout(self, node: TreeNode) -> float:
        """Performs a Monte Carlo Tree Search rollout.

        :param node: An TreeNode representing the root of the MCTS search.
        :return: The value (reward) of the rollout.
        """
        # Debugging
        if self.debug:
            print(f'Fragments = {node.fragments}')
            print(f'Num fragments = {node.num_fragments}')
            print(f'Num unique reagents = {len(node.unique_reagents)}')
            print(f'Num reactions = {node.num_reactions}')
            print(f'Score = {node.P}')
            print()

        # Stop the search if we've reached the maximum number of reactions
        if node.num_reactions >= self.max_reactions:
            return node.P

        # Expand the node both by running reactions with the current fragments and adding new fragments
        new_nodes = self.expand_node(node=node)

        if len(new_nodes) == 0:
            if node.num_fragments == 1:
                return node.P
            else:
                raise ValueError('Failed to expand a partially expanded node.')

        # Check the state map and merge with an existing node if available
        new_nodes = [self.state_map.get(new_node, new_node) for new_node in new_nodes]

        # If counting nodes, add node fragments to set
        if self.node_fragments_set is not None:
            self.node_fragments_set |= {node.fragments for node in new_nodes}

        # Add nodes with complete molecules to the state map
        for new_node in new_nodes:
            if new_node.num_fragments == 1 and new_node not in self.state_map:
                new_node.node_id = len(self.state_map)
                self.state_map[new_node] = new_node
                self.reagent_counts.update(new_node.unique_reagents)

        # Select a node based on the search type
        if self.search_type == 'random':
            selected_node = self.random_choice(new_nodes)

        elif self.search_type == 'greedy':
            # Randomly sample from top k
            # TODO: better method?
            top_k = 100
            sorted_nodes = sorted(new_nodes, key=lambda child: child.P, reverse=True)
            selected_node = self.random_choice(sorted_nodes[:top_k])

        elif self.search_type == 'mcts':
            total_visit_count = sum(new_node.N for new_node in new_nodes)
            # TODO: incorporate fragment counts into greedy search?
            selected_node = max(
                new_nodes,
                key=partial(self.compute_mcts_score, total_visit_count=total_visit_count)
            )

        else:
            raise ValueError(f'Search type "{self.search_type}" is not supported.')

        # Check the state map and merge with an existing node if available
        selected_node = self.state_map.setdefault(selected_node, selected_node)

        # Assign node ID as order in which the node was added
        if selected_node.node_id is None:
            selected_node.node_id = len(self.state_map)

        # Unroll the selected node
        v = self.rollout(node=selected_node)
        selected_node.W += v
        selected_node.N += 1

        return v

    def search(self) -> list[TreeNode]:
        """Runs the tree search, returning a list of TreeNode objects sorted from highest to lowest reward.

        NOTE: Only returns nodes with exactly one fragment, i.e., complete molecules, and at least one reaction.

        :return: A list of TreeNode objects sorted from highest to lowest reward.
        """
        for rollout_num in trange(self.n_rollout):
            self.rollout_num = rollout_num + 1

            if self.debug:
                print(f'Rollout {self.rollout_num}')

            self.rollout(node=self.root)

        # Get all the nodes representing fully constructed molecules that are not initial building blocks
        nodes = [node for _, node in self.state_map.items() if node.num_fragments == 1 and node.num_reactions > 0]

        # Sort by highest score and break ties by using node ID
        nodes = sorted(nodes, key=lambda node: (node.P, -node.node_id), reverse=True)

        return nodes

    @property
    def num_nodes(self) -> Optional[int]:
        """Returns the number of nodes in the search tree. Only available if count_nodes is True."""
        return len(self.node_fragments_set) if self.node_fragments_set is not None else None


def save_molecules(
        nodes: list[TreeNode],
        fragment_id_to_smiles: dict[int, str],
        save_path: Path
) -> None:
    # Convert construction logs from lists to dictionaries
    construction_dicts = []
    max_reaction_num = 0
    reaction_num_to_max_reagent_num = {}

    for node in nodes:
        construction_dict = {'num_reactions': len(node.construction_log)}
        max_reaction_num = max(max_reaction_num, len(node.construction_log))

        for reaction_index, reaction_log in enumerate(node.construction_log):
            reaction_num = reaction_index + 1
            construction_dict[f'reaction_{reaction_num}_id'] = reaction_log['reaction_id']

            reaction_num_to_max_reagent_num[reaction_num] = max(
                reaction_num_to_max_reagent_num.get(reaction_num, 0),
                len(reaction_log['reagent_ids'])
            )

            for reagent_index, reagent_id in enumerate(reaction_log['reagent_ids']):
                reagent_num = reagent_index + 1
                construction_dict[f'building_block_{reaction_num}_{reagent_num}_id'] = reagent_id
                construction_dict[f'building_block_{reaction_num}_{reagent_num}_smiles'] = fragment_id_to_smiles[reagent_id]

        construction_dicts.append(construction_dict)

    # Specify column order for CSV file
    columns = ['smiles', 'node_id', 'num_expansions', 'rollout_num', 'score', 'Q_value', 'num_reactions']

    for reaction_num in range(1, max_reaction_num + 1):
        columns.append(f'reaction_{reaction_num}_id')

        for reagent_num in range(1, reaction_num_to_max_reagent_num[reaction_num] + 1):
            columns.append(f'building_block_{reaction_num}_{reagent_num}_id')
            columns.append(f'building_block_{reaction_num}_{reagent_num}_smiles')

    # Save data
    # TODO: if incorporating train similarity, set it as a node property and then extract it here
    data = pd.DataFrame(
        data=[
            {
                'smiles': node.fragments[0],
                'node_id': node.node_id,
                'num_expansions': node.N,
                'rollout_num': node.rollout_num,
                'score': node.P,
                'Q_value': node.Q(),
                **construction_dict
            }
            for node, construction_dict in zip(nodes, construction_dicts)
        ],
        columns=columns
    )
    data.to_csv(save_path, index=False)


def create_model_scoring_fn(model_path: Path,
                            model_type: str,
                            fingerprint_type: Optional[str],
                            fragment_to_model_score: Optional[dict[str, float]],
                            binarize_scoring: float) -> Callable[[str], float]:
    """Creates a function that scores a molecule using a model."""
    # Check compatibility of model and fingerprint type
    if model_type != 'chemprop' and fingerprint_type is None:
        raise ValueError('Must define fingerprint_type if using sklearn model.')

    # Get model paths
    if model_path.is_dir():
        model_paths = list(model_path.glob('*.pt' if model_type == 'chemprop' else '*.pkl'))

        if len(model_paths) == 0:
            raise ValueError(f'Could not find any models in directory {model_path}.')
    else:
        model_paths = [model_path]

    # Load model and set up scoring function
    if model_type == 'chemprop':
        # Ensure reproducibility
        torch.manual_seed(0)
        torch.use_deterministic_algorithms(True)

        models = [load_checkpoint(path=model_path).eval() for model_path in model_paths]

        # Set up model scoring function for ensemble of chemprop models
        def model_scorer(smiles: str, fingerprint: Optional[np.ndarray]) -> float:
            fingerprint = [fingerprint] if fingerprint is not None else None

            return float(np.mean([
                model(batch=[[smiles]], features_batch=fingerprint).item() for model in models
            ]))
    else:
        models = []
        for model_path in model_paths:
            with open(model_path, 'rb') as f:
                models.append(pickle.load(f))

        # Set up model scoring function for ensemble of random forest or MLP models
        def model_scorer(smiles: str, fingerprint: np.ndarray) -> float:
            return float(np.mean([
                sklearn_predict(model=model, fingerprints=fingerprint.reshape(1, -1))[0] for model in models
            ]))

    @cache
    def model_scoring_fn(smiles: str) -> float:
        if fragment_to_model_score is None or smiles not in fragment_to_model_score:
            if fingerprint_type is not None:
                fingerprint = compute_fingerprint(smiles, fingerprint_type=fingerprint_type)
            else:
                fingerprint = None

            model_score = model_scorer(smiles=smiles, fingerprint=fingerprint)
        else:
            model_score = fragment_to_model_score[smiles]

        if binarize_scoring > 0:
            model_score = int(model_score >= binarize_scoring)

        return model_score

    return model_scoring_fn


def load_and_set_allowed_reaction_smiles(reaction_to_reagents_path: Path,
                                         fragment_id_to_smiles: dict[int, str]) -> None:
    """Sets the allowed reaction SMILES for each reaction in REAL_REACTIONS.

    :param reaction_to_reagents_path: Path to JSON file containing mapping from REAL reactions to allowed reagents.
    :param fragment_id_to_smiles: Dictionary mapping from fragment ID to fragment SMILES
    """
    # Load mapping from reactions to allowed reagent IDs
    with open(reaction_to_reagents_path) as f:
        # Dictionary mapping from reaction number to reagent number to reaction type (M or S) to set of reagent IDs
        reaction_to_reagent_ids: dict[int, dict[int, dict[str, set[int]]]] = json.load(f)

    # Map from reactions to reagent SMILES, merging reaction types
    reaction_to_reagent_smiles: dict[int, dict[int, set[str]]] = {
        int(reaction): {
            int(reagent): {
                fragment_id_to_smiles[reagent_id]
                for reagent_id_set in reaction_type_to_reagent_ids.values()
                for reagent_id in reagent_id_set
                if reagent_id in fragment_id_to_smiles
            }
            for reagent, reaction_type_to_reagent_ids in reagent_to_reaction_type.items()
        }
        for reaction, reagent_to_reaction_type in reaction_to_reagent_ids.items()
    }

    # Set the allowed SMILES for each reagent in each reaction
    set_allowed_reaction_smiles(
        reactions=REAL_REACTIONS,
        reaction_to_reagents=reaction_to_reagent_smiles
    )


def run_tree_search(args: Args) -> None:
    """Generate molecules combinatorially by performing a tree search."""
    # Validate arguments
    model_based_search = args.search_type in {'mcts', 'greedy'}

    if model_based_search:
        if args.model_type is None:
            raise ValueError('Must specify model_type when using model-based search.')

        if args.model_path is None:
            raise ValueError('Must specify model_path when using model-based search.')

        if args.fragment_to_model_score_path is None:
            print('For faster searching, please provide a fragment_to_model_score_path.')
    else:
        if args.model_type is not None:
            raise ValueError('Do not include model_type when using non-model-based search.')

        if args.model_path is not None:
            raise ValueError('Do not include model_path when using non-model-based search.')

        if args.train_similarity:
            raise ValueError('Do not use train similarity when using non-model-based search.')

        if args.noise:
            raise ValueError('Do not use noise when using non-model-based search.')

    if args.max_reactions > 1:
        raise NotImplementedError('Multiple reactions not yet implemented when using reaction_to_reagents.')

    # Create save directory and save arguments
    args.save_dir.mkdir(parents=True, exist_ok=True)
    args.save(args.save_dir / 'args.json')

    # Decide which reaction set to use
    if args.synnet_rxn:
        reactions = REAL_REACTIONS + SYNNET_REACTIONS
    else:
        reactions = REAL_REACTIONS

    # Load fragments
    # TODO: consider just using unique SMILES and ignoring IDs b/c multiple IDs per SMILES
    fragment_data = pd.read_csv(args.fragment_path)
    fragment_data.drop_duplicates(subset=args.smiles_column, inplace=True)

    # Map fragments to indices
    fragment_smiles_to_id = dict(zip(fragment_data[args.smiles_column], fragment_data[args.fragment_id_column]))
    fragment_id_to_smiles = dict(zip(fragment_data[args.fragment_id_column], fragment_data[args.smiles_column]))
    fragment_set = set(fragment_smiles_to_id)

    # Ensure unique fragment IDs
    if len(fragment_set) != len(fragment_smiles_to_id.values()):
        raise ValueError('Fragment IDs are not unique.')

    # Set the allowed reaction SMILES
    load_and_set_allowed_reaction_smiles(
        reaction_to_reagents_path=args.reaction_to_reagents_path,
        fragment_id_to_smiles=fragment_id_to_smiles
    )

    # Load mapping from SMILES to model scores
    if args.fragment_to_model_score_path is not None:
        with open(args.fragment_to_model_score_path) as f:
            fragment_to_model_score: dict[str, float] = json.load(f)

        if set(fragment_to_model_score) != fragment_set:
            raise ValueError('The fragments in fragment_to_model do not match the fragment set.')
    else:
        fragment_to_model_score = None

    # Define model scoring function
    if model_based_search:
        model_scoring_fn = create_model_scoring_fn(
            model_path=args.model_path,
            model_type=args.model_type,
            fingerprint_type=args.fingerprint_type,
            fragment_to_model_score=fragment_to_model_score,
            binarize_scoring=args.binarize_scoring
        )
    else:
        def model_scoring_fn(smiles: str) -> float:
            return 0.0

    if args.train_similarity:
        # Load mapping from SMILES to train similarity
        with open(args.fragment_to_train_similarity_path) as f:
            fragment_to_train_similarity: dict[str, float] = json.load(f)

        if set(fragment_to_train_similarity) != fragment_set:
            raise ValueError('The fragments in fragment_to_train_similarity do not match the fragment set.')

        # Load mapping from SMILES to train similarity
        with open(args.fragment_to_train_hits_similarity_path) as f:
            fragment_to_train_hits_similarity: dict[str, float] = json.load(f)

        if set(fragment_to_train_hits_similarity) != fragment_set:
            raise ValueError('The fragments in fragment_to_train_hits_similarity do not match the fragment set.')

        # Load train and compute Morgan fingerprints
        train = pd.read_csv(args.train_path)[args.smiles_column]
        train_morgans = compute_fingerprints(train, fingerprint_type='morgan')
        top_k = 10
        assert top_k <= len(train_morgans)

        # Load train hits and compute Morgan fingerprints
        train_hits = pd.read_csv(args.train_hits_path)[args.smiles_column]
        train_hits_morgans = compute_fingerprints(train_hits, fingerprint_type='morgan')
        train_hits_tversky = (train_hits_morgans.transpose() / train_hits_morgans.sum(axis=1))

        # Define train similarity scoring function (top 10 Tanimoto similarity)
        def train_similarity_scoring_fn(morgan_fingerprint: np.ndarray) -> float:
            train_distances = pairwise_distances(
                morgan_fingerprint.reshape(1, -1),
                train_morgans,
                metric='jaccard',
                n_jobs=-1
            )
            train_similarities = 1 - train_distances[0]
            argsort = np.argsort(train_similarities)
            top_k_index = argsort[-top_k]
            train_similarity_score = train_similarities[top_k_index]

            return train_similarity_score

        # TODO: maybe make this train similarity score exponential b/c
        # TODO: we only want to penalize significantly when the similarity is very high
        # TODO: but otherwise we don't care too much (i.e., not a linear relationship)
        # Define train hits similarity scoring function (top 1 Tversky similarity)
        def train_hits_similarity_scoring_fn(morgan_fingerprint: np.ndarray) -> float:
            train_hits_similarity_score = np.max(morgan_fingerprint @ train_hits_tversky)

            return train_hits_similarity_score

        # Define train scoring function including all comparisons to the training set
        @cache
        def train_scoring_fn(smiles: str) -> float:
            if smiles in fragment_set:
                train_similarity_score = fragment_to_train_similarity[smiles]
                train_hits_similarity_score = fragment_to_train_hits_similarity[smiles]
            else:
                morgan_fingerprint = compute_fingerprint(smiles, fingerprint_type='morgan')
                train_similarity_score = train_similarity_scoring_fn(morgan_fingerprint)
                train_hits_similarity_score = train_hits_similarity_scoring_fn(morgan_fingerprint)

            # Compute train score as a combination of the train similarity and train hits similarity
            # Strategy: encourage similarity to training set overall (top 10) but penalize similarity to nearest train hit
            # + 1 to ensure non-negative score since train hits similarity can be at most 1
            train_score = train_similarity_score - train_hits_similarity_score + 1

            return train_score

        # Define scoring function
        @cache
        def scoring_fn(smiles: str) -> float:
            return max(0, model_scoring_fn(smiles) + train_scoring_fn(smiles))
    else:
        # Define scoring function
        scoring_fn = model_scoring_fn

    # Define noisy scoring function
    rng = np.random.default_rng(seed=args.rng_seed)

    def noisy_scoring_fn(smiles: str) -> float:
        return scoring_fn(smiles) + rng.uniform(0, args.noise_limit)

    # Set up TreeSearcher
    tree_searcher = TreeSearcher(
        search_type=args.search_type,
        fragment_smiles_to_id=fragment_smiles_to_id,
        max_reactions=args.max_reactions,
        scoring_fn=noisy_scoring_fn if args.noise else scoring_fn,
        n_rollout=args.n_rollout,
        c_puct=args.c_puct,
        num_expand_nodes=args.num_expand_nodes,
        reactions=reactions,
        rng_seed=args.rng_seed,
        fragment_diversity=args.fragment_diversity,
        debug=args.debug,
        count_nodes=args.count_nodes
    )

    # Search for molecules
    start_time = datetime.now()
    nodes = tree_searcher.search()

    # Compute, print, and save stats
    stats = {
        'mcts_time': datetime.now() - start_time,
        'num_nodes_searched': tree_searcher.num_nodes,
        'num_nonzero_reaction_molecules': len(nodes)
    }

    print(f'MCTS time = {stats["mcts_time"]}')
    print(f'Number of full molecule, nonzero reaction nodes = {stats["num_nonzero_reaction_molecules"]:,}')

    if stats['num_nodes_searched'] is not None:
        print(f'Total number of nodes searched = {stats["num_nodes_searched"]:,}')

    pd.DataFrame(data=[stats]).to_csv(args.save_dir / 'mcts_stats.csv', index=False)

    # Save generated molecules
    save_molecules(
        nodes=nodes,
        fragment_id_to_smiles=fragment_id_to_smiles,
        save_path=args.save_dir / 'molecules.csv'
    )


if __name__ == '__main__':
    run_tree_search(Args().parse_args())
