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
from rdkit import Chem
from sklearn.metrics import pairwise_distances
from tap import Tap
from tqdm import trange

from chem_utils.molecular_fingerprints import compute_fingerprint, compute_fingerprints
from chemprop.utils import load_checkpoint
from real_reactions import Reaction, REAL_REACTIONS, SYNNET_REACTIONS


# TODO: make certain arguments optional if doing random search
class Args(Tap):
    model_path: Path
    """Path to a directory of model checkpoints or to a specific PKL or PT file containing a trained model."""
    fragment_path: Path
    """Path to CSV file containing molecular building blocks."""
    fragment_to_model_score_path: Path
    """Path to JSON file containing a dictionary mapping fragments to model prediction scores."""
    reagent_to_fragments_path: Path
    """Path to JSON file containing a dictionary mapping from reagents to fragments."""
    train_similarity: bool = False
    """Whether to incorporate similarity to the train set as part of the MCTS score."""
    fragment_to_train_similarity_path: Optional[Path] = None
    """Path to JSON file containing a dictionary mapping fragments to train similarities
     (top 10 nearest neighbor Tanimoto). Only used if train_similarity is True."""
    fragment_to_train_hits_similarity_path: Optional[Path] = None
    """Path to JSON file containing a dictionary mapping fragments to train hits similarities
     (nearest neighbor Tversky). Only used if train_similarity is True."""
    train_path: Optional[Path] = None
    """Path to CSV file containing SMILES for the molecules in the training set.
    Only used if train_similarity is True."""
    train_hits_path: Optional[Path] = None
    """Path to CSV file containing SMILES for the active molecules in the training set.
    Only used if train_similarity is True."""
    save_dir: Path
    """Path to directory where the generated molecules will be saved."""
    search_type: Literal['random', 'greedy', 'mcts']
    """Type of search to perform."""
    smiles_column: str = 'smiles'
    """Name of the column containing SMILES."""
    fragment_id_column: str = 'Reagent_ID'
    """Name of the column in fragment_path that contains fragment IDs."""
    max_reactions: int = 3
    """Maximum number of reactions that can be performed to expand fragments into molecules."""
    n_rollout: int = 100
    """The number of times to run the tree search."""
    c_puct: float = 10.0
    """The hyperparameter that encourages exploration."""
    num_expand_nodes: Optional[int] = None
    """The number of tree nodes to expand when extending the child nodes in the search tree."""
    model_type: Literal['rf', 'mlp', 'chemprop']
    """Type of model to train. 'rf' = random forest. 'mlp' = multilayer perceptron."""
    fingerprint_type: Literal['morgan', 'rdkit']
    """Type of fingerprints to use as input features."""
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


class TreeNode:
    """A node in a tree search representing a step in the molecule construction process."""

    def __init__(self,
                 c_puct: float,
                 scoring_fn: Callable[[str], float],
                 node_id: Optional[int] = None,
                 fragments: Optional[tuple[str]] = None,
                 reagent_counts: Optional[Counter] = None,
                 construction_log: Optional[tuple[dict[str, Any]]] = None,
                 rollout_num: Optional[int] = None) -> None:
        """Initializes the TreeNode object.

        :param c_puct: The hyperparameter that encourages exploration.
        :param scoring_fn: A function that takes as input a SMILES representing a molecule and returns a score.
        :param node_id: The ID of the node, which should correspond to the order in which the node was visited.
        :param fragments: A tuple of SMILES containing the fragments for the next reaction.
                         The first element is the currently constructed molecule while the remaining elements
                         are the fragments that are about to be added.
        :param reagent_counts: A Counter mapping reagent indices to the number of times they appear in this node.
        :param construction_log: A tuple of dictionaries containing information about each reaction.
        :param rollout_num: The number of the rollout on which this node was created.
        """
        self.c_puct = c_puct
        self.scoring_fn = scoring_fn
        self.node_id = node_id
        self.fragments = fragments if fragments is not None else tuple()
        self.reagent_counts = reagent_counts if reagent_counts is not None else Counter()
        self.construction_log = construction_log if construction_log is not None else tuple()
        # TODO: maybe change sum (really mean) to max since we don't care about finding the best leaf node, just the best node along the way?
        self.W = 0.0  # The sum of the leaf node values for leaf nodes that descend from this node.
        self.N = 0  # The number of times this node has been expanded.
        self.children: list[TreeNode] = []
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

    @cached_property
    def unique_reagents(self) -> set[int]:
        """A set of unique reagents in this node. (Note: The value is cached, so it assumes the node is immutable.)"""
        return set(self.reagent_counts)

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
                 fragment_to_id: dict[str, int],
                 reagent_to_fragments: dict[str, list[str]],
                 max_reactions: int,
                 scoring_fn: Callable[[str], float],
                 n_rollout: int,
                 c_puct: float,
                 num_expand_nodes: Optional[int],
                 reactions: list[Reaction],
                 rng_seed: int,
                 fragment_diversity: bool) -> None:
        """Creates the TreeSearcher object.

        :param search_type: Type of search to perform.
        :param fragment_to_id: A dictionary mapping fragment SMILES to their IDs.
        :param reagent_to_fragments: A dictionary mapping a reagent SMARTS to all fragment SMILES that match that reagent.
        :param max_reactions: The maximum number of reactions to use to construct a molecule.
        :param scoring_fn: A function that takes as input a SMILES representing a molecule and returns a score.
        :param n_rollout: The number of times to run the tree search.
        :param c_puct: The hyperparameter that encourages exploration.
        :param num_expand_nodes: The number of tree nodes to expand when extending the child nodes in the search tree.
        :param reactions: A list of chemical reactions that can be used to combine molecular fragments.
        :param rng_seed: Seed for the random number generator.
        :param fragment_diversity: Whether to encourage the use of diverse fragments by modifying the score.
        """
        self.search_type = search_type
        self.fragment_to_id = fragment_to_id
        self.reagent_to_fragments = reagent_to_fragments
        self.all_fragments = list(dict.fromkeys(
            fragment
            for reagent in reagent_to_fragments
            for fragment in self.reagent_to_fragments[reagent]
        ))
        self.max_reactions = max_reactions
        self.scoring_fn = scoring_fn
        self.n_rollout = n_rollout
        self.c_puct = c_puct
        self.num_expand_nodes = num_expand_nodes
        self.reactions = reactions
        self.rng = np.random.default_rng(seed=rng_seed)
        self.fragment_diversity = fragment_diversity

        self.rollout_num = 0
        self.TreeNodeClass = partial(TreeNode, c_puct=c_puct, scoring_fn=scoring_fn)
        self.root = self.TreeNodeClass(node_id=1, rollout_num=0)
        self.state_map: dict[TreeNode, TreeNode] = {self.root: self.root}
        self.reagent_counts = Counter()

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
                if reagent.has_substruct_match(fragment)
            ]
            for fragment in fragments
        ]

    # TODO: documentation
    def get_next_fragments(self, fragments: tuple[str]) -> list[str]:
        # For each reaction these fragments can participate in, get all the unfilled reagents
        unfilled_reagents = set()
        for reaction in self.reactions:
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
                    unfilled_reagents |= {reaction.reagents[index] for index in
                                          reagent_indices - matched_reagent_indices}

        if len(unfilled_reagents) == 0:
            return []

        # Get all the fragments that match the other reagents
        # TODO: cache this? or allow duplicates?
        available_fragments = list(dict.fromkeys(
            fragment
            for reagent in sorted(unfilled_reagents, key=lambda reagent: reagent.smarts)
            for fragment in self.reagent_to_fragments[reagent.smarts]
        ))

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
                'reagent_ids': [self.fragment_to_id.get(fragment, -1) for fragment in fragments],
                'reagent_smiles': [fragment for fragment in fragments]
            }

            product_nodes += [
                self.TreeNodeClass(
                    fragments=(product,),
                    reagent_counts=node.reagent_counts,
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
                reagent_counts=node.reagent_counts + Counter({self.fragment_to_id[next_fragment]: 1}),
                construction_log=node.construction_log,
                rollout_num=self.rollout_num
            )
            for next_fragment in next_fragments
        ]

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
        # Stop the search if we've reached the maximum number of reactions
        if node.num_reactions >= self.max_reactions:
            return node.P

        # TODO: prevent exploration of a subtree that has been fully explored
        # Expand if this node has never been visited
        if len(node.children) == 0:
            # Expand the node both by running reactions with the current fragments and adding new fragments
            new_nodes = self.expand_node(node=node)

            if len(new_nodes) == 0:
                if node.num_fragments == 1:
                    return node.P
                else:
                    raise ValueError('Failed to expand a partially expanded node.')

            # Add the expanded nodes as children
            child_set = set()
            for new_node in new_nodes:

                # Check whether this node was already added as a child
                if new_node not in child_set:

                    # Check the state map and merge with an existing node if available
                    new_node = self.state_map.setdefault(new_node, new_node)

                    # Assign node ID as order in which the node was added
                    if new_node.node_id is None:
                        new_node.node_id = len(self.state_map)

                        # Increment reagent counts if this node is a full molecule
                        if new_node.num_fragments == 1:
                            self.reagent_counts.update(new_node.unique_reagents)

                    # Add the new node as a child of the tree node
                    node.children.append(new_node)
                    child_set.add(new_node)

        # TODO: think more about balancing reaction nodes vs additional fragment nodes

        # Select a child node based on the search type
        if self.search_type == 'random':
            selected_node = self.random_choice(node.children)

        elif self.search_type == 'greedy':
            # Randomly sample from top k
            # TODO: better method?
            top_k = 100
            sorted_children = sorted(node.children, key=lambda child: child.P, reverse=True)
            selected_node = self.random_choice(sorted_children[:top_k])

        elif self.search_type == 'mcts':
            total_visit_count = sum(child.N for child in node.children)
            # TODO: incorporate fragment counts into greedy search?
            selected_node = max(
                node.children,
                key=partial(self.compute_mcts_score, total_visit_count=total_visit_count)
            )

        else:
            raise ValueError(f'Search type "{self.search_type}" is not supported.')

        # Unroll the child node
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
            self.rollout(node=self.root)

        # Get all the nodes representing fully constructed molecules that are not initial building blocks
        nodes = [node for _, node in self.state_map.items() if node.num_fragments == 1 and node.num_reactions > 0]

        # Sort by highest score and break ties by using node ID
        nodes = sorted(nodes, key=lambda node: (node.P, node.node_id), reverse=True)

        return nodes

    @property
    def num_nodes(self) -> int:
        """Gets the number of nodes in the search tree."""
        return len(self.state_map)


def save_molecules(nodes: list[TreeNode], save_path: Path) -> None:
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

            for reagent_index, (reagent_id, reagent_smiles) in enumerate(zip(reaction_log['reagent_ids'],
                                                                             reaction_log['reagent_smiles'])):
                reagent_num = reagent_index + 1
                construction_dict[f'reagent_{reaction_num}_{reagent_num}_id'] = reagent_id
                construction_dict[f'reagent_{reaction_num}_{reagent_num}_smiles'] = reagent_smiles

        construction_dicts.append(construction_dict)

    # Specify column order for CSV file
    columns = ['smiles', 'node_id', 'num_expansions', 'rollout_num', 'score', 'Q_value', 'num_reactions']

    for reaction_num in range(1, max_reaction_num + 1):
        columns.append(f'reaction_{reaction_num}_id')

        for reagent_num in range(1, reaction_num_to_max_reagent_num[reaction_num] + 1):
            columns.append(f'reagent_{reaction_num}_{reagent_num}_id')
            columns.append(f'reagent_{reaction_num}_{reagent_num}_smiles')

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
                            fingerprint_type: str,
                            fragment_to_model_score: dict[str, float],
                            binarize_scoring: float) -> Callable[[str], float]:
    """Creates a function that scores a molecule using a model."""
    # Get model paths
    if model_path.is_dir():
        model_paths = list(model_path.glob('*.pt' if model_type == 'chemprop' else '*.pkl'))

        if len(model_paths) == 0:
            raise ValueError(f'Could not find any models in directory {model_path}.')
    else:
        model_paths = [model_path]

    # Load model and set up scoring function
    if model_type == 'chemprop':
        models = [load_checkpoint(path=model_path).eval() for model_path in model_paths]

        def model_scorer(smiles: str, fingerprint: np.ndarray) -> float:
            return float(np.mean([model(batch=[[smiles]], features_batch=[fingerprint]).item() for model in models]))
    else:
        models = []
        for model_path in model_paths:
            with open(model_path, 'rb') as f:
                models.append(pickle.load(f))

        def model_scorer(smiles: str, fingerprint: np.ndarray) -> float:
            return float(np.mean([model.predict_proba([fingerprint])[0, 1] for model in models]))

    @cache
    def model_scoring_fn(smiles: str) -> float:
        if smiles not in fragment_to_model_score:
            fingerprint = compute_fingerprint(smiles, fingerprint_type=fingerprint_type)
            model_score = model_scorer(smiles=smiles, fingerprint=fingerprint)
        else:
            model_score = fragment_to_model_score[smiles]

        if binarize_scoring > 0:
            model_score = int(model_score >= binarize_scoring)

        return model_score

    return model_scoring_fn


def run_tree_search(args: Args) -> None:
    """Generate molecules combinatorially by performing a tree search."""
    # Create save directory and save arguments
    args.save_dir.mkdir(parents=True, exist_ok=True)
    args.save(args.save_dir / 'args.json')

    if args.synnet_rxn:
        reactions = REAL_REACTIONS + SYNNET_REACTIONS
    else:
        reactions = REAL_REACTIONS

    # Load fragments
    fragment_data = pd.read_csv(args.fragment_path)
    fragment_data.drop_duplicates(subset=[args.smiles_column], inplace=True)

    # Map fragments to indices
    fragment_to_id = dict(zip(fragment_data[args.smiles_column], fragment_data[args.fragment_id_column]))
    fragment_set = set(fragment_to_id)

    # Ensure unique fragment IDs
    if len(fragment_set) != len(fragment_to_id.values()):
        raise ValueError('Fragment IDs are not unique.')

    # Load mapping from SMILES to model scores
    with open(args.fragment_to_model_score_path) as f:
        fragment_to_model_score: dict[str, float] = json.load(f)

    if set(fragment_to_model_score) != fragment_set:
        raise ValueError('The fragments in fragment_to_model do not match the fragment set.')

    # Load mapping from reagents to fragments
    with open(args.reagent_to_fragments_path) as f:
        reagent_to_fragments: dict[str, list[str]] = json.load(f)

    if not ({fragment for fragments in reagent_to_fragments.values() for fragment in fragments} <= fragment_set):
        raise ValueError('The fragments in reagent_to_fragments is not a subset of the fragment set.')

    # Define model scoring function
    model_scoring_fn = create_model_scoring_fn(
        model_path=args.model_path,
        model_type=args.model_type,
        fingerprint_type=args.fingerprint_type,
        fragment_to_model_score=fragment_to_model_score,
        binarize_scoring=args.binarize_scoring
    )

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
        fragment_to_id=fragment_to_id,
        reagent_to_fragments=reagent_to_fragments,
        max_reactions=args.max_reactions,
        scoring_fn=noisy_scoring_fn if args.noise else scoring_fn,
        n_rollout=args.n_rollout,
        c_puct=args.c_puct,
        num_expand_nodes=args.num_expand_nodes,
        reactions=reactions,
        rng_seed=args.rng_seed,
        fragment_diversity=args.fragment_diversity
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
    print(f'Total number of nodes searched = {stats["num_nodes_searched"]:,}')
    print(f'Number of full molecule, nonzero reaction nodes = {stats["num_nonzero_reaction_molecules"]:,}')

    pd.DataFrame(data=[stats]).to_csv(args.save_dir / 'stats.csv', index=False)

    # Save generated molecules
    save_molecules(nodes=nodes, save_path=args.save_dir / 'molecules.csv')


if __name__ == '__main__':
    run_tree_search(Args().parse_args())
