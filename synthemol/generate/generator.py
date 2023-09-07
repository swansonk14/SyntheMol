"""Contains the Generator class, which generates molecules."""
import itertools
import time
from collections import Counter
from functools import partial
from typing import Callable, Literal

import numpy as np
import wandb
from chemfunc import get_fingerprint_generator
from rdkit import Chem
from scipy.special import softmax
from sklearn.metrics import pairwise_distances
from tqdm import trange

from synthemol.constants import OPTIMIZATION_TYPES
from synthemol.generate.node import Node
from synthemol.models import RLModel
from synthemol.reactions import Reaction
from synthemol.utils import random_choice


morgan_fingerprint_generator = get_fingerprint_generator('morgan')


class Generator:
    """A class that generates molecules."""

    def __init__(
            self,
            search_type: Literal['mcts', 'rl'],
            building_block_smiles_to_id: dict[str, int],
            max_reactions: int,
            scoring_fn: Callable[[str], float],
            explore_weight: float,
            num_expand_nodes: int | None,
            rl_temperature: float,
            rl_train_frequency: int,
            optimization: OPTIMIZATION_TYPES,
            reactions: tuple[Reaction],
            rng_seed: int,
            no_building_block_diversity: bool,
            store_nodes: bool,
            verbose: bool,
            rl_model: RLModel | None = None,
            replicate: bool = False,
            wandb_log: bool = False
    ) -> None:
        """Creates the Generator.

        :param search_type: The type of search to perform. 'mcts' = Monte Carlo tree search. 'rl' = Reinforcement learning.
        :param building_block_smiles_to_id: A dictionary mapping building block SMILES to their IDs.
        :param max_reactions: The maximum number of reactions to use to construct a molecule.
        :param scoring_fn: A function that takes as input a SMILES representing a molecule and returns a score.
        :param explore_weight: The hyperparameter that encourages exploration.
        :param num_expand_nodes: The number of tree nodes to expand when extending the child nodes in the search tree.
                                  If None, then all nodes are expanded.
        :param rl_temperature: The temperature parameter for the softmax function used to select building blocks.
                               Higher temperature means more exploration.
        :param rl_train_frequency: The number of rollouts between each training step of the RL model.
        :param optimization: Whether to maximize or minimize the score.
        :param reactions: A tuple of reactions that combine molecular building blocks.
        :param rng_seed: Seed for the random number generator.
        :param no_building_block_diversity: Whether to turn off the score modification that encourages diverse building blocks.
        :param store_nodes: Whether to store the child nodes of each node in the search tree.
                            This doubles the speed of the search but significantly increases
                            the memory usage (e.g., 450GB for 20,000 rollouts instead of 600 MB).
        :param verbose: Whether to print out additional statements during generation.
        :param rl_model: The RL model to use for the RL search. Must be provided if and only if search_type is 'rl'.
        :param replicate: This is necessary to replicate the results from the paper, but otherwise should not be used
                          since it limits the potential choices of building blocks.
        :param wandb_log: Whether to log results to Weights & Biases.
        """
        self.search_type = search_type
        self.building_block_smiles_to_id = building_block_smiles_to_id
        self.max_reactions = max_reactions
        self.scoring_fn = scoring_fn
        self.explore_weight = explore_weight
        self.num_expand_nodes = num_expand_nodes
        self.rl_temperature = rl_temperature
        self.rl_train_frequency = rl_train_frequency
        self.optimization = optimization
        self.reactions = reactions
        self.rng = np.random.default_rng(seed=rng_seed)
        self.building_block_diversity = not no_building_block_diversity
        self.store_nodes = store_nodes
        self.verbose = verbose
        self.replicate = replicate
        self.rl_model = rl_model
        self.wandb_log = wandb_log

        # Check that the search type is valid
        if (self.search_type == 'rl') != (self.rl_model is not None):
            raise ValueError('RL model must be provided if and only if search_type is "rl"')

        # Get all building blocks that are used in at least one reaction
        if self.replicate:
            self.all_building_blocks = list(dict.fromkeys(
                building_block
                for reaction in self.reactions
                for reactant in reaction.reactants
                for building_block in reactant.allowed_building_blocks
            ))
        else:
            self.all_building_blocks = sorted(
                building_block
                for building_block in self.building_block_smiles_to_id
                if any(reactant.has_match(building_block) for reaction in reactions for reactant in reaction.reactants)
            )

        # Get the function to use for optimization
        if self.optimization == 'maximize':
            self.optimization_fn = max
            self.optimization_sign = 1
            self.ascending_scores = False
        elif self.optimization == 'minimize':
            self.optimization_fn = min
            self.optimization_sign = -1
            self.ascending_scores = True
        else:
            raise ValueError(f'Invalid optimization type: {self.optimization}')

        # Initialize the root node
        self.rollout_num = 0
        self.root = Node(
            explore_weight=explore_weight,
            scoring_fn=scoring_fn,
            node_id=0,
            rollout_num=self.rollout_num
        )

        # Initialize the rollout num, node map, building block counts, and node to children
        self.node_map: dict[Node, Node] = {self.root: self.root}
        self.building_block_counts = Counter()
        self.node_to_children: dict[Node, list[Node]] = {}

        # Initialize the array of Morgan fingerprints of full molecules for diversity calculations
        self.full_molecule_morgan_fingerprints: np.ndarray | None = None

        # Initialize dictionary mapping tuple of molecules to cached RL score
        self.molecules_to_rl_score = {}

    def get_next_building_blocks(self, molecules: tuple[str]) -> list[str]:
        """Get the next building blocks that can be added to the given molecules.

        :param molecules: A tuple of SMILES strings representing the molecules to get the next building blocks for.
        :return: A list of SMILES strings representing the next building blocks that can be added to the given molecules.
        """
        # Initialize list of allowed building blocks
        available_building_blocks = []

        # Loop through each reaction
        for reaction in self.reactions:
            # Get indices of the reactants in this reaction
            reactant_indices = set(range(reaction.num_reactants))

            # Skip reaction if there's no room to add more reactants
            if len(molecules) >= reaction.num_reactants:
                continue

            # For each molecule, get a list of indices of reactants it matches
            reactant_matches_per_molecule = [
                reaction.get_reactant_matches(smiles=molecule)
                for molecule in molecules
            ]

            # Loop through products of reactant indices that the molecules match to
            # and for each product, if it matches to all separate reactants,
            # then include the missing reactants in the set of unfilled reactants
            for matched_reactant_indices in itertools.product(*reactant_matches_per_molecule):
                matched_reactant_indices = set(matched_reactant_indices)

                if len(matched_reactant_indices) == len(molecules):
                    for index in sorted(reactant_indices - matched_reactant_indices):
                        available_building_blocks += reaction.reactants[index].allowed_building_blocks

        # Remove duplicates but maintain order for reproducibility
        available_building_blocks = list(dict.fromkeys(available_building_blocks))

        return available_building_blocks

    def get_reactions_for_molecules(self, molecules: tuple[str]) -> list[tuple[Reaction, dict[str, int]]]:
        """Get all reactions that can be run on the given molecules.

        :param molecules: A tuple of SMILES strings representing the molecules to run reactions on.
        :return: A list of tuples, where each tuple contains a reaction and a dictionary mapping
                 the molecules to the indices of the reactants they match.
        """
        matching_reactions = []

        # Check each reaction to see if it can be run on the given molecules
        for reaction in self.reactions:
            # Skip reaction if the number of molecules doesn't match the number of reactants
            if len(molecules) != reaction.num_reactants:
                continue

            # For each molecule, get a list of indices of reactants it matches
            reactant_matches_per_molecule = [
                reaction.get_reactant_matches(smiles=molecule)
                for molecule in molecules
            ]

            # Include every assignment of molecules to reactants that fills all the reactants
            for matched_reactant_indices in itertools.product(*reactant_matches_per_molecule):
                if len(set(matched_reactant_indices)) == reaction.num_reactants:
                    molecule_to_reactant_index = dict(zip(molecules, matched_reactant_indices))
                    matching_reactions.append((reaction, molecule_to_reactant_index))

        return matching_reactions

    def run_all_reactions(self, node: Node) -> list[Node]:
        """Run all possible reactions for the molecules in the Node and return the resulting product Nodes.

        :param node: A Node to run reactions for.
        :return: A list of Nodes for the products of the reactions.
        """
        # Get all reactions that are possible for the molecules in the Node
        matching_reactions = self.get_reactions_for_molecules(molecules=node.molecules)

        # Run all possible reactions and create Nodes for the products
        product_nodes = []
        product_set = set()
        for reaction, molecule_to_reactant_index in matching_reactions:
            # Put molecules in the right order for the reaction
            molecules = sorted(node.molecules, key=lambda frag: molecule_to_reactant_index[frag])

            # Run reaction
            products = reaction.run_reactants(molecules)

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
                'building_block_ids': tuple(
                    self.building_block_smiles_to_id.get(molecule, -1)
                    for molecule in molecules
                ),
            }

            product_nodes += [
                Node(
                    explore_weight=self.explore_weight,
                    scoring_fn=self.scoring_fn,
                    molecules=(product,),
                    unique_building_block_ids=node.unique_building_block_ids,
                    construction_log=node.construction_log + (reaction_log,),
                    rollout_num=self.rollout_num
                )
                for product in products
            ]

        return product_nodes

    def get_child_nodes(self, node: Node) -> list[Node]:
        """Get the child Nodes of a given Node.

        Child Nodes are created in two ways:
        1. By running all possible reactions on the molecules in the Node and creating Nodes with the products.
        2. By adding all possible next building blocks to the molecules in the Node and creating a new Node for each.

        :param node: The Node to get the child Nodes of.
        :return: A list of child Nodes of the given Node.
        """
        # Run all valid reactions on the current molecules to combine them into new molecules
        new_nodes = self.run_all_reactions(node=node)

        # Add all possible next building blocks to the current molecules in the Node
        if node.num_molecules == 0:
            next_building_blocks = self.all_building_blocks
        else:
            next_building_blocks = self.get_next_building_blocks(molecules=node.molecules)

        # Optionally, limit the number of next nodes
        if self.num_expand_nodes is not None and len(next_building_blocks) > self.num_expand_nodes:
            next_building_blocks = random_choice(
                rng=self.rng,
                array=next_building_blocks,
                size=self.num_expand_nodes,
                replace=False
            )

        # Convert next node molecule tuples into Node objects
        new_nodes += [
            Node(
                explore_weight=self.explore_weight,
                scoring_fn=self.scoring_fn,
                molecules=node.molecules + (next_building_block,),
                unique_building_block_ids=node.unique_building_block_ids | {next_building_block},
                construction_log=node.construction_log,
                rollout_num=self.rollout_num
            )
            for next_building_block in next_building_blocks
        ]

        # Remove duplicates but maintain order for reproducibility
        new_nodes = list(dict.fromkeys(new_nodes))

        return new_nodes

    def compute_mcts_score(self, node: Node, total_visit_count: int) -> float:
        """Computes the MCTS score of a Node.

        :param node: A Node.
        :param total_visit_count: The total number of visits to all nodes at the same level as this Node.
        :return: The MCTS score of the Node.
        """
        # Compute initial MCTS score
        mcts_score = node.Q() + node.U(n=total_visit_count)

        # Optionally encourage building block diversity by reducing MCTS score based on building block usage
        if self.building_block_diversity and node.num_molecules > 0:
            # Determine the maximum number of times any building block in the molecule has been used
            max_building_block_count = max(
                self.building_block_counts[building_block_id]
                for building_block_id in node.unique_building_block_ids
            )

            # Reduce MCTS score based on maximum building block usage
            # Note: count - 1 is used because every building block appears once as its own node
            mcts_score /= np.exp((max_building_block_count - 1) / 100)

        return mcts_score

    def rollout(self, node: Node) -> float:
        """Performs a generation rollout.

        :param node: A Node representing the root of the generation.
        :return: The value (reward) of the rollout.
        """
        if self.verbose:
            print(f'Node {node.node_id} (rollout {self.rollout_num})')
            print(f'Molecules = {node.molecules}')
            print(f'Num molecules = {node.num_molecules}')
            print(f'Num unique building blocks = {len(node.unique_building_block_ids)}')
            print(f'Num reactions = {node.num_reactions}')
            print(f'Score = {node.P}')
            print()

        # Stop the search if we've reached the maximum number of reactions
        if node.num_reactions >= self.max_reactions:
            return node.P

        # If this node has already been visited and the children have been stored, get its children from the dictionary
        if node in self.node_to_children:
            child_nodes = self.node_to_children[node]

        # Otherwise, expand the node to get its children
        else:
            # Expand the node both by running reactions with the current molecules and adding new building blocks
            child_nodes = self.get_child_nodes(node=node)

            # Check the node map and merge with an existing node if available
            child_nodes = [self.node_map.get(new_node, new_node) for new_node in child_nodes]

            # Add nodes with complete molecules to the node map
            for child_node in child_nodes:
                if child_node.num_molecules == 1 and child_node not in self.node_map:
                    child_node.node_id = len(self.node_map)
                    self.node_map[child_node] = child_node
                    self.building_block_counts.update(child_node.unique_building_block_ids)

            # Save the number of children in order to maintain a total node count
            node.num_children = len(child_nodes)

            # If storing nodes, store the children
            if self.store_nodes:
                self.node_to_children[node] = child_nodes

        # If no new nodes were generated, return the current node's value
        if len(child_nodes) == 0:
            if node.num_molecules == 1:
                return node.P
            else:
                raise ValueError('Failed to expand a partially expanded node.')

        # Select a node based on the search type
        if self.search_type == 'mcts':
            # Select node with the highest MCTS score
            total_visit_count = sum(child_node.N for child_node in child_nodes)
            selected_node = self.optimization_fn(
                child_nodes,
                key=partial(self.compute_mcts_score, total_visit_count=total_visit_count)
            )
        elif self.search_type == 'rl':
            # Compute RL scores for any nodes missing RL scores
            child_node_molecules = [
                child_node.molecules for child_node in child_nodes
            ]
            child_node_molecules_missing_scores = [
                molecules for molecules in child_node_molecules
                if molecules not in self.molecules_to_rl_score
            ]

            if len(child_node_molecules_missing_scores) > 0:
                rl_scores = self.rl_model.predict(child_node_molecules_missing_scores)

                for molecules, rl_score in zip(child_node_molecules_missing_scores, rl_scores):
                    self.molecules_to_rl_score[molecules] = rl_score

            # Convert RL scores to temperature-scaled probabilities
            child_node_scores = np.array([self.molecules_to_rl_score[molecules] for molecules in child_node_molecules])
            child_node_probs = softmax(self.optimization_sign * child_node_scores / self.rl_temperature)

            # Select node proportional to the temperature-scaled RL score
            selected_node = self.rng.choice(child_nodes, p=child_node_probs)
        else:
            raise ValueError(f'Invalid search type: {self.search_type}')

        # Check the node map and merge with an existing node if available
        if selected_node in self.node_map:
            selected_node = self.node_map[selected_node]
        # Otherwise, assign node ID and add to node map
        else:
            selected_node.node_id = len(self.node_map)
            self.node_map[selected_node] = selected_node

        # Unroll the selected node
        v = self.rollout(node=selected_node)

        # Get max whole molecule (non-building block) score across rollouts as feedback
        if selected_node.num_molecules == 1 and node.num_reactions > 0:
            v = max(v, selected_node.P)

        # Update exploit score and visit count
        selected_node.W += v
        selected_node.N += 1

        # Save example to RL buffer
        if self.rl_model is not None and node.num_molecules > 0:
            self.rl_model.buffer(molecules=node.molecules, reward=v)

        return v

    def update_similarity(self) -> float:
        """Computes the Tanimoto similarity of the current generation and updates for future comparisons.

        :return: The Tanimoto similarity of the current generation compared to previous generations.
        """
        # Get full molecule nodes from this generation
        new_full_molecule_nodes = self.get_full_molecule_nodes(
            rollout_start=self.rollout_num
        )

        # If there is no new full molecule, then a duplicate was found
        if len(new_full_molecule_nodes) == 0:
            return 1.0

        # Get full molecule SMILES from this generation
        new_full_molecule_smiles = [
            node.molecules[0]
            for node in new_full_molecule_nodes
        ]

        # Compute the Morgan fingerprints of the new nodes
        new_full_molecule_morgan_fingerprints = np.array([
            morgan_fingerprint_generator(smiles)
            for smiles in new_full_molecule_smiles
        ])

        # If this is the first generation, save the fingerprints and return
        if self.full_molecule_morgan_fingerprints is None:
            self.full_molecule_morgan_fingerprints = new_full_molecule_morgan_fingerprints
            return 0.0

        # Compute the Tanimoto similarity between the new molecules and the previous molecules
        tanimoto_distances = pairwise_distances(
            new_full_molecule_morgan_fingerprints,
            self.full_molecule_morgan_fingerprints,
            metric='jaccard',
            n_jobs=-1
        )
        tanimoto_similarities = 1 - tanimoto_distances

        # Compute average maximum similarity across new molecules
        avg_max_similarity = float(np.mean(np.max(tanimoto_similarities, axis=1)))

        # Update the full molecule fingerprints
        self.full_molecule_morgan_fingerprints = np.concatenate(
            [self.full_molecule_morgan_fingerprints, new_full_molecule_morgan_fingerprints],
            axis=0
        )

        return avg_max_similarity

    def get_full_molecule_nodes(
            self,
            rollout_start: int | None = None,
            rollout_end: int | None = None
    ) -> list[Node]:
        """Returns a list of all nodes with complete molecules.

        :param rollout_start: The first rollout to include. If None, starts from the first rollout.
        :param rollout_end: The last rollout to include. If None, ends at the last rollout.
        :return: A list of all nodes with complete molecules, sorted by score.
        """
        # Get all the Nodes representing fully constructed molecules
        nodes = [
            node
            for node in self.node_map
            if node.num_molecules == 1 and node.num_reactions > 0
        ]

        # Set rollout start and end
        if rollout_start is None:
            rollout_start = min(node.rollout_num for node in nodes)

        if rollout_end is None:
            rollout_end = max(node.rollout_num for node in nodes)

        # Get all the Nodes within the provided rollouts
        nodes = [
            node
            for node in nodes
            if rollout_start <= node.rollout_num <= rollout_end
        ]

        # Sort Nodes by score and break ties by using Node ID
        nodes = sorted(
            nodes,
            key=lambda node: (node.P, (1 if self.ascending_scores else -1) * node.node_id),
            reverse=not self.ascending_scores
        )

        return nodes

    def generate(self, n_rollout: int) -> list[Node]:
        """Generate molecules for the specified number of rollouts.

        NOTE: Only returns Nodes with exactly one molecule and at least one reaction.

        :param n_rollout: The number of rollouts to perform.
        :return: A list of Node objects sorted from best to worst score from these rollouts.
        """
        # Set up rollout bounds
        rollout_start = self.rollout_num + 1
        rollout_end = rollout_start + n_rollout

        # Run the generation algorithm for the specified number of rollouts
        for rollout_num in trange(rollout_start, rollout_end):
            # Record rollout number
            self.rollout_num = rollout_num

            # Set up rollout stats
            rollout_stats = {'Rollout Number': rollout_num}

            # Run rollout
            start_time = time.time()
            rollout_stats['Rollout Score'] = self.rollout(node=self.root)
            rollout_stats['Rollout Time'] = time.time() - start_time

            # Compute similarity of new molecules compared to previous molecules
            start_time = time.time()
            rollout_stats['Rollout Similarity'] = self.update_similarity()
            rollout_stats['Similarity Time'] = time.time() - start_time

            # Determine number of unique full molecules found
            rollout_stats['Unique Molecules'] = sum(
                node.num_molecules == 1 and node.num_reactions > 0
                for node in self.node_map
            )

            # Train and evaluate RL model
            if self.rl_model is not None and rollout_num % self.rl_train_frequency == 0:
                # Reset RL scores since RL model is being updated
                self.molecules_to_rl_score = {}

                # Train model
                start_time = time.time()
                self.rl_model.train()
                rollout_stats['RL Train Time'] = time.time() - start_time
                rollout_stats['RL Train Examples'] = self.rl_model.buffer_size

                # Evaluate model
                if self.wandb_log:
                    start_time = time.time()
                    rollout_stats |= self.rl_model.evaluate()
                    rollout_stats['RL Eval Time'] = time.time() - start_time

            # Log rollout stats
            if self.wandb_log:
                wandb.log(rollout_stats)

        # Get all the Nodes representing fully constructed molecules within these rollouts sorted by score
        nodes = self.get_full_molecule_nodes(
            rollout_start=rollout_start,
            rollout_end=rollout_end
        )

        return nodes

    @property
    def approx_num_nodes_searched(self) -> int:
        """Gets the approximate number of nodes seen during the search.

        Note: This will over count any node that appears as a child node of multiple parent nodes.
        """
        return 1 + sum(node.num_children for node in self.node_map)

    @property
    def num_nodes_searched(self) -> int:
        """Gets the precise number of nodes seen during the search. Only possible if store_nodes is True."""
        if not self.store_nodes:
            raise ValueError('Cannot get the precise number of nodes searched if store_nodes is False.'
                             'Use approx_num_nodes_searched instead.')

        # Get a set of all nodes and child nodes that have been visited
        visited_nodes = set()
        for node, children in self.node_to_children.items():
            visited_nodes.add(node)
            visited_nodes.update(children)

        return len(visited_nodes)
