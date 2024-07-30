"""Contains the Generator class, which generates molecules."""
import itertools
import pickle
import time
from collections import Counter
from functools import partial
from pathlib import Path
from typing import Any, Callable, Literal

import numpy as np
import wandb
from chemfunc import get_fingerprint_generator
from scipy.special import softmax
from sklearn.metrics import pairwise_distances
from tqdm import trange

from synthemol.generate.logs import ConstructionLog, ReactionLog
from synthemol.generate.score_weights import ScoreWeights
from synthemol.generate.node import Node
from synthemol.generate.scorer import MoleculeScorer
from synthemol.generate.rl_models import RLModel
from synthemol.reactions import Reaction
from synthemol.utils import random_choice


morgan_fingerprint_generator = get_fingerprint_generator("morgan")


class Generator:
    """A class that generates molecules."""

    def __init__(
        self,
        search_type: Literal["mcts", "rl"],
        chemical_space_to_building_block_smiles_to_id: dict[str, dict[str, str]],
        max_reactions: int,
        scorer: MoleculeScorer,
        score_weights: ScoreWeights,
        success_comparators: tuple[Callable[[float], bool], ...] | None,
        explore_weight: float,
        num_expand_nodes: int | None,
        rl_base_temperature: float,
        rl_temperature_similarity_target: float | None,
        rl_train_frequency: int,
        reactions: tuple[Reaction, ...],
        rng_seed: int,
        no_building_block_diversity: bool,
        store_nodes: bool,
        verbose: bool,
        rl_model: RLModel | None = None,
        rolling_average_weight: float = 0.98,
        min_temperature: float = 0.001,
        max_temperature: float = 10.0,
        min_score_weight: float = 0.001,
        replicate: bool = False,
        wandb_log: bool = False,
        log_path: Path | None = None,
    ) -> None:
        """Creates the Generator.

        :param search_type: The type of search to perform. 'mcts' = Monte Carlo tree search. 'rl' = Reinforcement learning.
        :param chemical_space_to_building_block_smiles_to_id: A dictionary mapping chemical space to building block
            SMILES to IDs.
        :param max_reactions: The maximum number of reactions to use to construct a molecule.
        :param scorer: A callable object that takes as input a SMILES representing a molecule and returns a score.
        :param score_weights: The weights to use for the scoring function.
        :param success_comparators: A tuple of functions that take as input a score and return whether the score
            indicates a successful molecule. If provided, then the score weights are dynamically set.
        :param explore_weight: The hyperparameter that encourages exploration.
        :param num_expand_nodes: The number of tree nodes to expand when extending the child nodes in the search tree.
            If None, then all nodes are expanded.
        :param rl_base_temperature: The initial temperature parameter for the softmax function used to select building blocks.
            Higher temperature means more exploration. If rl_temperature_similarity_target is provided,
            the temperature is adjusted based on generated molecule diversity.
        :param rl_temperature_similarity_target: If provided, adjusts the temperature to obtain the maximally scoring molecules
            that are at most this similar to previously generated molecules. Starts with
            the temperature provided by rl_base_temperature. If None, the temperature is not adjusted.
        :param rl_train_frequency: The number of rollouts between each training step of the RL model.
        :param reactions: A tuple of reactions that combine molecular building blocks.
        :param rng_seed: Seed for the random number generator.
        :param no_building_block_diversity: Whether to turn off the score modification that encourages diverse building blocks.
        :param store_nodes: Whether to store the child nodes of each node in the search tree.
            This doubles the speed of the search but significantly increases
            the memory usage (e.g., 450GB for 20,000 rollouts instead of 600 MB).
        :param verbose: Whether to print out additional statements during generation.
        :param rl_model: The RL model to use for the RL search. Must be provided if and only if search_type is 'rl'.
        :param rolling_average_weight: The weight to use for the rolling average of similarity (dynamic temperature)
            and success (dynamic model weights).
        :param min_temperature: The minimum temperature when using dynamic temperature.
        :param max_temperature: The maximum temperature when using dynamic temperature.
        :param min_score_weight: The minimum score weight (pre-normalization) when using dynamic score weights.
        :param replicate: This is necessary to replicate the results from the paper, but otherwise should not be used
            since it limits the potential choices of building blocks.
        :param wandb_log: Whether to log results to Weights & Biases.
        :param log_path: Path to a PKL file to save logs to. If None, logs are not saved locally.
        """
        self.search_type = search_type
        self.chemical_space_to_building_block_smiles_to_id = (
            chemical_space_to_building_block_smiles_to_id
        )
        self.max_reactions = max_reactions
        self.scorer = scorer
        self.score_weights = score_weights
        self.success_comparators = success_comparators
        self.explore_weight = explore_weight
        self.num_expand_nodes = num_expand_nodes
        self.rl_temperature = rl_base_temperature
        self.rl_temperature_similarity_target = rl_temperature_similarity_target
        self.rl_train_frequency = rl_train_frequency
        self.reactions = reactions
        self.rng = np.random.default_rng(seed=rng_seed)
        self.building_block_diversity = not no_building_block_diversity
        self.store_nodes = store_nodes
        self.verbose = verbose
        self.replicate = replicate
        self.rl_model = rl_model
        self.rolling_average_weight = rolling_average_weight
        self.min_temperature = min_temperature
        self.max_temperature = max_temperature
        self.min_score_weight = min_score_weight
        self.wandb_log = wandb_log
        self.log_path = log_path

        # Set up list of rollout stats
        self.rollout_stats_record: list[dict[str, Any]] = []

        # Check that the search type is valid
        if (self.search_type == "rl") != (self.rl_model is not None):
            raise ValueError(
                'RL model must be provided if and only if search_type is "rl"'
            )

        # Check that score weights is a list if success comparators is provided
        if self.success_comparators is not None and self.score_weights.immutable:
            raise ValueError(
                "Score weights must be mutable if success comparators are provided"
            )

        # Get all building blocks that are used in at least one reaction
        if self.replicate:
            self.all_building_blocks = list(
                dict.fromkeys(
                    building_block
                    for reaction in self.reactions
                    for reactant in reaction.reactants
                    for building_block in reactant.allowed_building_blocks
                )
            )
        else:
            self.all_building_blocks = sorted(
                building_block
                for building_block_smiles_to_id in self.chemical_space_to_building_block_smiles_to_id.values()
                for building_block in building_block_smiles_to_id
                if any(
                    reactant.has_match(building_block)
                    for reaction in reactions
                    for reactant in reaction.reactants
                )
            )

        # Initialize the root node
        self.rollout_num = 0
        self.root = Node(
            explore_weight=explore_weight,
            scorer=scorer,
            node_id=0,
            rollout_num=self.rollout_num,
        )

        # Initialize the rollout num, node map, building block counts, and node to children
        self.node_map: dict[Node, Node] = {self.root: self.root}
        self.building_block_counts = Counter()
        self.node_to_children: dict[Node, list[Node]] = {}

        # Initialize the array of Morgan fingerprints of full molecules for diversity calculations
        self.full_molecule_morgan_fingerprints: np.ndarray | None = None

        # Initialize dictionary mapping tuple of molecules to cached RL score
        self.molecules_to_rl_score = {}

        # Set up rolling average similarity for temperature adjustment
        if self.rl_temperature_similarity_target is not None:
            self.rolling_average_similarity = self.rl_temperature_similarity_target
        else:
            self.rolling_average_similarity = 0.0

        # Set up rolling average successes for score weight adjustment
        self.rolling_average_success_rate = np.zeros(self.score_weights.num_weights)

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
                reaction.get_reactant_matches(smiles=molecule) for molecule in molecules
            ]

            # Loop through products of reactant indices that the molecules match to
            # and for each product, if it matches to all separate reactants,
            # then include the missing reactants in the set of unfilled reactants
            for matched_reactant_indices in itertools.product(
                *reactant_matches_per_molecule
            ):
                matched_reactant_indices = set(matched_reactant_indices)

                if len(matched_reactant_indices) == len(molecules):
                    for index in sorted(reactant_indices - matched_reactant_indices):
                        available_building_blocks += reaction.reactants[
                            index
                        ].allowed_building_blocks

        # Remove duplicates but maintain order for reproducibility
        available_building_blocks = list(dict.fromkeys(available_building_blocks))

        return available_building_blocks

    def get_reactions_for_molecules(
        self, molecules: tuple[str]
    ) -> list[tuple[Reaction, dict[str, int]]]:
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
                reaction.get_reactant_matches(smiles=molecule) for molecule in molecules
            ]

            # Include every assignment of molecules to reactants that fills all the reactants
            for matched_reactant_indices in itertools.product(
                *reactant_matches_per_molecule
            ):
                if len(set(matched_reactant_indices)) == reaction.num_reactants:
                    molecule_to_reactant_index = dict(
                        zip(molecules, matched_reactant_indices)
                    )
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
            molecules = sorted(
                node.molecules, key=lambda frag: molecule_to_reactant_index[frag]
            )

            # Run reaction
            products = reaction.run_reactants(molecules)

            # Filter out products that have already been created and deduplicate
            products = list(
                dict.fromkeys(
                    product for product in products if product not in product_set
                )
            )
            product_set |= set(products)

            # Create reaction log
            reaction_log = ReactionLog(
                chemical_space=reaction.chemical_space,
                reaction_id=reaction.reaction_id,
                reactant_ids=tuple(
                    self.chemical_space_to_building_block_smiles_to_id[
                        reaction.chemical_space
                    ].get(molecule, "-1")
                    for molecule in molecules
                ),
            )

            product_nodes += [
                Node(
                    explore_weight=self.explore_weight,
                    scorer=self.scorer,
                    molecules=(product,),
                    unique_building_block_ids=node.unique_building_block_ids,
                    construction_log=ConstructionLog(
                        node.construction_log.reaction_logs + (reaction_log,)
                    ),
                    rollout_num=self.rollout_num,
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
            next_building_blocks = self.get_next_building_blocks(
                molecules=node.molecules
            )

        # Optionally, limit the number of next nodes
        if (
            self.num_expand_nodes is not None
            and len(next_building_blocks) > self.num_expand_nodes
        ):
            next_building_blocks = random_choice(
                rng=self.rng,
                array=next_building_blocks,
                size=self.num_expand_nodes,
                replace=False,
            )

        # Convert next node molecule tuples into Node objects
        new_nodes += [
            Node(
                explore_weight=self.explore_weight,
                scorer=self.scorer,
                molecules=node.molecules + (next_building_block,),
                unique_building_block_ids=node.unique_building_block_ids
                | {next_building_block},
                construction_log=node.construction_log,
                rollout_num=self.rollout_num,
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
        # Compute components of MCTS score
        exploit_score = node.exploit_score
        property_score = node.property_score
        explore_score = node.explore_score(
            n=total_visit_count, sign=1 if property_score >= 0 else -1
        )

        # Compute initial MCTS score (without building block diversity)
        mcts_score = exploit_score + property_score * explore_score

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

    def rollout(self, node: Node) -> Node:
        """Performs a generation rollout.

        :param node: A Node representing the current state of the rollout.
        :return: A Node containing a molecule with the maximum reward on the rollout.
        """
        if self.verbose:
            print(f"Node {node.node_id} (rollout {self.rollout_num})")
            print(f"Molecules = {node.molecules}")
            print(f"Num molecules = {node.num_molecules}")
            print(f"Num unique building blocks = {len(node.unique_building_block_ids)}")
            print(f"Num reactions = {node.num_reactions}")
            print(f"Score = {node.property_score}")
            print()

        # Stop the search if we've reached the maximum number of reactions
        if node.num_reactions >= self.max_reactions:
            return node

        # If this node has already been visited and the children have been stored, get its children from the dictionary
        if node in self.node_to_children:
            child_nodes = self.node_to_children[node]

        # Otherwise, expand the node to get its children
        else:
            # Expand the node both by running reactions with the current molecules and adding new building blocks
            child_nodes = self.get_child_nodes(node=node)

            # Check the node map and merge with an existing node if available
            child_nodes = [
                self.node_map.get(new_node, new_node) for new_node in child_nodes
            ]

            # Add nodes with complete molecules to the node map
            for child_node in child_nodes:
                if child_node.num_molecules == 1 and child_node not in self.node_map:
                    child_node.node_id = len(self.node_map)
                    self.node_map[child_node] = child_node
                    self.building_block_counts.update(
                        child_node.unique_building_block_ids
                    )

            # Save the number of children in order to maintain a total node count
            node.num_children = len(child_nodes)

            # If storing nodes, store the children
            if self.store_nodes:
                self.node_to_children[node] = child_nodes

        # If no new nodes were generated, return the current node
        if len(child_nodes) == 0:
            if node.num_molecules == 1:
                return node
            else:
                raise ValueError("Failed to expand a partially expanded node.")

        # Select a node based on the search type
        if self.search_type == "mcts":
            # Select node with the highest MCTS score
            total_visit_count = sum(child_node.num_visits for child_node in child_nodes)
            selected_node = max(
                child_nodes,
                key=partial(
                    self.compute_mcts_score, total_visit_count=total_visit_count
                ),
            )
        elif self.search_type == "rl":
            # Determine child nodes that are missing RL scores
            child_node_molecules = [child_node.molecules for child_node in child_nodes]
            child_node_molecules_missing_scores = [
                molecules
                for molecules in child_node_molecules
                if molecules not in self.molecules_to_rl_score
            ]

            # Compute RL scores for nodes that are missing RL scores
            if len(child_node_molecules_missing_scores) > 0:
                # Compute RL scores
                rl_scores = self.rl_model.predict(
                    molecule_tuples=child_node_molecules_missing_scores
                ).tolist()

                # Cache RL scores
                for molecules, rl_score in zip(
                    child_node_molecules_missing_scores, rl_scores
                ):
                    self.molecules_to_rl_score[molecules] = rl_score

            # Look up child node scores
            child_node_scores = np.array(
                [
                    self.molecules_to_rl_score[molecules]
                    for molecules in child_node_molecules
                ]
            )

            # Convert RL scores to temperature-scaled probabilities
            child_node_probs = softmax(child_node_scores / self.rl_temperature)

            # Select node proportional to the temperature-scaled RL score
            selected_node = self.rng.choice(child_nodes, p=child_node_probs)
        else:
            raise ValueError(f"Invalid search type: {self.search_type}")

        # Check the node map and merge with an existing node if available
        if selected_node in self.node_map:
            selected_node = self.node_map[selected_node]
        # Otherwise, assign node ID and add to node map
        else:
            selected_node.node_id = len(self.node_map)
            self.node_map[selected_node] = selected_node

        # Unroll the selected node
        best_node = self.rollout(node=selected_node)

        # Get full molecule (non-building block) with max score across rollouts
        if selected_node.num_molecules == 1 and node.num_reactions > 0:
            best_node = max(best_node, selected_node, key=lambda n: n.property_score)

        # Update exploit score and visit count
        selected_node.total_best_molecule_scores += best_node.property_score
        selected_node.num_visits += 1

        # Add RL training example
        if self.rl_model is not None:
            self.rl_model.buffer(source_node=selected_node, target_node=best_node)

        return best_node

    def compute_similarity(self, node: Node) -> float:
        """Computes the maximum Tanimoto similarity of the node to previously generated nodes.

        :param node: A Node with a single molecule.
        :return: The maximum Tanimoto similarity of the node to previously generated nodes.
        """
        # If there are no previously generated nodes, then the similarity is 0
        if self.full_molecule_morgan_fingerprints is None:
            return 0.0

        # Compute the Morgan fingerprint of the node
        node_morgan_fingerprint = morgan_fingerprint_generator(node.molecules[0])

        # Compute the Tanimoto similarity between the node and the previous molecules
        tanimoto_distances = pairwise_distances(
            node_morgan_fingerprint.reshape(1, -1),
            self.full_molecule_morgan_fingerprints,
            metric="jaccard",
            n_jobs=-1,
        )
        tanimoto_similarities = 1 - tanimoto_distances

        # Compute maximum similarity
        max_similarity = float(np.max(tanimoto_similarities))

        return max_similarity

    def update_full_molecule_fingerprints(self) -> None:
        """Updates the fingerprints of the full molecules with the new molecules from this generation."""
        # Get new full molecule nodes from this generation
        new_nodes = self.get_full_molecule_nodes(rollout_start=self.rollout_num)

        # If no new full molecules, then a duplicate was found
        if len(new_nodes) == 0:
            return

        # Get full molecule SMILES from this generation
        new_full_molecule_smiles = [node.molecules[0] for node in new_nodes]

        # Compute the Morgan fingerprints of the new nodes
        new_full_molecule_morgan_fingerprints = np.array(
            [
                morgan_fingerprint_generator(smiles)
                for smiles in new_full_molecule_smiles
            ]
        )

        # Set full molecule fingerprints the first time
        if self.full_molecule_morgan_fingerprints is None:
            self.full_molecule_morgan_fingerprints = (
                new_full_molecule_morgan_fingerprints
            )
            return

        # Update the full molecule fingerprints
        self.full_molecule_morgan_fingerprints = np.concatenate(
            [
                self.full_molecule_morgan_fingerprints,
                new_full_molecule_morgan_fingerprints,
            ],
            axis=0,
        )

    def get_full_molecule_nodes(
        self, rollout_start: int | None = None, rollout_end: int | None = None
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
            node for node in nodes if rollout_start <= node.rollout_num <= rollout_end
        ]

        # Sort Nodes by score and break ties by using Node ID
        nodes = sorted(
            nodes,
            key=lambda node: (node.property_score, -1 * node.node_id,),
            reverse=True,
        )

        return nodes

    def update_temperature(self, new_similarity: float) -> None:
        """Update the RL temperature based on the rolling average similarity.

        :param new_similarity: The similarity of the new molecules compared to previous molecules.
        """
        # Update rolling average similarity with weighted combination of new and old similarity
        self.rolling_average_similarity = (
            self.rolling_average_weight * self.rolling_average_similarity
            + (1 - self.rolling_average_weight) * new_similarity
        )

        # Determine percent difference between rolling average similarity and max similarity
        percent_similarity_difference = (
            self.rolling_average_similarity - self.rl_temperature_similarity_target
        ) / self.rl_temperature_similarity_target

        # Compute new temperature based on percent similarity difference
        new_temperature = (
            self.rl_temperature + percent_similarity_difference * self.rl_temperature
        )

        # Update temperature based on rolling average of new and old temperatures
        self.rl_temperature = (
            self.rolling_average_weight * self.rl_temperature
            + (1 - self.rolling_average_weight) * new_temperature
        )

        # Clip temperature within min/max bounds
        self.rl_temperature = max(
            self.min_temperature, min(self.rl_temperature, self.max_temperature)
        )

    def compute_successes(self, node: Node) -> list[int]:
        """Compute the successes of a generated molecule with respect to property success thresholds.

        :param node: A Node with the new generated molecule.
        :return: A list of successes (1/0) for each property success threshold.
        """
        # Compute scores of generated molecule and determine success rates
        scores = self.scorer.compute_individual_scores(smiles=node.molecules[0])
        successes = [
            int(success_comparator(score))
            for success_comparator, score in zip(self.success_comparators, scores)
        ]

        return successes

    def update_score_weights(self, successes: list[int]) -> None:
        """Update the score weights based on the rolling average success rate.

        :param successes: The successes (1/0) of the new molecule with respect to property success thresholds.
        """
        # Convert successes to NumPy
        successes = np.array(successes)

        # Update rolling average successes with weighted combination of new and old successes
        self.rolling_average_success_rate = (
            self.rolling_average_weight * self.rolling_average_success_rate
            + (1 - self.rolling_average_weight) * successes
        )  # (num_properties,)

        # Compute average success across properties
        average_success_rate = np.mean(self.rolling_average_success_rate)

        # Compute percent deviation from average success
        if average_success_rate > 0:
            percent_deviation_from_average = (
                self.rolling_average_success_rate - average_success_rate
            ) / average_success_rate
        else:
            percent_deviation_from_average = np.zeros(
                self.rolling_average_success_rate.shape
            )

        # Get current score weights
        weights = np.array(self.score_weights.weights)

        # Compute new score weights based on percent deviation from average success
        new_weights = weights - percent_deviation_from_average * weights

        # Compute score weights as weighted average of new and old weights
        weights = (
            self.rolling_average_weight * weights
            + (1 - self.rolling_average_weight) * new_weights
        )

        # Normalize weights
        weights /= np.sum(weights)

        # Ensure score weights exceed min bound
        weights = np.maximum(weights, self.min_score_weight)

        # Renormalize weights
        weights /= np.sum(weights)

        # Update score weights
        self.score_weights.weights = weights

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
            rollout_stats = {"Rollout Number": rollout_num}

            # Run rollout
            start_time = time.time()
            best_node = self.rollout(node=self.root)
            rollout_stats["Rollout Score"] = best_node.property_score
            rollout_stats["Rollout Time"] = time.time() - start_time

            # Log individual scores of new molecule
            for score_name, average_score in zip(
                self.score_weights.score_names, best_node.individual_scores
            ):
                rollout_stats[f"{score_name} Score"] = average_score

            # Compute similarity of new molecule compared to previous molecules
            start_time = time.time()
            new_similarity = self.compute_similarity(node=best_node)
            self.update_full_molecule_fingerprints()
            rollout_stats["Rollout Similarity"] = new_similarity
            rollout_stats["Similarity Time"] = time.time() - start_time

            # Determine number of unique full molecules found
            rollout_stats["Unique Molecules"] = sum(
                node.num_molecules == 1 and node.num_reactions > 0
                for node in self.node_map
            )

            # Optionally, update score weights based on success rate
            if self.success_comparators is not None:
                # Compute successes for new node
                successes = self.compute_successes(node=best_node)

                # Add successes to rollout stats
                for score_name, success in zip(
                    self.score_weights.score_names, successes
                ):
                    rollout_stats[f"{score_name} Success"] = success

                rollout_stats[f"Joint Success"] = int(int(all(successes)))

                # Update score weights
                self.update_score_weights(successes=successes)

            # Add score weights to rollout stats
            for score_name, score_weight in zip(
                self.score_weights.score_names, self.score_weights.weights
            ):
                rollout_stats[f"{score_name} Weight"] = score_weight

            # RL-specific updates and training
            if self.search_type == "rl":
                # Optionally, update temperature based on similarity of new molecules to previous molecules
                if self.rl_temperature_similarity_target is not None:
                    self.update_temperature(new_similarity=new_similarity)

                # Add RL temperature to rollout stats
                rollout_stats["RL Temperature"] = self.rl_temperature

                # Train and evaluate RL model
                if rollout_num % self.rl_train_frequency == 0:
                    # Reset RL scores since RL model is being updated
                    self.molecules_to_rl_score = {}

                    # Evaluate model on test set
                    start_time = time.time()
                    rollout_stats |= self.rl_model.evaluate(split="test")
                    rollout_stats["RL Test Eval Time"] = time.time() - start_time
                    rollout_stats["RL Test Examples"] = self.rl_model.test_size

                    # Move test set to train
                    self.rl_model.test_to_train()

                    # Train model
                    start_time = time.time()
                    self.rl_model.train()
                    rollout_stats["RL Train Time"] = time.time() - start_time

                    # Evaluate model on train set
                    start_time = time.time()
                    rollout_stats |= self.rl_model.evaluate(split="train")
                    rollout_stats["RL Train Eval Time"] = time.time() - start_time
                    rollout_stats["RL Train Examples"] = self.rl_model.train_size

            # Add rollout stats to record
            self.rollout_stats_record.append(rollout_stats)

            # Log rollout stats to W&B
            if self.wandb_log:
                wandb.log(rollout_stats)

        # Log rollout stats to file
        if self.log_path is not None:
            with open(self.log_path, "wb") as f:
                pickle.dump(self.rollout_stats_record, f)

        # Get all the Nodes representing fully constructed molecules within these rollouts sorted by score
        nodes = self.get_full_molecule_nodes(
            rollout_start=rollout_start, rollout_end=rollout_end
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
            raise ValueError(
                "Cannot get the precise number of nodes searched if store_nodes is False."
                "Use approx_num_nodes_searched instead."
            )

        # Get a set of all nodes and child nodes that have been visited
        visited_nodes = set()
        for node, children in self.node_to_children.items():
            visited_nodes.add(node)
            visited_nodes.update(children)

        return len(visited_nodes)
