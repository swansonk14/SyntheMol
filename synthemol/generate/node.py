"""Contains the Node class, which represents a step in the combinatorial molecule construction process."""
import math
from functools import cached_property
from typing import Any, Literal

import numpy as np

from synthemol.generate.logs import ConstructionLog
from synthemol.generate.scorer import MoleculeScorer


class Node:
    """A Node represents a step in the combinatorial molecule construction process."""

    def __init__(
        self,
        explore_weight: float,
        scorer: MoleculeScorer,
        node_id: int | None = None,
        molecules: tuple[str] | None = None,
        unique_building_block_ids: set[int] | None = None,
        construction_log: ConstructionLog = None,
        rollout_num: int | None = None,
    ) -> None:
        """Initializes the Node.

        :param explore_weight: The hyperparameter that encourages exploration.
        :param scorer: A callable object that takes as input a SMILES representing a molecule and returns a score.
        :param node_id: The ID of the Node, which should correspond to the order in which the Mode was visited.
        :param molecules: A tuple of SMILES. The first element is the currently constructed molecule
                          while the remaining elements are the building blocks that are about to be added.
        :param unique_building_block_ids: A set of building block IDS used in this Node.
        :param construction_log: A ConstructionLog containing information about each reaction
                                 used to construct the molecules in this Node.
        :param rollout_num: The number of the rollout on which this Node was created.
        """
        self.explore_weight = explore_weight
        self.scorer = scorer
        self.node_id = node_id
        self.molecules = molecules if molecules is not None else tuple()
        self.unique_building_block_ids = (
            unique_building_block_ids
            if unique_building_block_ids is not None
            else set()
        )
        self.construction_log = (
            construction_log if construction_log is not None else ConstructionLog()
        )
        self.total_best_molecule_scores = np.zeros(self.scorer.num_scores)
        self.num_visits = 0
        self.rollout_num = rollout_num
        self.num_children = 0

    @cached_property
    def individual_scores(self) -> list[float]:
        """Computes the average individual scores of the Node's molecules.

        Note: Individual scores are assumed to be immutable.

        :return: The average individual scores of the Node's molecules.
        """
        if self.num_molecules == 0:
            return [0.0] * self.scorer.num_scores

        # Compute individual scores for each molecule
        individual_scores = np.array(
            [
                self.scorer.compute_individual_scores(molecule)
                for molecule in self.molecules
            ]
        )  # (num_molecules, num_scores)

        # Average individual scores for the molecules
        average_individual_scores = np.mean(
            individual_scores, axis=0
        ).tolist()  # (num_scores,)

        return average_individual_scores

    @property
    def property_score(self) -> float:
        """Computes the average property score of the Node's molecules (weighted average of individual property scores).

        :return: The average property score of the Node's molecules (weighted average of individual property scores).
        """
        return sum(
            weight * score
            for weight, score in zip(
                self.scorer.score_weights.weights, self.individual_scores
            )
        )

    @property
    def exploit_score(self) -> float:
        """Value that encourages exploitation of Nodes with high reward."""
        # Return 0 if no visits
        if self.num_visits == 0:
            return 0.0

        # Compute average best molecule scores
        average_best_molecule_scores = (
            self.total_best_molecule_scores / self.num_visits
        )  # (num_scores,)

        # Get score weights
        score_weights = np.array(self.scorer.score_weights.weights)  # (num_scores,)

        # Compute weighted average of model scores
        weighted_average_best_molecule_scores = np.dot(
            average_best_molecule_scores, score_weights
        )

        return weighted_average_best_molecule_scores

    def explore_score(self, n: int, sign: Literal[1, -1] = 1) -> float:
        """Value that encourages exploration of Nodes with few visits.

        :param n: The total number of times that nodes on the same level have been visited.
        :param sign: The sign of the property score, which determines whether the explore
                     score or its reciprocal should be used since the two are multiplied.
        """
        return (self.explore_weight * math.sqrt(1 + n) / (1 + self.num_visits)) ** sign

    @property
    def num_molecules(self) -> int:
        """Gets the number of building blocks in the Node."""
        return len(self.molecules)

    @property
    def num_reactions(self) -> int:
        """Gets the number of reactions used so far to generate the molecule in the Node."""
        return self.construction_log.num_reactions

    @property
    def num_building_blocks(self) -> int:
        """Gets the number of building blocks used to generate the molecule in the Node."""
        return (
            self.construction_log.num_building_blocks
            + self.num_molecules
            - (self.num_reactions > 0)
        )

    def __hash__(self) -> int:
        """Hashes the Node based on the building blocks."""
        return hash(self.molecules)

    def __eq__(self, other: Any) -> bool:
        """Checks if the Node is equal to another Node based on the building blocks."""
        if not isinstance(other, Node):
            return False

        return self.molecules == other.molecules
