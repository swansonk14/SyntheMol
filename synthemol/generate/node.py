"""Contains the Node class, which represents a step in the combinatorial molecule construction process."""
import math
from functools import cached_property
from typing import Any

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
        self.W = 0.0  # The sum of the leaf Node values for leaf Nodes that descend from this Node.
        self.N = 0  # The number of times this Node has been expanded.
        self.rollout_num = rollout_num
        self.num_children = 0

    @classmethod
    def compute_score(cls, molecules: tuple[str], scorer: MoleculeScorer) -> float:
        """Computes the score of the molecules.

        :param molecules: A tuple of SMILES. The first element is the currently constructed molecule
                          while the remaining elements are the building blocks that are about to be added.
        :param scorer: A callable object that takes as input a SMILES representing a molecule and returns a score.
        :return: The score of the molecules.
        """
        return (
            sum(scorer(molecule) for molecule in molecules) / len(molecules)
            if len(molecules) > 0
            else 0.0
        )

    @property
    def P(self) -> float:
        """The property score of this Node."""
        return self.compute_score(molecules=self.molecules, scorer=self.scorer)

    def Q(self) -> float:
        """Value that encourages exploitation of Nodes with high reward."""
        return self.W / self.N if self.N > 0 else 0.0

    def U(self, n: int) -> float:
        """Value that encourages exploration of Nodes with few visits."""
        return self.explore_weight * self.P * math.sqrt(1 + n) / (1 + self.N)

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
