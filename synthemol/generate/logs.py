"""Contains the ReactionLog and ConstructionLog classes, which are used to log the construction of molecules."""
from typing import Iterator


class ReactionLog:
    """A ReactionLog contains information about a reaction and the building blocks combined by that reaction."""

    def __init__(self, chemical_space: str, reaction_id: int, reactant_ids: tuple[int, ...]) -> None:
        """Initializes the ReactionLog.

        :param chemical_space: The chemical space the reaction belongs to.
        :param reaction_id: The ID of the reaction.
        :param reactant_ids: A tuple of reactant IDs, which are building block IDs except for -1,
                             which indicates a non-building block molecule (i.e., a molecule already
                             containing multiple building blocks).
        """
        self.chemical_space = chemical_space
        self.reaction_id = reaction_id
        self.reactant_ids = reactant_ids

    @property
    def num_reactants(self) -> int:
        """Gets the number of reactants in the reaction."""
        return len(self.reactant_ids)

    @property
    def num_building_blocks(self) -> int:
        """Gets the number of building blocks in the reaction (i.e., all reactants with ID != -1)."""
        return sum(reactant_id != -1 for reactant_id in self.reactant_ids)


class ConstructionLog:
    """A ConstructionLog contains information about the reactions and building blocks used to construct a molecule.

    Note: The ConstructionLog is immutable.
    """

    def __init__(self, reaction_logs: tuple[ReactionLog, ...] | None = None) -> None:
        """Initializes the ConstructionLog.

        :param reaction_logs: A tuple of ReactionLogs.
        """
        self.reaction_logs = reaction_logs if reaction_logs is not None else tuple()

    @property
    def num_reactions(self) -> int:
        """Gets the number of reactions used to construct the molecule."""
        return len(self.reaction_logs)

    @property
    def num_building_blocks(self) -> int:
        """Gets the number of building blocks used to construct the molecule."""
        return sum(reaction_log.num_building_blocks for reaction_log in self.reaction_logs)

    def __len__(self) -> int:
        """Gets the number of reactions used to construct the molecule."""
        return len(self.reaction_logs)

    def __iter__(self) -> Iterator[ReactionLog]:
        """Iterates through the ReactionLogs."""
        return iter(self.reaction_logs)

    def __getitem__(self, item: int | slice) -> ReactionLog:
        """Gets the ReactionLog(s) at the specified index or slice."""
        return self.reaction_logs[item]
