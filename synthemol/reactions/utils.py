"""Utility functions for synthemol reactions."""
import pickle
from pathlib import Path

from synthemol.reactions.reaction import Reaction


def set_all_building_blocks(
    reactions: tuple[Reaction, ...], building_blocks: set[str]
) -> None:
    """Sets the allowed building block SMILES for all Reactions in a list of Reactions.

    Note: Modifies Reactions in place.

    :param reactions: A tuple of Reactions whose allowed building block SMILES will be set.
    :param building_blocks: A set of allowed building block SMILES.
    """
    for reaction in reactions:
        for reactant in reaction.reactants:
            reactant.all_building_blocks = building_blocks


def load_and_set_allowed_reaction_building_blocks(
    reactions: tuple[Reaction, ...],
    chemical_spaces: tuple[str, ...],
    reaction_to_building_blocks_paths: tuple[Path, ...],
) -> None:
    """Loads a mapping of allowed building blocks for each reaction and sets the allowed SMILES for each reaction.

    :param reactions: A tuple of Reactions whose allowed SMILES will be set.
    :param chemical_spaces: A tuple of chemical spaces that the reaction_to_building_blocks_paths belong to.
    :param reaction_to_building_blocks_paths: A tuple of paths to a PKL files mapping from reaction ID
                                              to reactant index to a set of allowed building block SMILES.
    """
    # Ensure that number of chemical_spaces matches number of reaction_to_building_blocks_paths
    if len(chemical_spaces) != len(reaction_to_building_blocks_paths):
        raise ValueError(
            "Number of chemical_spaces does not match number of reaction_to_building_blocks_paths."
        )

    # Load allowed building blocks for each reaction in each chemical space
    chemical_space_to_reaction_to_reactant_to_building_blocks: dict[
        str, dict[str, dict[int, set[str]]]
    ] = {}
    for chemical_space, reaction_to_building_blocks_path in zip(
        chemical_spaces, reaction_to_building_blocks_paths
    ):
        with open(reaction_to_building_blocks_path, "rb") as f:
            chemical_space_to_reaction_to_reactant_to_building_blocks[
                chemical_space
            ] = pickle.load(f)

    # Set allowed building blocks for each reactant in each reaction
    for reaction in reactions:
        for reactant_index, reactant in enumerate(reaction.reactants):
            reactant.allowed_building_blocks = chemical_space_to_reaction_to_reactant_to_building_blocks[
                reaction.chemical_space
            ][
                reaction.id
            ][
                reactant_index
            ]
