"""Utility functions for synthemol reactions."""
import pickle
from pathlib import Path

from tqdm import tqdm

from synthemol.reactions.reaction import Reaction


def set_all_building_blocks(
        reactions: list[Reaction],
        building_blocks: set[str]
) -> None:
    """Sets the allowed building block SMILES for all Reactions in a list of Reactions.

    Note: Modifies Reactions in place.

    :param reactions: A list of Reactions whose allowed building block SMILES will be set.
    :param building_blocks: A set of allowed building block SMILES.
    """
    for reaction in reactions:
        for reactant in reaction.reactants:
            reactant.all_building_blocks = building_blocks


def set_allowed_reaction_building_blocks(
        reactions: list[Reaction],
        reaction_to_reactant_to_building_blocks: dict[int, dict[int, set[str]]]
) -> None:
    """Sets the allowed building block SMILES for each reactant in each Reaction in a list of Reactions.

    Note: Modifies Reactions in place.

    :param reactions: A list of Reactions whose allowed building block SMILES will be set.
    :param reaction_to_reactant_to_building_blocks: A dictionary mapping from reaction ID
                                                    to reactant index to a set of allowed SMILES.
    """
    for reaction in tqdm(reactions):
        for reactant_index, reactant in enumerate(reaction.reactants):
            reactant.allowed_building_blocks = reaction_to_reactant_to_building_blocks[reaction.id][reactant_index]


def load_and_set_allowed_reaction_building_blocks(
        reactions: list[Reaction],
        reaction_to_reactant_to_building_blocks_path: Path,
        building_block_id_to_smiles: dict[int, str],
) -> None:
    """Loads a mapping of allowed building blocks for each reaction and sets the allowed SMILES for each reaction.

    :param reactions: A list of Reactions whose allowed SMILES will be set.
    :param reaction_to_reactant_to_building_blocks_path: Path to a PKL file mapping from reaction ID
                                                            to reactant index to a set of allowed building block IDs.
    :param building_block_id_to_smiles: A dictionary mapping from building block ID to SMILES.
    """
    # Load allowed building blocks for each reaction
    with open(reaction_to_reactant_to_building_blocks_path, 'rb') as f:
        reaction_to_reactant_to_building_block_ids: dict[int, dict[int, set[int]]] = pickle.load(f)

    # Convert building block IDs to SMILES
    reaction_to_reactant_to_building_blocks = {
        reaction: {
            reactant: {
                building_block_id_to_smiles[building_block_id]
                for building_block_id in building_block_ids
                if building_block_id in building_block_id_to_smiles
            }
            for reactant, building_block_ids in reactant_to_building_block_ids.items()
        } for reaction, reactant_to_building_block_ids in reaction_to_reactant_to_building_block_ids.items()
    }

    # Set allowed building blocks for each reaction
    set_allowed_reaction_building_blocks(
        reactions=reactions,
        reaction_to_reactant_to_building_blocks=reaction_to_reactant_to_building_blocks
    )
