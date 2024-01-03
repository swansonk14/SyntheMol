"""Filter out REAL building blocks that do not match the REAL reaction templates and map BB IDs to SMILES."""
import pickle
from pathlib import Path

import pandas as pd
from tqdm import tqdm

from synthemol.reactions import Reaction, REAL_REACTIONS
from synthemol.constants import (
    BUILDING_BLOCKS_PATH,
    REAL_BUILDING_BLOCK_ID_COL,
    SMILES_COL,
)


def filter_real_reactions_to_building_blocks(
    reaction_to_building_blocks_path: Path,
    save_path: Path,
    reactions: tuple[Reaction, ...] = REAL_REACTIONS,
    building_blocks_path: Path = BUILDING_BLOCKS_PATH,
    building_blocks_id_column: str = REAL_BUILDING_BLOCK_ID_COL,
    building_blocks_smiles_column: str = SMILES_COL,
) -> None:
    """Filter out REAL building blocks that do not match the REAL reaction templates and map BB IDs to SMILES.

    :param reaction_to_building_blocks_path: Path to a PKL file mapping from reaction ID
        to reactant index to a set of allowed building block IDs.
    :param save_path: Path to a PKL file where the filtered reaction to building blocks will be saved.
    :param reactions: A tuple of Reactions whose allowed building blocks will be determined.
    :param building_blocks_path: Path to a CSV file containing building block IDs and SMILES.
    :param building_blocks_id_column: The column name for building block IDs.
    :param building_blocks_smiles_column: The column name for building block SMILES.
    """
    # Load allowed building blocks for each reaction
    with open(reaction_to_building_blocks_path, "rb") as f:
        reaction_to_reactant_to_building_block_ids: dict[
            int, dict[int, set[int]]
        ] = pickle.load(f)

    # Load building blocks
    building_block_data = pd.read_csv(building_blocks_path)

    # Get building block IDs and SMILES
    building_block_id_to_smiles = dict(
        zip(
            building_block_data[building_blocks_id_column],
            building_block_data[building_blocks_smiles_column],
        )
    )

    # Convert building block IDs to SMILES
    reaction_to_reactant_to_building_blocks = {
        reaction: {
            reactant: {
                building_block_id_to_smiles[building_block_id]
                for building_block_id in building_block_ids
                if building_block_id in building_block_id_to_smiles
            }
            for reactant, building_block_ids in reactant_to_building_block_ids.items()
        }
        for reaction, reactant_to_building_block_ids in reaction_to_reactant_to_building_block_ids.items()
    }

    # Filter out building blocks and reactions that do not match the given reaction templates
    reaction_to_reactant_to_building_blocks_filtered = {
        reaction.id: {
            reactant_index: {
                smiles
                for smiles in reaction_to_reactant_to_building_blocks[reaction.id][
                    reactant_index
                ]
                if reactant.has_substruct_match(smiles)
            }
            for reactant_index, reactant in enumerate(reaction.reactants)
        }
        for reaction in tqdm(reactions)
    }

    # Save filtered reaction to building blocks
    with open(save_path, "wb") as f:
        pickle.dump(reaction_to_reactant_to_building_blocks_filtered, f)


if __name__ == "__main__":
    from tap import tapify

    tapify(filter_real_reactions_to_building_blocks)
