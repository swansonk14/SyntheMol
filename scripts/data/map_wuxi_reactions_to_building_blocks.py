"""Determines which WuXi building blocks can be used in which WuXi reactions."""
import pickle
from pathlib import Path

import pandas as pd
from tqdm import tqdm

from synthemol.constants import (
    WUXI_BUILDING_BLOCK_SMILES_COL,
    WUXI_BUILDING_BLOCK_SUBSET_COL
)
from synthemol.reactions import WUXI_REACTIONS


def map_wuxi_reactions_to_building_blocks(
        building_blocks_path: Path,
        save_path: Path,
        smiles_column: str = WUXI_BUILDING_BLOCK_SMILES_COL,
        subset_column: str = WUXI_BUILDING_BLOCK_SUBSET_COL
) -> None:
    """Determines which WuXi building blocks can be used in which WuXi reactions.

    :param building_blocks_path: Path to a CSV file containing building blocks.
    :param save_path: Path to a PKL file where the mapping will be saved.
    :param smiles_column: The column name for building block SMILES.
    :param subset_column: The column name for the building block subset.
    """
    # Load data
    data = pd.read_csv(building_blocks_path)

    # Create mapping from reaction ID to reactant ID to list of building block SMILES
    # TODO: implement reactant.subsets
    reaction_to_reactant_to_building_blocks: dict[int, dict[int, set[str]]] = {
        reaction.id: {
            reactant_index: sorted({
                building_block
                for building_block in data[data[subset_column].isin(reactant.subsets)][smiles_column]
                if reactant.has_substruct_match(building_block)
            })
            for reactant_index, reactant in enumerate(reaction.reactants)
        }
        for reaction in tqdm(WUXI_REACTIONS)
    }

    # Save mapping
    with open(save_path, 'wb') as f:
        pickle.dump(reaction_to_reactant_to_building_blocks, f)


if __name__ == '__main__':
    from tap import tapify

    tapify(map_wuxi_reactions_to_building_blocks)
