"""Test Enamine REAL reaction SMARTS for building blocks and products."""
from pathlib import Path

import pandas as pd

from synthemol.reactions import REAL_REACTIONS
from synthemol.constants import (
    REAL_BUILDING_BLOCK_COLS,
    REAL_BUILDING_BLOCK_ID_COL,
    REAL_REACTION_COL,
    REAL_SMILES_COL,
)


def test_real_reactions(data_path: Path, building_blocks_path: Path) -> None:
    """Test Enamine REAL reaction SMARTS for building blocks and products.

    :param data_path: Path to a CSV file containing a sample of REAL data with building block IDs and product SMILES.
    :param building_blocks_path: Path to a CSV file containing building blocks with IDs and SMILES.
    """
    # Load data
    data = pd.read_csv(data_path)
    building_blocks = pd.read_csv(building_blocks_path)

    # Set building block IDs as index
    building_blocks.set_index(REAL_BUILDING_BLOCK_ID_COL, inplace=True)

    # Map building block IDs to SMILES in data
    for building_block_col in REAL_BUILDING_BLOCK_COLS:
        data[building_block_col] = data[building_block_col].map(
            building_blocks[REAL_SMILES_COL]
        )

    # Test each reaction
    for reaction in REAL_REACTIONS:
        print(f"Reaction {reaction.reaction_id}")

        # Get reaction data
        reaction_data = data[data[REAL_REACTION_COL] == reaction.reaction_id]

        # For each reactant in the reaction, test all the building blocks against each reactant SMARTS
        for building_block_num, building_block_col in enumerate(
            REAL_BUILDING_BLOCK_COLS
        ):
            # Get reactant building blocks
            reactant_building_blocks = reaction_data[building_block_col].dropna()

            # Skip if no building blocks for the reactant
            if reactant_building_blocks.empty:
                continue

            # For each reactant in the reaction, test all the building blocks
            for reactant_num, reactant in enumerate(reaction.reactants):
                num_matches = sum(
                    reactant.has_substruct_match(building_block)
                    for building_block in reactant_building_blocks
                )

                print(
                    f"Reactant {reactant_num + 1} building block {building_block_num + 1}: "
                    f"{num_matches} ({100 * num_matches / len(reactant_building_blocks)}%) matches"
                )


if __name__ == "__main__":
    from tap import tapify

    tapify(test_real_reactions)
