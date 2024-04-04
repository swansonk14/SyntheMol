"""Test Enamine REAL reaction SMARTS for building blocks and products."""
from pathlib import Path

import numpy as np
import pandas as pd

from synthemol.reactions import REAL_REACTIONS
from synthemol.constants import (
    REAL_BUILDING_BLOCK_COLS,
    REAL_BUILDING_BLOCK_ID_COL,
    REAL_REACTION_COL,
    REAL_SMILES_COL,
)


def test_real_reactions(
    data_path: Path,
    building_blocks_path: Path,
    sample_num_per_reaction: int | None = None,
    show_mismatches_num: int = 0,
) -> None:
    """Test Enamine REAL reaction SMARTS for building blocks and products.

    :param data_path: Path to a CSV file containing a sample of REAL data with building block IDs and product SMILES.
    :param building_blocks_path: Path to a CSV file containing building blocks with IDs and SMILES.
    :param sample_num_per_reaction: Number of rows to sample from the data for each reaction. If None, use all rows.
    :param show_mismatches_num: Number of mismatches to show for each reactant and building block.
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
        print(f"Reaction {reaction.reaction_id}\n")

        # Get reaction data
        reaction_data = data[data[REAL_REACTION_COL] == reaction.reaction_id]

        # Sample reaction data
        if sample_num_per_reaction is not None and sample_num_per_reaction < len(
            reaction_data
        ):
            reaction_data = reaction_data.sample(sample_num_per_reaction, replace=False)

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
                # Determine which building blocks match the reactant
                match_mask = np.array(
                    [
                        reactant.has_substruct_match(building_block)
                        for building_block in reactant_building_blocks
                    ]
                )
                percent_matches = 100 * np.sum(match_mask) / len(match_mask)

                # Print statistics about the matches
                print(
                    f"Reactant {reactant_num + 1}, building block {building_block_num + 1}: "
                    f"{percent_matches:.2f}% matches"
                )

                # Optionally show mismatching building blocks
                if (
                    show_mismatches_num > 0
                    and building_block_num == reactant_num
                    and not np.all(match_mask)
                ):
                    print("Reactant SMARTS")
                    print(f"  {reactant.smarts}")
                    print(f"Example mismatching building blocks:")
                    mismatches = reactant_building_blocks[~match_mask]
                    for mismatch in mismatches.head(show_mismatches_num):
                        print(f"  {mismatch}")
            print()
        print()


if __name__ == "__main__":
    from tap import tapify

    tapify(test_real_reactions)
