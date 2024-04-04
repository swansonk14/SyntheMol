"""Test Enamine REAL reaction SMARTS for building blocks and products."""
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem

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
    building_blocks_data = pd.read_csv(building_blocks_path)

    # Set building block IDs as index
    building_blocks_data.set_index(REAL_BUILDING_BLOCK_ID_COL, inplace=True)

    # Map building block IDs to SMILES in data
    for building_block_col in REAL_BUILDING_BLOCK_COLS:
        data[building_block_col] = data[building_block_col].map(
            building_blocks_data[REAL_SMILES_COL]
        )

    # Test each reaction
    for reaction in REAL_REACTIONS:
        print(f"Reaction {reaction.reaction_id}\n")

        # Get reaction data
        reaction_data = data[data[REAL_REACTION_COL] == reaction.reaction_id]

        # Get building block columns
        building_block_cols = REAL_BUILDING_BLOCK_COLS[: len(reaction.reactants)]

        # Identify any NaNs in building blocks
        nan_mask = reaction_data[building_block_cols].isna().any(axis=1)

        if nan_mask.any():
            print(
                f"Warning: {nan_mask.sum()} rows contain NaNs in building blocks and will be skipped."
            )
            reaction_data = reaction_data[~nan_mask]

        # Sample reaction data
        if sample_num_per_reaction is not None and sample_num_per_reaction < len(
            reaction_data
        ):
            reaction_data = reaction_data.sample(
                sample_num_per_reaction, random_state=0, replace=False
            )

        # Check if there is any data for this reaction
        if len(reaction_data) == 0:
            print("No data for this reaction.\n")
            continue

        # For each reactant in the reaction, test all the building blocks against each reactant SMARTS
        for building_block_num, building_block_col in enumerate(building_block_cols):
            # Get reactant building blocks
            reactant_building_blocks = reaction_data[building_block_col]

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
                    num_show = min(show_mismatches_num, len(mismatches))
                    for mismatch in mismatches.sample(
                        num_show, random_state=0, replace=False
                    ):
                        print(f"  {mismatch}")
            print()
        print()

        # Get all the reaction data where every building block matches the reactant SMARTS
        match_mask = np.array(
            [
                reaction_data[building_block_col].apply(reactant.has_substruct_match)
                for reactant, building_block_col in zip(
                    reaction.reactants, REAL_BUILDING_BLOCK_COLS
                )
            ]
        )

        matching_reaction_data = reaction_data[np.all(match_mask, axis=0)]

        # For each row of reaction data, compare the real product to the generated product from running the reaction
        for building_blocks, product in zip(
            matching_reaction_data[REAL_BUILDING_BLOCK_COLS].itertuples(index=False),
            matching_reaction_data[REAL_SMILES_COL],
        ):
            # Canonicalize product
            product = Chem.MolToSmiles(Chem.MolFromSmiles(product))

            # Run the reaction on the building blocks
            generated_products = reaction.run_reactants(building_blocks)

            # Check if the generated product matches the real product
            if product not in generated_products:
                print("Mismatch")
                print(f"Real product: {product}")
                print(f"Generated products: {generated_products}")
                print()


if __name__ == "__main__":
    from tap import tapify

    tapify(test_real_reactions)
