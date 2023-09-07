"""Plots the average model score of building blocks vs the model score of the full molecule."""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import r2_score

from synthemol.constants import (
    BUILDING_BLOCKS_PATH,
    REAL_BUILDING_BLOCK_COLS,
    REAL_BUILDING_BLOCK_ID_COL,
    SCORE_COL
)


def plot_building_block_vs_molecule_scores(
        data_path: Path,
        save_dir: Path,
        title: str = 'Average Building Block Score vs Full Molecule Score',
        score_column: str = SCORE_COL,
        data_building_block_id_columns: list[str] = REAL_BUILDING_BLOCK_COLS,
        building_blocks_path: Path = BUILDING_BLOCKS_PATH,
        building_blocks_id_column: str = REAL_BUILDING_BLOCK_ID_COL,
        building_blocks_score_column: str = SCORE_COL
) -> None:
    """Plots the average model score of building blocks vs the model score of the full molecule.

    :param data_path: Path to a CSV file containing molecules, their building blocks, and full molecule scores.
    :param save_dir: Path to directory where the plot will be saved.
    :param title: Title of the plot to generate.
    :param score_column: Name of the model whose scores will be used.
    :param data_building_block_id_columns: Names of the columns in the data file containing building block IDs.
    :param building_blocks_path: Path to CSV file containing building block IDs and SMILES.
    :param building_blocks_id_column: Name of the column in the building blocks file containing building block IDs.
    :param building_blocks_score_column: Name of the column in the building blocks file containing building block scores.
    """
    # Load data
    data = pd.read_csv(data_path)
    print(f'Data size = {len(data):,}')

    # Load building block data
    building_block_data = pd.read_csv(building_blocks_path)

    # Map from building block ID to scores
    building_block_id_to_score = dict(zip(
        building_block_data[building_blocks_id_column],
        building_block_data[building_blocks_score_column]
    ))

    # Skip molecules with any missing building block SMILES
    missing_building_blocks = np.zeros(len(data), dtype=bool)

    for building_block_id_column in data_building_block_id_columns:
        missing_building_blocks |= np.array(
            [(building_block_id not in building_block_id_to_score and not np.isnan(building_block_id))
             for building_block_id in data[building_block_id_column]]
        )

    data = data[~missing_building_blocks]
    print(f'Data size after removing missing building blocks = {len(data):,}')

    # Determine number of building blocks in each molecule
    data['num_building_blocks'] = data[data_building_block_id_columns].notnull().sum(axis=1)

    # Analyze correlation stratified by number of building blocks in the full molecule
    save_dir.mkdir(parents=True, exist_ok=True)
    for num_building_blocks_in_full_molecule in [2, 3]:
        # Get data with the specified number of building blocks
        data_with_num_building_blocks = data[data['num_building_blocks'] == num_building_blocks_in_full_molecule]

        # Get the full molecule scores
        full_molecule_scores = data_with_num_building_blocks[score_column]

        # Analyze correlation stratified by number of building blocks predicting the full molecule score
        for num_building_blocks_to_average in range(1, num_building_blocks_in_full_molecule + 1):
            # Get the average building block score using the specified number of building blocks
            building_block_cols = data_building_block_id_columns[:num_building_blocks_to_average]

            # Get the relevant building blocks
            building_blocks = data_with_num_building_blocks[building_block_cols]

            # Look up the building block scores
            building_block_scores = building_blocks.applymap(building_block_id_to_score.get).mean(axis=1)

            # Compute R2
            r2 = r2_score(full_molecule_scores, building_block_scores)

            # Plot full vs average building block scores
            plt.clf()
            plt.scatter(full_molecule_scores, building_block_scores, s=3)
            plt.xlabel(f'Full Molecule Score (#BB = {num_building_blocks_in_full_molecule})')
            plt.ylabel(f'Average Building Block Score (#BB = {num_building_blocks_to_average})')
            plt.title(title)
            plt.text(0.98, 0.98, f'$R^2 = {r2:.3f}$',
                     horizontalalignment='right', verticalalignment='top',
                     transform=plt.gca().transAxes)

            # Save plot
            save_path = save_dir / (f'score_correlation_{num_building_blocks_to_average}_bbs_'
                                    f'vs_{num_building_blocks_in_full_molecule}_bb_full_molecule.pdf')
            plt.savefig(save_path, bbox_inches='tight')

            # Save data
            fig_data = pd.DataFrame({
                'full_molecule_score': full_molecule_scores,
                'building_block_score': building_block_scores,
            })
            fig_data.to_csv(save_path.with_suffix('.csv'), index=False)


if __name__ == '__main__':
    from tap import tapify

    tapify(plot_building_block_vs_molecule_scores)
