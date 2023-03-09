"""Plots the average model score of building blocks vs the model score of the full molecule."""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import r2_score
from tqdm import tqdm

from SyntheMol.constants import (
    REAL_BUILDING_BLOCK_COLS,
    REAL_BUILDING_BLOCK_ID_COL,
    SCORE_COL
)


def plot_building_block_vs_molecule_scores(
        data_path: Path,
        score_column: str,
        building_blocks_path: Path,
        title: str,
        save_dir: Path,
        building_blocks_id_column: str = REAL_BUILDING_BLOCK_ID_COL,
        building_blocks_score_column: str = SCORE_COL
) -> None:
    """Plots the average model score of building blocks vs the model score of the full molecule.

    :param data_path: Path to a CSV file containing molecules, their building blocks, and full molecule scores.
    :param score_column: Name of the model whose scores will be used.
    :param building_blocks_path: Path to CSV file containing building block IDs and SMILES.
    :param title: Title of the plot to generate.
    :param save_dir: Path to directory where the plot will be saved.
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

    for building_block_id_column in REAL_BUILDING_BLOCK_COLS:
        missing_building_blocks |= np.array(
            [(building_block_id not in building_block_id_to_score and not np.isnan(building_block_id))
             for building_block_id in data[building_block_id_column]]
        )

    data = data[~missing_building_blocks]
    print(f'Data size after removing missing building blocks = {len(data):,}')

    # Get full molecule scores
    full_molecule_scores = data[score_column]

    # Get average building block scores
    building_block_scores = [
        np.mean([
            building_block_id_to_score[building_block]
            for building_block in row[REAL_BUILDING_BLOCK_COLS].dropna()
        ])
        for _, row in tqdm(data.iterrows(), total=len(data), desc='Average building block scores')
    ]

    # Compute R2
    r2 = r2_score(full_molecule_scores, building_block_scores)

    # Plot full vs average building block scores
    plt.scatter(full_molecule_scores, building_block_scores, s=3)
    plt.xlabel('Full Molecule Score')
    plt.ylabel('Average Building Block Score')
    plt.title(title)
    plt.text(0.98, 0.98, f'$R^2 = {r2:.3f}$',
             horizontalalignment='right', verticalalignment='top',
             transform=plt.gca().transAxes)

    # Save plot
    save_dir.mkdir(parents=True, exist_ok=True)
    plt.savefig(save_dir / 'building_block_vs_full_scores.pdf', bbox_inches='tight')

    # Save data
    fig_data = pd.DataFrame({
        'full_molecule_score': full_molecule_scores,
        'building_block_score': building_block_scores,
    })
    fig_data.to_csv(save_dir / 'building_block_vs_full_scores.csv', index=False)


if __name__ == '__main__':
    from tap import tapify

    tapify(plot_building_block_vs_molecule_scores)
