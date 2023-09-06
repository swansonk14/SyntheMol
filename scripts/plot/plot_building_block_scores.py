"""Plots the distribution of model scores over molecular building blocks."""
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

from synthemol.constants import SCORE_COL, SMILES_COL


def plot_building_block_scores(
        building_blocks_path: Path,
        title: str,
        save_dir: Path,
        building_blocks_smiles_column: str = SMILES_COL,
        building_blocks_score_column: str = SCORE_COL
) -> None:
    """Plots the distribution of model scores over molecular building blocks.

    :param building_blocks_path: Path to CSV file containing building blocks and scores.
    :param title: Title of the plot.
    :param save_dir: Path to a directory where the plot will be saved.
    :param building_blocks_smiles_column: Name of the column in the building blocks file containing SMILES.
    :param building_blocks_score_column: Name of the column in the building blocks file containing building block scores.

    """
    # Load mapping from building blocks to scores
    data = pd.read_csv(building_blocks_path)
    building_block_to_score = dict(zip(data[building_blocks_smiles_column], data[building_blocks_score_column]))

    # Plot distribution of building block scores
    scores = list(building_block_to_score.values())

    plt.hist(scores, bins=100)
    plt.xlabel('Model Score')
    plt.ylabel('Count')
    plt.title(title)

    # Save plot
    save_dir.mkdir(parents=True, exist_ok=True)
    plt.savefig(save_dir / 'building_block_scores.pdf', bbox_inches='tight')

    # Save data
    fig_data = pd.DataFrame({SCORE_COL: scores})
    fig_data.to_csv(save_dir / 'building_block_scores.csv', index=False)


if __name__ == '__main__':
    from tap import tapify

    tapify(plot_building_block_scores)
