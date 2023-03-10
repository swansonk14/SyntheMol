"""Plots the distribution of model scores over molecular building blocks."""
import pickle
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def plot_building_block_scores(
        building_block_to_score_path: Path,
        title: str,
        save_dir: Path
) -> None:
    """Plots the distribution of model scores over molecular building blocks.

    :param building_block_to_score_path: Path to PKL file containing a dictionary mapping
                                         from building block SMILES to model scores.
    :param title: Title of the plot.
    :param save_dir: Path to a directory where the plot will be saved.
    """
    # Load mapping from building blocks to scores
    with open(building_block_to_score_path, 'rb') as f:
        building_block_to_score: dict[str, float] = pickle.load(f)

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
    fig_data = pd.DataFrame({'score': scores})
    fig_data.to_csv(save_dir / 'building_block_scores.csv', index=False)


if __name__ == '__main__':
    from tap import tapify

    tapify(plot_building_block_scores)
