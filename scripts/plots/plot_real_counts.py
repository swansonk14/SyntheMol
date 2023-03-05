"""Plot REAL reaction and reactant counts."""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from tap import tapify


def plot_real_counts(
        reaction_counts_path: Path,
        building_block_counts_path: Path,
        save_dir: Path,
        count_column: str = 'count'
) -> None:
    """Plot REAL reaction and reactant counts.

    :param reaction_counts_path: Path to a CSV file containing reaction counts.
    :param building_block_counts_path: Path to a CSV file containing reagent counts.
    :param save_dir: Path to a directory where the plots will be saved.
    :param count_column: Name of the column containing counts.
    """
    # Load data
    reaction_counts = pd.read_csv(reaction_counts_path)
    building_block_counts = pd.read_csv(building_block_counts_path)

    # Ensure counts are sorted
    reaction_counts.sort_values(count_column, ascending=False, inplace=True)
    building_block_counts.sort_values(count_column, ascending=False, inplace=True)

    # Create save directory
    save_dir.mkdir(parents=True, exist_ok=True)

    # Plot reaction counts
    plt.clf()
    plt.scatter(np.arange(len(reaction_counts)), np.cumsum(reaction_counts[count_column]), s=3)
    plt.xlabel('Reaction Index')
    plt.ylabel('Cumulative Molecule Count')
    plt.title('REAL Space Reaction Counts')
    plt.savefig(save_dir / 'reaction_counts.pdf', bbox_inches='tight')

    # Save reaction counts
    reaction_counts.to_csv(save_dir / 'reaction_counts.csv', index=False)

    # Plot building block counts
    plt.clf()
    plt.scatter(np.arange(len(building_block_counts)), building_block_counts[count_column], s=3)
    plt.xlabel('Building Block Index')
    plt.ylabel('Number of Molecules with Building Block')
    plt.title('REAL Space Building Block Counts')
    plt.savefig(save_dir / 'building_block_counts.pdf', bbox_inches='tight')

    # Save building block counts
    building_block_counts.to_csv(save_dir / 'building_block_counts.csv', index=False)


if __name__ == '__main__':
    tapify(plot_real_counts)
