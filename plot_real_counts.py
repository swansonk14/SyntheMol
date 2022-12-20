"""Plot REAL reaction and reactant counts."""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from tap import Tap


class Args(Tap):
    reaction_counts_path: Path  # Path to a CSV file containing reaction counts.
    # reactant_counts_path: Path  # Path to a CSV file containing reactant counts.
    count_column: str = 'count'  # Name of the column containing counts.
    save_dir: Path  # Path to a directory where the plots will be saved.

    def process_args(self) -> None:
        self.save_dir.mkdir(parents=True, exist_ok=True)


def plot_real_counts(args: Args) -> None:
    """Plot REAL reaction and reactant counts."""
    # Load data
    reaction_counts = pd.read_csv(args.reaction_counts_path)
    # reactant_counts = pd.read_csv(args.reactant_counts_path)

    # Ensure counts are sorted
    reaction_counts.sort_values(args.count_column, ascending=False, inplace=True)
    # reactant_counts.sort_values(args.count_column, ascending=False, inplace=True)

    # Plot reaction counts
    plt.clf()
    plt.scatter(np.arange(len(reaction_counts)), np.cumsum(reaction_counts[args.count_column]), s=3)
    plt.xlabel('Reaction Index')
    plt.ylabel('Cumulative Molecule Count')
    plt.title('Cumulative Molecule Counts by Reaction')
    plt.savefig(args.save_dir / 'reaction_counts.pdf', bbox_inches='tight')

    # TODO: Plot reactant counts


if __name__ == '__main__':
    plot_real_counts(Args().parse_args())
