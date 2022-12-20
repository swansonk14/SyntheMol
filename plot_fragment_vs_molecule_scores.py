"""Plots the average model score of molecular fragments vs the model score of the full molecule."""
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from sklearn.metrics import r2_score
from tap import Tap


class Args(Tap):
    data_path: Path  # Path to a CSV file containing molecules, their molecular fragments, and full molecule scores.
    score_name: str  # Name of the model whose scores will be used.
    title: str  # Title of the plot to generate
    save_path: Path  # Path to PDF or PNG file where the plot will be saved.

    def process_args(self) -> None:
        self.save_path.parent.mkdir(parents=True, exist_ok=True)


def plot_fragment_vs_molecule_scores(args: Args) -> None:
    """Plots the average model score of molecular fragments vs the model score of the full molecule."""
    # Load data
    data = pd.read_csv(args.data_path)

    # Get full molecule scores
    full_molecule_scores = data[args.score_name]

    # Get average fragment scores
    fragment_score_columns = [
        column
        for column in data.columns
        if column.startswith('reagent') and column.endswith(args.score_name)
    ]
    fragment_scores = data[fragment_score_columns].mean(axis=1)

    # Compute R2
    r2 = r2_score(full_molecule_scores, fragment_scores)

    # Plot full vs average fragment scores
    plt.scatter(full_molecule_scores, fragment_scores, s=3)
    plt.xlabel('Full Molecule Score')
    plt.ylabel('Average Fragment Score')
    plt.title(args.title)
    plt.text(0.98, 0.98, f'$R^2 = {r2:.3f}$',
             horizontalalignment='right', verticalalignment='top',
             transform=plt.gca().transAxes)
    plt.savefig(args.save_path, bbox_inches='tight')


if __name__ == '__main__':
    plot_fragment_vs_molecule_scores(Args().parse_args())
