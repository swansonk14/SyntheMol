"""Plots the average model score of molecular fragments vs the model score of the full molecule."""
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import r2_score
from tap import Tap


class Args(Tap):
    data_path: Path  # Path to a CSV file containing molecules, their molecular fragments, and full molecule scores.
    score_column: str  # Name of the model whose scores will be used.
    fragment_to_score_path: Path  # Path to JSON file containing a dictionary mapping from fragment SMILES to model scores.
    title: str  # Title of the plot to generate
    save_path: Path  # Path to PDF or PNG file where the plot will be saved.

    def process_args(self) -> None:
        self.save_path.parent.mkdir(parents=True, exist_ok=True)


def plot_fragment_vs_molecule_scores(args: Args) -> None:
    """Plots the average model score of molecular fragments vs the model score of the full molecule."""
    # Load data
    data = pd.read_csv(args.data_path)

    # Load mapping from fragments to scores
    with open(args.fragment_to_score_path) as f:
        fragment_to_score: dict[str, float] = json.load(f)

    # Get full molecule scores
    full_molecule_scores = data[args.score_column]

    # Get average fragment scores
    fragment_columns = [
        column
        for column in data.columns
        if column.startswith('reagent') and column.endswith('smiles')
    ]
    fragment_scores = [
        np.mean([fragment_to_score[fragment] for fragment in row[fragment_columns].dropna()])
        for _, row in data.iterrows()
    ]

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
