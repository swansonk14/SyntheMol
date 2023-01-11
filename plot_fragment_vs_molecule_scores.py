"""Plots the average model score of molecular fragments vs the model score of the full molecule."""
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import r2_score
from tap import Tap
from tqdm import tqdm

from constants import REAL_REAGENT_COLS


class Args(Tap):
    data_path: Path  # Path to a CSV file containing molecules, their molecular fragments, and full molecule scores.
    score_column: str  # Name of the model whose scores will be used.
    fragment_path: Path  # Path to CSV file containing fragment IDs and SMILES.
    fragment_to_score_path: Path  # Path to JSON file containing a dictionary mapping from fragment SMILES to model scores.
    title: str  # Title of the plot to generate
    save_dir: Path  # Path to directory where the plot will be saved.

    def process_args(self) -> None:
        self.save_dir.mkdir(parents=True, exist_ok=True)


def plot_fragment_vs_molecule_scores(args: Args) -> None:
    """Plots the average model score of molecular fragments vs the model score of the full molecule."""
    # Load data
    data = pd.read_csv(args.data_path)
    print(f'Data size = {len(data):,}')

    # Load fragment data
    fragment_data = pd.read_csv(args.fragment_path)

    # Load mapping from fragments to scores
    with open(args.fragment_to_score_path) as f:
        fragment_to_score: dict[str, float] = json.load(f)

    # Map from fragment ID to scores
    fragment_id_to_score = {
        fragment_id: fragment_to_score[fragment_smiles]
        for fragment_id, fragment_smiles in zip(fragment_data['Reagent_ID'], fragment_data['smiles'])
    }

    # Skip molecules with any missing fragment SMILES
    missing_fragments = np.zeros(len(data), dtype=bool)

    for reagent_column in REAL_REAGENT_COLS:
        missing_fragments |= np.array(
            [(reagent not in fragment_id_to_score and not np.isnan(reagent)) for reagent in data[reagent_column]]
        )

    data = data[~missing_fragments]
    print(f'Data size after removing missing fragments = {len(data):,}')

    # Get full molecule scores
    full_molecule_scores = data[args.score_column]

    # Get average fragment scores
    fragment_scores = [
        np.mean([fragment_id_to_score[fragment] for fragment in row[REAL_REAGENT_COLS].dropna()])
        for _, row in tqdm(data.iterrows(), total=len(data), desc='Average fragment scores')
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
    plt.savefig(args.save_dir / 'fragment_vs_full_scores.pdf', bbox_inches='tight')

    # Save data
    fig_data = pd.DataFrame({
        'full_molecule_score': full_molecule_scores,
        'fragment_score': fragment_scores,
    })
    fig_data.to_csv(args.save_dir / 'fragment_vs_full_scores.csv', index=False)


if __name__ == '__main__':
    plot_fragment_vs_molecule_scores(Args().parse_args())
