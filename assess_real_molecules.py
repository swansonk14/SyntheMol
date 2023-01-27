"""Assess the quality and diversity of random REAL molecules."""
from pathlib import Path
from typing import Optional

import pandas as pd
from tap import Tap

from assess_generated_molecules import plot_scores, plot_similarity


class Args(Tap):
    data_path: Path  # Path to CSV file containing scores.
    save_dir: Path  # Path to directory where the plot will be saved.
    train_hits_path: Optional[Path] = None  # Optional path to CSV file containing hits from the training set for computing novelty.
    smiles_column: str = 'smiles'  # The name of the column containing SMILES in data_path.
    train_hits_smiles_column: str = 'smiles'  # The name of the column containing SMILES in train_hits_path.
    score_columns: list[str]  # Names of the columns containing scores.

    def process_args(self) -> None:
        self.save_dir.mkdir(parents=True, exist_ok=True)


def assess_real_molecules(args: Args) -> None:
    """Assess the quality and diversity of random REAL molecules."""
    # Load data
    data = pd.read_csv(args.data_path)

    # Count molecules
    print(f'Number of molecules = {len(data):,}')

    # Plot score distributions
    for score_column in args.score_columns:
        plot_scores(
            scores=data[score_column],
            save_dir=args.save_dir,
            score_name=score_column
        )

    # Similarity within REAL molecules
    plot_similarity(
        smiles=data[args.smiles_column],
        similarity_type='tanimoto',
        save_dir=args.save_dir
    )

    if args.train_hits_path is not None:
        # Load train hits
        train_hits = pd.read_csv(args.train_hits_path)

        # Similarity between REAL molecules and train hits
        plot_similarity(
            smiles=data[args.smiles_column],
            similarity_type='tversky',
            save_dir=args.save_dir,
            reference_smiles=train_hits[args.train_hits_smiles_column],
            reference_name='Train Hits'
        )


if __name__ == '__main__':
    assess_real_molecules(Args().parse_args())
