"""Assess the novelty, scores, and diversity of molecules."""
from pathlib import Path

import pandas as pd

from SyntheMol.constants import SMILES_COL
from plot_generated_molecule_analysis import plot_scores, plot_similarity


def plot_molecule_analysis(
        data_path: Path,
        save_dir: Path,
        score_columns: list[str],
        train_hits_path: Path | None = None,
        smiles_column: str = SMILES_COL,
        train_hits_smiles_column: str = SMILES_COL,
) -> None:
    """Assess the novelty, scores, and diversity of molecules.

    :param data_path: Path to CSV file containing scores.
    :param save_dir: Path to directory where the plot will be saved.
    :param score_columns: Names of the columns containing scores.
    :param train_hits_path: Optional path to CSV file containing hits from the training set for computing novelty.
    :param smiles_column: The name of the column containing SMILES in data_path.
    :param train_hits_smiles_column: The name of the column containing SMILES in train_hits_path.
    """
    # Load data
    data = pd.read_csv(data_path)

    # Count molecules
    print(f'Number of molecules = {len(data):,}')

    # Create save directory
    save_dir.mkdir(parents=True, exist_ok=True)

    # Plot score distributions
    for score_column in score_columns:
        plot_scores(
            scores=data[score_column],
            save_dir=save_dir,
            score_name=score_column
        )

    # Similarity within REAL molecules
    plot_similarity(
        smiles=data[smiles_column],
        similarity_type='tanimoto',
        save_dir=save_dir
    )

    if train_hits_path is not None:
        # Load train hits
        train_hits = pd.read_csv(train_hits_path)

        # Similarity between REAL molecules and train hits
        plot_similarity(
            smiles=data[smiles_column],
            similarity_type='tversky',
            save_dir=save_dir,
            reference_smiles=train_hits[train_hits_smiles_column],
            reference_name='Train Hits'
        )


if __name__ == '__main__':
    from tap import tapify

    tapify(plot_molecule_analysis)
