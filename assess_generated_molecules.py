"""Assess the quality and diversity of generated molecules."""
from pathlib import Path
from typing import Optional, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from rdkit import Chem
# TODO: sklearn intelex
from sklearn.metrics import pairwise_distances
from tap import Tap

from morgan_fingerprint import compute_morgan_fingerprints


class Args(Tap):
    data_path: Path  # Path to CSV file containing generated molecules.
    save_dir: Path  # Path to directory where plots and results will be saved.
    train_hits_path: Optional[Path] = None  # Path to CSV file containing hits from the training set for computing novelty.
    smiles_column: str = 'smiles'  # The name of the column containing SMILES in data_path.
    train_hits_smiles_column: str = 'smiles'  # The name of the column containing SMILES in train_hits_path.
    score_column: str = 'score'  # The name of the column containing scores.

    def process_args(self) -> None:
        self.save_dir.mkdir(parents=True, exist_ok=True)


# TODO: load this from chem_utils
def compute_pairwise_tanimoto_distances(mols_1: list[Union[str, Chem.Mol]],
                                        mols_2: Optional[list[Union[str, Chem.Mol]]] = None) -> np.ndarray:
    """
    Computes pairwise Tanimoto distances between the molecules in mols_1 and mols_2.

    :param mols_1: A list of molecules, either SMILES strings or RDKit molecules.
    :param mols_2: A list of molecules, either SMILES strings or RDKit molecules.
                   If None, copies mols_1 list.
    :return: A 2D numpy array of pairwise distances.
    """
    # Compute Morgan fingerprints
    fps_1 = np.array(compute_morgan_fingerprints(mols_1), dtype=bool)
    fps_2 = np.array(compute_morgan_fingerprints(mols_2), dtype=bool) if mols_2 is not None else fps_1

    # Compute pairwise distances
    tanimoto_distances = pairwise_distances(fps_1, fps_2, metric='jaccard', n_jobs=-1)

    return tanimoto_distances


# TODO: load this from chem_utils
def compute_min_tanimoto_distances(smiles: list[str]) -> np.ndarray:
    """Computes the minimum Tanimoto distance between each molecule and every other molecule.

    Note: Does NOT remove duplicate SMILES before computing pairwise distances.

    :param smiles: A list of SMILES.
    :return: The minimum Tanimoto distance between each molecule and every other molecule.
    """
    # Compute pairwise Tanimoto distances
    pairwise_tanimoto_distances = compute_pairwise_tanimoto_distances(mols_1=smiles)

    # Compute average minimum Tanimoto distance between each molecule and the other molecules
    np.fill_diagonal(pairwise_tanimoto_distances, np.inf)
    min_tanimoto_distances = np.min(pairwise_tanimoto_distances, axis=1)

    return min_tanimoto_distances


def compute_novelty(smiles: list[str], reference_smiles: set[str]) -> float:
    """Computes the novelty of a set of molecules compared to a reference set (i.e., the proportion of new molecules).

    :param smiles: A list of new SMILES.
    :param reference_smiles: A reference set of SMILES against which to compute novelty.
    :return: The novelty (proportion of new molecules) of the set of molecules.
    """
    novelty = float(np.mean([smile not in reference_smiles for smile in smiles]))

    return novelty


def assess_generated_molecules(args: Args) -> None:
    """Assess the quality and diversity of generated molecules."""
    # Load data
    data = pd.read_csv(args.data_path)
    print(f'Number of molecules = {len(data):,}')

    # Get SMILES and scores
    smiles = data[args.smiles_column]
    scores = data[args.score_column]

    # Assess score distribution
    score_mean = np.mean(scores)
    score_std = np.std(scores)
    print(f'Scores = {score_mean:.3f} +/- {score_std:.3f}')

    # Plot score distribution
    plt.clf()
    plt.hist(scores, bins=100)
    plt.xlabel('Score')
    plt.ylabel('Count')
    plt.title('Score Distribution')
    plt.tight_layout()
    plt.savefig(args.save_dir / 'scores.pdf')

    # Assess diversity
    min_tanimoto_distances = compute_min_tanimoto_distances(smiles)
    diversity_mean = np.mean(min_tanimoto_distances)
    diversity_std = np.std(min_tanimoto_distances)
    print(f'Diversity = {diversity_mean:.3f} +/- {diversity_std:.3f}')

    # Plot diversity distribution
    plt.clf()
    plt.hist(min_tanimoto_distances, bins=100)
    plt.xlabel('Minimum Tanimoto Distance')
    plt.ylabel('Count')
    plt.title('Minimum Tanimoto Distance Distribution')
    plt.tight_layout()
    plt.savefig(args.save_dir / 'diversity.pdf')

    # Assess novelty
    if args.train_hits_path is not None:
        train_hits = pd.read_csv(args.train_hits_path)
        train_hits_smiles = set(train_hits[args.train_hits_smiles_column])
        novelty = compute_novelty(smiles=smiles, reference_smiles=train_hits_smiles)
        print(f'Novelty = {novelty:.3f}')
    else:
        novelty = None

    # Save data
    evaluation = pd.DataFrame(data=[{
        'score_mean': score_mean,
        'score_std': score_std,
        'diversity_mean': diversity_mean,
        'diversity_std': diversity_std,
        'novelty': novelty
    }])
    evaluation.to_csv(args.save_dir / 'evaluation.csv', index=False)


if __name__ == '__main__':
    assess_generated_molecules(Args().parse_args())
