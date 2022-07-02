"""Assess the quality and diversity of generated molecules."""
from pathlib import Path
from typing import Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from tap import Tap

from chem_utils.molecular_similarities import compute_max_similarities


class Args(Tap):
    data_path: Path  # Path to CSV file containing generated molecules.
    save_dir: Path  # Path to directory where plots and results will be saved.
    train_hits_path: Optional[Path] = None  # Path to CSV file containing hits from the training set for computing novelty.
    min_score: Optional[float] = None  # If provided, only molecules with scores >= this threshold are assessed.
    smiles_column: str = 'smiles'  # The name of the column containing SMILES in data_path.
    train_hits_smiles_column: str = 'smiles'  # The name of the column containing SMILES in train_hits_path.
    score_column: str = 'score'  # The name of the column containing scores.

    def process_args(self) -> None:
        self.save_dir.mkdir(parents=True, exist_ok=True)


SCORE_TYPES = ['score', 'model_score']
SIMILARITY_TYPES = ['tanimoto', 'tversky']


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
    # Results dictionary
    results = {}

    # Load data
    data = pd.read_csv(args.data_path)

    # Threshold molecules
    if args.min_score is not None:
        data = data[data[args.score_column] >= args.min_score]
        data.to_csv(args.save_dir / 'molecules.csv', index=False)

    # Count molecules
    results['num_molecules'] = len(data)
    print(f'Number of molecules = {results["num_molecules"]:,}')

    # Get SMILES
    smiles = data[args.smiles_column]

    for score_type in SCORE_TYPES:
        scores = data[score_type]

        # Assess score distribution
        results['score_mean'] = np.mean(scores)
        results['score_std'] = np.std(scores)
        print(f'{score_type.title()} = {results["score_mean"]:.3f} +/- {results["score_std"]:.3f}')

        # Plot score distribution
        plt.clf()
        plt.hist(scores, bins=100)
        plt.xlabel(score_type.title())
        plt.ylabel('Count')
        plt.title('Score Distribution')
        plt.tight_layout()
        plt.savefig(args.save_dir / f'{score_type}.pdf')

    # Assess diversity within generated molecules
    for similarity_type in SIMILARITY_TYPES:
        # Compute similarities within generated molecules
        generated_max_similarities = compute_max_similarities(
            similarity_type=similarity_type,
            mols=data[args.smiles_column]
        )
        results[f'generated_diversity_{similarity_type}_mean'] = np.mean(generated_max_similarities)
        results[f'generated_diversity_{similarity_type}_std'] = np.std(generated_max_similarities)
        print(f'Generated diversity {similarity_type} = '
              f'{results[f"generated_diversity_{similarity_type}_mean"]:.3f} +/- '
              f'{results[f"generated_diversity_{similarity_type}_std"]:.3f}')

        # Plot diversity distribution within generated molecules
        plt.clf()
        plt.hist(generated_max_similarities, bins=100)
        plt.xlabel(f'Maximum {similarity_type.title()} Similarity between Generated Molecules')
        plt.ylabel('Count')
        plt.title(f'Generated Maximum {similarity_type.title()} Similarity Distribution')
        plt.tight_layout()
        plt.savefig(args.save_dir / f'generated_diversity_{similarity_type}.pdf')

    # Compare to train hits
    if args.train_hits_path is not None:
        # Load train hits
        train_hits = pd.read_csv(args.train_hits_path)
        train_hits_smiles = train_hits[args.train_hits_smiles_column]

        # Assess diversity compared to train
        for similarity_type in SIMILARITY_TYPES:
            train_max_similarities = compute_max_similarities(
                similarity_type=similarity_type,
                mols=smiles,
                reference_mols=train_hits_smiles
            )
            results[f'train_diversity_{similarity_type}_mean'] = np.mean(train_max_similarities)
            results[f'train_diversity_{similarity_type}_std'] = np.std(train_max_similarities)
            print(f'Train diversity {similarity_type} = '
                  f'{results[f"train_diversity_{similarity_type}_mean"]:.3f} +/- '
                  f'{results[f"train_diversity_{similarity_type}_std"]:.3f}')

            # Plot diversity distribution compared to train
            plt.clf()
            plt.hist(train_max_similarities, bins=100)
            plt.xlabel(f'Maximum {similarity_type.title()} Similarity from Generated to Train Molecules')
            plt.ylabel('Count')
            plt.title(f'Train Maximum {similarity_type.title()} Similarity Distribution')
            plt.tight_layout()
            plt.savefig(args.save_dir / f'train_diversity_{similarity_type}.pdf')

        # Assess novelty
        results['novelty'] = compute_novelty(smiles=smiles, reference_smiles=train_hits_smiles)
        print(f'Novelty = {results["novelty"]:.3f}')

    # Save data
    evaluation = pd.DataFrame(data=[results])
    evaluation.to_csv(args.save_dir / 'evaluation.csv', index=False)


if __name__ == '__main__':
    assess_generated_molecules(Args().parse_args())
