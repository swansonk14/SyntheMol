"""Analyze the similarities between two sets of molecules."""
from pathlib import Path
from typing import Literal, Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from tap import Tap

from chem_utils.molecular_similarities import get_similarity_function


class Args(Tap):
    data_path: Path  # Path to CSV file containing a set of molecules.
    reference_data_path: Path  # Path to CSV file containing a reference set of molecules against which molecules from data_path will be compared.
    smiles_column: str = 'smiles'  # Name of the column in data_path containing SMILES.
    reference_smiles_column: Optional[str] = None  # Name of the column in reference_data_path containing SMILES. Defaults to smiles_column.
    similarity_type: Literal['tanimoto', 'tversky']  # Type of similarity to compute.
    save_path: Path  # Path to PDF file where the results plot will be saved.

    def process_args(self) -> None:
        if self.reference_smiles_column is None:
            self.reference_smiles_column = self.smiles_column

        self.save_path.parent.mkdir(parents=True, exist_ok=True)


def analyze_molecular_similarities(args: Args) -> None:
    """Analyze the similarities between two sets of molecules."""
    # Load data
    data = pd.read_csv(args.data_path)
    reference_data = pd.read_csv(args.reference_data_path)

    # Get SMILES
    smiles = data[args.smiles_column].iloc[:500]
    reference_smiles = reference_data[args.reference_smiles_column].iloc[:400]

    # Compute similarities
    similarity_function = get_similarity_function(args.similarity_type)
    pairwise_similarities = similarity_function(smiles, reference_smiles)

    # Analyze similarities at different percentiles
    percentiles = np.arange(0, 101, 5)
    percentile_similarities = np.percentile(pairwise_similarities, percentiles, axis=1)

    mean_percentile_similarities = np.mean(percentile_similarities, axis=1)
    std_percentile_similarities = np.std(percentile_similarities, axis=1)

    # Plot similarities at different percentiles
    plt.clf()
    plt.bar(percentiles, mean_percentile_similarities, yerr=std_percentile_similarities, width=3, capsize=3)
    plt.xlabel('Percentile')
    plt.ylabel(f'{args.similarity_type.title()} Similarity')
    plt.title(f'{args.similarity_type.title()} Similarity Percentiles')
    plt.savefig(args.save_path, bbox_inches='tight')


if __name__ == '__main__':
    analyze_molecular_similarities(Args().parse_args())
