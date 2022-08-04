"""Analyze the similarities between two sets of molecules."""
from pathlib import Path
from typing import Iterable, Literal, Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from tap import Tap

from chem_utils.constants import Molecule
from chem_utils.molecular_similarities import get_similarity_function


class Args(Tap):
    data_path: Path  # Path to CSV file containing a set of molecules.
    reference_data_path: Path  # Path to CSV file containing a reference set of molecules against which molecules from data_path will be compared.
    smiles_column: str = 'smiles'  # Name of the column in data_path containing SMILES.
    reference_smiles_column: Optional[str] = None  # Name of the column in reference_data_path containing SMILES. Defaults to smiles_column.
    similarity_type: Literal['tanimoto', 'tversky']  # Type of similarity to compute.
    save_path: Path  # Path to PDF file where the results plot will be saved.
    percentile_range: tuple[float, ...] = [0.0, 100.1, 5.0]  # The [start, ]stop, [step] parameters of np.arange for the similarity percentiles to measure.

    def process_args(self) -> None:
        if self.reference_smiles_column is None:
            self.reference_smiles_column = self.smiles_column


def plot_molecular_similarities(smiles: Iterable[Molecule],
                                reference_smiles: Iterable[Molecule],
                                similarity_type: str,
                                save_path: Path,
                                percentile_range: tuple[float, ...] = (0.0, 100.1, 5.0)) -> None:
    """Plots a violin plot of similarities between two sets of molecules."""
    if not (1 <= len(percentile_range) <= 3):
        raise ValueError('Percentile range must have between 1 and 3 values (for [start, ]stop, [step]).')

    # Compute similarities
    similarity_function = get_similarity_function(similarity_type)
    pairwise_similarities = similarity_function(smiles, reference_smiles)

    # Analyze similarities at different percentiles
    percentiles = np.arange(*percentile_range)
    percentile_similarities = np.percentile(pairwise_similarities, percentiles, axis=1)

    step = percentile_range[-1] if len(percentile_range) == 3 else 1.0

    # Plot similarities at different percentiles
    plt.clf()
    plt.violinplot(percentile_similarities.transpose(), positions=percentiles, widths=step / 2, showmedians=True)
    plt.xlabel('Percentile')
    plt.ylabel(f'{similarity_type.title()} Similarity')
    plt.title(f'{similarity_type.title()} Similarity Percentiles')
    save_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(save_path, bbox_inches='tight')


def analyze_molecular_similarities(args: Args) -> None:
    """Analyze the similarities between two sets of molecules."""
    # Load data
    data = pd.read_csv(args.data_path)
    reference_data = pd.read_csv(args.reference_data_path)

    # Get SMILES
    smiles = data[args.smiles_column]
    reference_smiles = reference_data[args.reference_smiles_column]

    # Plot similarities
    plot_molecular_similarities(
        smiles=smiles,
        reference_smiles=reference_smiles,
        similarity_type=args.similarity_type,
        percentile_range=args.percentile_range,
        save_path=args.save_path
    )


if __name__ == '__main__':
    analyze_molecular_similarities(Args().parse_args())
