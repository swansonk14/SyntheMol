"""Compare model predictions on molecular fragments vs whole molecules."""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from tap import Tap


class Args(Tap):
    data_path: Path  # Path to CSV file containing prediction scores for fragments and whole molecules.
    save_dir: Path  # Path to directory where plots will be saved.
    models: list[str] = ['chemprop', 'chemprop_rdkit', 'rf_rdkit']  # Names of models to compare.

    def process_args(self) -> None:
        self.save_dir.mkdir(parents=True, exist_ok=True)


def compare_fragment_vs_whole_preds(args: Args) -> None:
    """Compare model predictions on molecular fragments vs whole molecules."""
    # Load data
    data = pd.read_csv(args.data_path)

    # For now, limit to only those molecules with two fragments
    num_fragments = 2
    data = data[data['reagent_1_3_smiles'].isna()]

    # Compare fragment vs whole molecule predictions
    for model in args.models:
        # Get fragment scores
        fragment_scores = np.sum(
            [data[f'reagent_1_{i}_{model}_ensemble_preds'].to_numpy() for i in range(1, num_fragments + 1)],
            axis=0
        )

        # Get whole molecule scores
        whole_mol_scores = data[f'{model}_ensemble_preds'].to_numpy()

        # Plot
        plt.clf()
        plt.scatter(fragment_scores, whole_mol_scores, s=2)
        plt.xlabel('Fragment score')
        plt.ylabel('Whole molecule score')
        plt.title(f'{model} fragment vs whole molecule predictions')
        plt.savefig(args.save_dir / f'{model}_fragment_vs_whole_mol_preds.pdf', bbox_inches='tight')


if __name__ == '__main__':
    compare_fragment_vs_whole_preds(Args().parse_args())
