"""Plot a matrix showing the generalization of a property prediction model across datasets."""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from tap import Tap

# TODO: update/fill in these numbers
DATASETS = ['AB_2560_normalized', 'Mar27_normalized', 'For_gen_AB_DRH']
METRIC_TO_MODEL_TO_MATRIX = {
    'PRC-AUC': {
        'random_forest': np.ndarray([
            [],
            [],
            []
        ]),
        'chemprop': np.ndarray([
            [],
            [],
            []
        ]),
        'chemprop_rdkit': np.ndarray([
            [],
            [],
            []
        ])
    },
    'ROC-AUC': {
        'random_forest': np.ndarray([
            [],
            [],
            []
        ]),
        'chemprop': np.ndarray([
            [],
            [],
            []
        ]),
        'chemprop_rdkit': np.ndarray([
            [],
            [],
            []
        ])
    }
}


class Args(Tap):
    save_dir: Path  # Path to directory where plots will be saved.

    def process_args(self) -> None:
        self.save_dir.mkdir(parents=True, exist_ok=True)


def plot_model_generalization(args: Args) -> None:
    """Plot a matrix showing the generalization of a property prediction model across datasets."""
    for metric, model_to_matrix in METRIC_TO_MODEL_TO_MATRIX.items():
        for model, matrix in model_to_matrix.items():
            plt.clf()
            plt.imshow(matrix, cmap='Blues', vmin=0, vmax=1)
            for i in range(matrix.shape[0]):
                for j in range(matrix.shape[1]):
                    plt.text(j, i, f'{matrix[i, j]:.3f}', ha='center', va='center', color='black')
            plt.xticks(np.arange(len(DATASETS)), DATASETS, rotation=45)
            plt.yticks(np.arange(len(DATASETS)), DATASETS)
            plt.xlabel('Test Dataset')
            plt.ylabel('Train Dataset')
            plt.title(f'{model} {metric} Generalization')
            plt.colorbar()
            plt.savefig(args.save_dir / f'generalization_{model}_{metric}.pdf', bbox_inches='tight')


if __name__ == '__main__':
    plot_model_generalization(Args().parse_args())
