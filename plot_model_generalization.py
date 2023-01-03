"""Plot a matrix showing the generalization of a property prediction model across datasets."""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from tap import Tap

# Train in rows, predict in columns
DATASETS = ['AB_2560_normalized', 'Mar27_normalized', 'For_gen_AB_DRH']
METRIC_TO_MODEL_TO_MATRIX = {
    'PRC-AUC': {
        'random_forest': np.array([
            [0.546, 0.082, 0.361],
            [0.175, 0.181, 0.125],
            [0.571, 0.074, 0.389]
        ]),
        'chemprop': np.array([
            [0.528, 0.055, 0.298],
            [0.098, 0.069, 0.067],
            [0.567, 0.051, 0.328]
        ]),
        'chemprop_rdkit': np.array([
            [0.520, 0.070, 0.365],
            [0.105, 0.118, 0.072],
            [0.633, 0.058, 0.382]
        ])
    },
    'ROC-AUC': {
        'random_forest': np.array([
            [0.874, 0.590, 0.772],
            [0.747, 0.765, 0.705],
            [0.888, 0.605, 0.789]
        ]),
        'chemprop': np.array([
            [0.865, 0.534, 0.752],
            [0.605, 0.661, 0.600],
            [0.855, 0.624, 0.793]
        ]),
        'chemprop_rdkit': np.array([
            [0.873, 0.574, 0.770],
            [0.606, 0.702, 0.609],
            [0.899, 0.690, 0.797]
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
            plt.imshow(matrix, cmap='Blues', vmin=0.5 if metric == 'ROC-AUC' else 0, vmax=1)
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
