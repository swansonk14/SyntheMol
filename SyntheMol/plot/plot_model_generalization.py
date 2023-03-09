"""Plot a matrix showing the generalization of a property prediction model across datasets."""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Train in rows, predict in columns
DATASETS = ['library_1', 'library_2', 'library_3']
METRIC_TO_MODEL_TO_MATRIX = {
    'PRC-AUC': {
        'random_forest': np.array([
            [0.546, 0.361, 0.082],
            [0.571, 0.389, 0.074],
            [0.175, 0.125, 0.181]
        ]),
        'chemprop': np.array([
            [0.528, 0.298, 0.055],
            [0.567, 0.328, 0.051],
            [0.098, 0.067, 0.069]
        ]),
        'chemprop_rdkit': np.array([
            [0.52, 0.365, 0.07],
            [0.633, 0.382, 0.058],
            [0.105, 0.072, 0.118]
        ])
    },
    'ROC-AUC': {
        'random_forest': np.array([
            [0.874, 0.772, 0.59],
            [0.888, 0.789, 0.605],
            [0.747, 0.705, 0.765]
        ]),
        'chemprop': np.array([
            [0.865, 0.752, 0.534],
            [0.855, 0.793, 0.624],
            [0.605, 0.6, 0.661]
        ]),
        'chemprop_rdkit': np.array([
            [0.873, 0.77, 0.574],
            [0.899, 0.797, 0.69],
            [0.606, 0.609, 0.702]
        ])
    }
}


def plot_model_generalization(
        save_dir: Path
) -> None:
    """Plot a matrix showing the generalization of a property prediction model across datasets."""
    save_dir.mkdir(parents=True, exist_ok=True)

    for metric, model_to_matrix in METRIC_TO_MODEL_TO_MATRIX.items():
        for model, matrix in model_to_matrix.items():
            plt.clf()

            # Plot AUC confusion matrix
            plt.imshow(matrix, cmap='Blues', vmin=0.5 if metric == 'ROC-AUC' else 0, vmax=1)
            for i in range(matrix.shape[0]):
                for j in range(matrix.shape[1]):
                    plt.text(j, i, f'{matrix[i, j]:.3f}', ha='center', va='center', color='black')

            # Label plot
            plt.xticks(np.arange(len(DATASETS)), DATASETS, rotation=45)
            plt.yticks(np.arange(len(DATASETS)), DATASETS)
            plt.xlabel('Test Dataset')
            plt.ylabel('Train Dataset')
            plt.title(f'{model} {metric} Generalization')
            plt.colorbar()
            plt.savefig(save_dir / f'generalization_{model}_{metric}.pdf', bbox_inches='tight')

            # Save matrix data
            fig_data = pd.DataFrame(
                data=matrix,
                index=[f'Train {dataset}' for dataset in DATASETS],
                columns=[f'Predict {dataset}' for dataset in DATASETS]
            )
            fig_data.to_csv(save_dir / f'generalization_{model}_{metric}.csv')


if __name__ == '__main__':
    from tap import tapify

    tapify(plot_model_generalization)
