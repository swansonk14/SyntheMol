"""Plot a matrix showing the generalization of a property prediction model across datasets."""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from tap import Tap

# TODO: update/fill in these numbers
DATASETS = ['AB_2560_normalized', 'Mar27_normalized', 'For_gen_AB_DRH']
METRIC_TO_MODEL_TO_MATRIX = {
    'PRC-AUC': {
        # 'random_forest': np.ndarray([[]]),
        'chemprop': np.ndarray([
            [0.579, 0.084, 0.365],
            [0.135, 0.159, 0.106],
            [0.583, 0.074, 0.306]
        ]),
        # 'chemprop_rdkit': np.ndarray([[]])
    },
    'ROC-AUC': {
        # 'random_forest': np.ndarray([[]]),
        'chemprop': np.ndarray([
            [0.822, 0.535, 0.737],
            [0.631, 0.679, 0.627],
            [0.882, 0.562, 0.764]
        ]),
        # 'chemprop_rdkit': np.ndarray([[]])
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
            # TODO: confusion matrix
            plt.savefig(args.save_dir / f'generalization_{model}_{metric}.pdf', bbox_inches='tight')


if __name__ == '__main__':
    plot_model_generalization(Args().parse_args())
