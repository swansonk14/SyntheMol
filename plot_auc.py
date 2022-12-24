"""Plot an AUC curve given a set of predictions and true values."""
from pathlib import Path
from typing import Literal

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
from sklearn.metrics import precision_recall_curve, roc_auc_score, roc_curve, average_precision_score
from tap import Tap


class Args(Tap):
    data_dir: Path  # Path to directory containing CSV files with test set predictions and true values.
    save_path: Path  # Path to PDF or PNG file where the plot will be saved.
    model_name: str  # Name of the model whose predictions are being plotted.
    curve_type: Literal['ROC', 'PRC']  # Type of curve to plot (ROC or PR).
    activity_column: str = 'activity'  # Name of the column containing the true values.
    prediction_column: str = 'prediction'  # Name of the column containing the predictions.

    def process_args(self) -> None:
        self.save_path.parent.mkdir(parents=True, exist_ok=True)


def compute_curve(
        activities: pd.Series,
        predictions: pd.Series,
        curve_type: Literal['ROC', 'PRC']
) -> tuple[np.ndarray, np.ndarray, float]:
    """Plot an individual ROC or PRC curve."""
    if curve_type == 'ROC':
        fpr, tpr, _ = roc_curve(activities, predictions)
        x_values, y_values = fpr, tpr
        auc_score = roc_auc_score(activities, predictions)
    elif curve_type == 'PRC':
        precision, recall, _ = precision_recall_curve(activities, predictions)
        x_values, y_values = recall, precision
        auc_score = average_precision_score(activities, predictions)
    else:
        raise ValueError(f'Invalid curve type: {curve_type}')

    return x_values, y_values, auc_score


def plot_auc(args: Args) -> None:
    """Plot an AUC curve given a set of predictions and true values."""
    # Plot ROC or PRC curve for each dataset
    auc_scores, all_y_values = [], []
    base_x_values = np.linspace(0, 1, 101)
    for test_path in args.data_dir.rglob('*_test_preds.csv'):
        # Load predictions and true values
        data = pd.read_csv(test_path)

        # Plot curve
        x_values, y_values, auc_score = compute_curve(
            data[args.activity_column],
            data[args.prediction_column],
            args.curve_type
        )
        plt.plot(x_values, y_values, color='blue', alpha=0.25)

        # Save AUC score
        auc_scores.append(auc_score)

        # Interpolate y_values
        all_y_values.append(interp1d(x_values, y_values)(base_x_values))

    # Plot average curve
    mean_y_values = np.mean(all_y_values, axis=0)
    plt.plot(base_x_values, mean_y_values, color='red',
             label=rf'{args.curve_type}-AUC = ${np.mean(auc_scores):.3f} \pm {np.std(auc_scores):.3f}$')

    # Label plot
    if args.curve_type == 'ROC':
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
    elif args.curve_type == 'PRC':
        plt.xlabel('Recall')
        plt.ylabel('Precision')
    else:
        raise ValueError(f'Invalid curve type: {args.curve_type}')

    plt.legend()
    plt.title(f'{args.model_name} {args.curve_type} Curve')
    plt.savefig(args.save_path, bbox_inches='tight')


if __name__ == '__main__':
    plot_auc(Args().parse_args())
