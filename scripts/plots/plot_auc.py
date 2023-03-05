"""Plot ROC or precision recall curves and compute AUCs given a set of predictions and true values."""
from pathlib import Path
from typing import Literal

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
from sklearn.metrics import precision_recall_curve, roc_auc_score, roc_curve, average_precision_score
from tap import tapify


def compute_curve(
        activities: pd.Series,
        predictions: pd.Series,
        curve_type: Literal['ROC', 'PRC']
) -> tuple[np.ndarray, np.ndarray, float]:
    """Plot an individual ROC or precision recall curve and compute AUC.

    :param activities: True binary activities.
    :param predictions: Predicted probabilities of activity.
    :param curve_type: Type of curve to plot (ROC or PRC).
    :return: A tuple of x values, y values, and AUC.
    """
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


def plot_auc(
        data_dir: Path,
        save_dir: Path,
        model_name: str,
        curve_type: Literal['ROC', 'PRC'],
        activity_column: str = 'activity',
        prediction_column: str = 'prediction'
) -> None:
    """Plot ROC or precision recall curves and compute AUCs given a set of predictions and true values.

    :param data_dir: Path to directory containing CSV files with test set predictions and true values.
    :param save_dir: Path to a directory where the plot will be saved.
    :param model_name: Name of the model whose predictions are being plotted.
    :param curve_type: Type of curve to plot (ROC or PRC).
    :param activity_column: Name of the column containing the true values.
    :param prediction_column: Name of the column containing the predictions.
    """
    # Label setup
    if curve_type == 'ROC':
        x_label = 'False Positive Rate'
        y_label = 'True Positive Rate'
    elif curve_type == 'PRC':
        x_label = 'Recall'
        y_label = 'Precision'
    else:
        raise ValueError(f'Invalid curve type: {curve_type}')

    # Plot ROC or PRC curve for each dataset
    coordinates = {}
    auc_scores, all_y_values = [], []
    base_x_values = np.linspace(0, 1, 101)
    for i, test_path in enumerate(data_dir.rglob('*_test_preds.csv')):
        # Load predictions and true values
        data = pd.read_csv(test_path)

        # Plot curve
        x_values, y_values, auc_score = compute_curve(
            data[activity_column],
            data[prediction_column],
            curve_type
        )
        plt.plot(x_values, y_values, color='blue', alpha=0.25)

        # Save curve data
        coordinates[f'{x_label} {i}'] = x_values
        coordinates[f'{y_label} {i}'] = y_values

        # Save AUC score
        auc_scores.append(auc_score)

        # Interpolate y_values
        all_y_values.append(interp1d(x_values, y_values)(base_x_values))

    # Plot average curve
    mean_y_values = np.mean(all_y_values, axis=0)
    plt.plot(base_x_values, mean_y_values, color='red',
             label=rf'{curve_type}-AUC = ${np.mean(auc_scores):.3f} \pm {np.std(auc_scores):.3f}$')

    coordinates[f'{x_label} mean'] = base_x_values
    coordinates[f'{y_label} mean'] = mean_y_values

    # Label plot
    plt.xlabel(x_label)
    plt.ylabel(y_label)

    plt.legend()
    plt.title(f'{model_name} {curve_type} Curve')

    # Save plot
    save_dir.mkdir(parents=True, exist_ok=True)
    save_name = save_dir / f'{model_name.lower().replace(" ", "_")}_{curve_type.lower()}'
    plt.savefig(f'{save_name}.pdf', bbox_inches='tight')

    # Save data
    max_len = max(len(values) for values in coordinates.values())
    fig_data = pd.DataFrame({
        key: np.pad(values, (0, max_len - len(values)), constant_values=np.nan)
        for key, values in coordinates.items()
    })
    fig_data.to_csv(f'{save_name}.csv', index=False)


if __name__ == '__main__':
    tapify(plot_auc)
