"""Evaluate the confidence of an ensemble for selecting molecules."""
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score
from tap import Tap


class Args(Tap):
    data_path: Path  # Path to CSV file containing ensemble predictions on a test set.
    save_path: Path  # Path to PDF or PNG file where results plot will be saved.
    activity_column: str = 'activity'  # Name of the column containing true activity values.
    preds_column_suffix: str = 'preds'  # Suffix of the names of the columns containing predictions.

    def process_args(self) -> None:
        self.save_path.parent.mkdir(parents=True, exist_ok=True)


def evaluate_ensemble_confidence(args: Args) -> None:
    """Evaluate the confidence of an ensemble for selecting molecules."""
    # Load data
    data = pd.read_csv(args.data_path)

    # Get preds
    preds_columns = [column for column in data.columns if column.endswith(args.preds_column_suffix)]
    all_preds = data[preds_columns].to_numpy()

    # Compute variance
    preds_std = np.std(all_preds, axis=1)

    # Use first model's predictions
    pred_column = preds_columns[0]

    # Compute performance of first model as a function of std percentile among ensemble
    percentiles = np.arange(0, 101, 10)
    percentile_stds = np.percentile(preds_std, percentiles)

    for i in range(len(percentiles) - 1):
        print(f'Std percentiles {percentiles[i]} to {percentiles[i + 1]}')

        selected_data = data[(preds_std > percentile_stds[i]) & (preds_std <= percentile_stds[i + 1])]
        print(selected_data[args.activity_column].value_counts())

        activity_classes = selected_data[args.activity_column].unique()
        if len(activity_classes) == 1:
            continue

        roc_auc = roc_auc_score(selected_data[args.activity_column], selected_data[pred_column])
        prc_auc = average_precision_score(selected_data[args.activity_column], selected_data[pred_column])

        print(f'ROC-AUC = {roc_auc:.3}')
        print(f'PRC-AUC = {prc_auc:.3}')
        print()


if __name__ == '__main__':
    evaluate_ensemble_confidence(Args().parse_args())
