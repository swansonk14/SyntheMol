"""Compute the ROC-AUC and PRC-AUC of a set of predicted and true values."""
from pathlib import Path

import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score
from tap import Tap


class Args(Tap):
    data_path: Path  # Path to CSV file containing predicted and true values.
    pred_column: str  # Name of the column containing the predictions.
    true_column: str  # Name of the column containing the true values.


def compute_auc(args: Args) -> None:
    """Compute the ROC-AUC and PRC-AUC of a set of predicted and true values."""
    # Load data
    data = pd.read_csv(args.data_path)

    # Compute AUCs
    roc_auc = roc_auc_score(data[args.true_column], data[args.pred_column])
    prc_auc = average_precision_score(data[args.true_column], data[args.pred_column])

    # Print AUCs
    print(f'ROC-AUC: {roc_auc:.3f}')
    print(f'PRC-AUC: {prc_auc:.3f}')


if __name__ == '__main__':
    compute_auc(Args().parse_args())
