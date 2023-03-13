"""Compute the ROC-AUC and PRC-AUC of a set of predicted and true values."""
from pathlib import Path

import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score


def compute_auc(
        data_path: Path,
        pred_column: str,
        true_column: str
) -> None:
    """Compute the ROC-AUC and PRC-AUC of a set of predicted and true values.

    :param data_path: Path to CSV file containing predicted and true values.
    :param pred_column: Name of the column containing the predictions.
    :param true_column: Name of the column containing the true values.
    """
    # Load data
    data = pd.read_csv(data_path)

    # Compute AUCs
    roc_auc = roc_auc_score(data[true_column], data[pred_column])
    prc_auc = average_precision_score(data[true_column], data[pred_column])

    # Print AUCs
    print(f'ROC-AUC: {roc_auc:.3f}')
    print(f'PRC-AUC: {prc_auc:.3f}')


if __name__ == '__main__':
    from tap import tapify

    tapify(compute_auc)
