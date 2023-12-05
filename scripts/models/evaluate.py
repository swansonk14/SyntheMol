"""Evaluate predictions from a model."""
import numpy as np
from sklearn.metrics import (
    average_precision_score,
    mean_absolute_error,
    r2_score,
    roc_auc_score,
)

from synthemol.constants import DATASET_TYPES


def evaluate_classification(true: np.ndarray, preds: np.ndarray) -> dict[str, float]:
    """Evaluates classification predictions.

    :param true: A 1D array of true values (num_molecules,).
    :param preds: A 1D array of predicted values (num_molecules,).
    :return: A dictionary of scores.
    """
    return {
        "roc_auc": roc_auc_score(true, preds),
        "prc_auc": average_precision_score(true, preds),
    }


def evaluate_regression(true: np.ndarray, preds: np.ndarray) -> dict[str, float]:
    """Evaluates regression predictions.

    :param true: A 1D array of true values (num_molecules,).
    :param preds: A 1D array of predicted values (num_molecules,).
    :return: A dictionary of scores.
    """
    return {
        "mae": mean_absolute_error(true, preds),
        "r2": r2_score(true, preds),
    }


def evaluate(
    true: np.ndarray, preds: np.ndarray, dataset_type: DATASET_TYPES
) -> dict[str, float]:
    """Evaluates predictions.

    :param true: A 1D array of true values (num_molecules,).
    :param preds: A 1D array of predicted values (num_molecules,).
    :param dataset_type: The type of dataset (classification or regression).
    :return: A dictionary of scores.
    """
    if dataset_type == "classification":
        return evaluate_classification(true=true, preds=preds)
    elif dataset_type == "regression":
        return evaluate_regression(true=true, preds=preds)
    else:
        raise ValueError(f'Dataset type "{dataset_type}" is not supported.')
