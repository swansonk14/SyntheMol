"""Contains training and predictions functions for scikit-learn models."""
import pickle
from pathlib import Path

import numpy as np
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor

from synthemol.constants import SKLEARN_MODEL_TYPES


def sklearn_load(model_path: Path) -> SKLEARN_MODEL_TYPES:
    """Loads a scikit-learn model.

    :param model_path: Path to a scikit-learn model.
    :return: A scikit-learn model.
    """
    with open(model_path, "rb") as f:
        model = pickle.load(f)

    return model


def sklearn_predict(model: SKLEARN_MODEL_TYPES, fingerprints: np.ndarray) -> np.ndarray:
    """Predicts molecular properties using a scikit-learn model.

    :param model: A scikit-learn model.
    :param fingerprints: A 2D array of molecular fingerprints (num_molecules, num_features).
    :return: A 1D array of predicted properties (num_molecules,).
    """
    if isinstance(model, RandomForestClassifier):
        preds = model.predict_proba(fingerprints)[:, 1]
    elif isinstance(model, RandomForestRegressor):
        preds = model.predict(fingerprints)
    else:
        raise ValueError(f"Model type {type(model)} is not supported.")

    return preds


def sklearn_predict_on_molecule(
    model: SKLEARN_MODEL_TYPES, fingerprint: np.ndarray | None = None
) -> np.ndarray:
    """Predicts the property of a molecule using a scikit-learn model.

    :param model: A scikit-learn model.
    :param fingerprint: A 1D array of molecular fingerprints (if applicable).
    :return: The model prediction on the molecule (1,).
    """
    return sklearn_predict(model=model, fingerprints=fingerprint.reshape(1, -1))


def sklearn_predict_on_molecule_ensemble(
    models: list[SKLEARN_MODEL_TYPES], fingerprint: np.ndarray | None = None
) -> np.ndarray:
    """Predicts the property of a molecule using an ensemble of scikit-learn models.

    :param models: An ensemble of scikit-learn models.
    :param fingerprint: A 1D array of molecular fingerprints (if applicable).
    :return: The ensemble prediction on the molecule (1,).
    """
    return np.mean(
        [
            sklearn_predict_on_molecule(model=model, fingerprint=fingerprint)
            for model in models
        ],
        axis=0,
    )
