"""Contains training and predictions functions for scikit-learn models."""
import pickle
from pathlib import Path

import numpy as np
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.neural_network import MLPClassifier, MLPRegressor

from synthemol.constants import DATASET_TYPES, SKLEARN_MODEL_NAME_TYPES, SKLEARN_MODEL_TYPES


def sklearn_save(
        model: SKLEARN_MODEL_TYPES,
        save_path: Path
) -> None:
    """Saves a scikit-learn model.

    :param model: A scikit-learn model.
    :param save_path: The path to save the model to.
    """
    save_path.parent.mkdir(parents=True, exist_ok=True)

    with open(save_path, 'wb') as f:
        pickle.dump(model, f)


def sklearn_build_model(
        model_type: SKLEARN_MODEL_NAME_TYPES,
        dataset_type: DATASET_TYPES
) -> SKLEARN_MODEL_TYPES:
    """Builds a scikit-learn model.

    :param model_type: The type of model (random forest or multilayer perceptron).
    :param dataset_type: The type of dataset (classification or regression).
    :return: A scikit-learn model.
    """
    if model_type == 'random_forest':
        if dataset_type == 'classification':
            random_forest_class = RandomForestClassifier
        elif dataset_type == 'regression':
            random_forest_class = RandomForestRegressor
        else:
            raise ValueError(f'Dataset type "{dataset_type}" is not supported.')

        model = random_forest_class(n_jobs=-1, random_state=0)
    elif model_type == 'mlp':
        if dataset_type == 'classification':
            mlp_class = MLPClassifier
        elif dataset_type == 'regression':
            mlp_class = MLPRegressor
        else:
            raise ValueError(f'Dataset type "{dataset_type}" is not supported.')

        model = mlp_class(hidden_layer_sizes=(100, 100, 100), random_state=0)
    else:
        raise ValueError(f'Model type "{model_type}" is not supported.')

    return model


def sklearn_train(
        model_type: SKLEARN_MODEL_NAME_TYPES,
        dataset_type: DATASET_TYPES,
        fingerprints: np.ndarray,
        properties: np.ndarray,
        save_path: Path
) -> SKLEARN_MODEL_TYPES:
    """Trains a scikit-learn model.

    :param model_type: The type of model (random forest or multilayer perceptron).
    :param dataset_type: The type of dataset (classification or regression).
    :param fingerprints: A 2D array of molecular fingerprints (num_molecules, num_features).
    :param properties: A 1D array of molecular properties (num_molecules,).
    :param save_path: The path to save the model to.
    :return: A trained scikit-learn model.
    """
    # Build model
    model = sklearn_build_model(
        model_type=model_type,
        dataset_type=dataset_type
    )

    # Train model
    model.fit(fingerprints, properties)

    # Save model
    sklearn_save(
        model=model,
        save_path=save_path
    )

    return model
