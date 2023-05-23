"""Make predictions with a model or ensemble of models and save them to a file."""
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from chemfunc import compute_fingerprints
from tqdm import tqdm

from SyntheMol.constants import FINGERPRINT_TYPES, MODEL_TYPES, SMILES_COL
from SyntheMol.models.chemprop import chemprop_load, chemprop_predict
from SyntheMol.models.sklearn import sklearn_load, sklearn_predict


def predict(
        data_path: Path,
        model_path: Path,
        model_type: MODEL_TYPES,
        save_path: Path | None = None,
        smiles_column: str = SMILES_COL,
        preds_column_prefix: str | None = None,
        fingerprint_type: FINGERPRINT_TYPES | None = None,
        average_preds: bool = False
) -> None:
    """Make predictions with a model or ensemble of models and save them to a file.

    :param data_path: Path to a CSV file containing SMILES.
    :param model_path: Path to a directory of model checkpoints or to a specific PKL or PT file containing a trained model.
    :param model_type: Type of model to use.
    :param save_path: Path to a CSV file where model predictions will be saved. If None, defaults to data_path.
    :param smiles_column: Name of the column containing SMILES.
    :param preds_column_prefix: Prefix for the column containing model predictions.
    :param fingerprint_type: Type of fingerprints to use as input features.
    :param average_preds: Whether to average predictions across models for an ensemble model.
    """
    # Load SMILES
    data = pd.read_csv(data_path)
    smiles = list(data[smiles_column])

    # Check compatibility of model and fingerprint type
    if model_type != 'chemprop' and fingerprint_type is None:
        raise ValueError('Must define fingerprint_type if using sklearn model.')

    # Compute fingerprints
    if fingerprint_type is not None:
        fingerprints = compute_fingerprints(smiles, fingerprint_type=fingerprint_type)
    else:
        fingerprints = None

    # Get model paths
    if model_path.is_dir():
        model_paths = list(model_path.glob('**/*.pt' if model_type == 'chemprop' else '**/*.pkl'))

        if len(model_paths) == 0:
            raise ValueError(f'Could not find any models in directory {model_path}.')
    else:
        model_paths = [model_path]

    # Load models
    if model_type == 'chemprop':
        # Ensure reproducibility
        torch.manual_seed(0)
        torch.use_deterministic_algorithms(True)

        models = [chemprop_load(model_path=model_path) for model_path in model_paths]
    else:
        models = [sklearn_load(model_path=model_path) for model_path in model_paths]

    # Make predictions
    if model_type == 'chemprop':
        preds = np.array([
            chemprop_predict(
                model=model,
                smiles=smiles,
                fingerprints=fingerprints,
            ) for model in tqdm(models, desc='models')
        ])
    else:
        preds = np.array([
            sklearn_predict(
                model=model,
                fingerprints=fingerprints,
            ) for model in tqdm(models, desc='models')
        ])

    if average_preds:
        preds = np.mean(preds, axis=0)

    # Define model string
    model_string = f'{model_type}{f"_{fingerprint_type}" if fingerprint_type is not None else ""}'
    preds_string = f'{f"{preds_column_prefix}_" if preds_column_prefix is not None else ""}{model_string}'

    if average_preds:
        data[f'{preds_string}_ensemble_preds'] = preds
    else:
        for model_num, model_preds in enumerate(preds):
            data[f'{preds_string}_model_{model_num}_preds'] = model_preds

    # Save predictions
    if save_path is None:
        save_path = data_path

    save_path.parent.mkdir(parents=True, exist_ok=True)
    data.to_csv(save_path, index=False)


if __name__ == '__main__':
    from tap import tapify

    tapify(predict)
