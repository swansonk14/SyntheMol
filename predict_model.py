"""Make predictions with a model."""
import pickle
from pathlib import Path
from typing import Literal, Optional

import numpy as np
import pandas as pd
import torch
from sklearn.ensemble import RandomForestClassifier
from sklearn.neural_network import MLPClassifier
from tap import Tap
from tqdm import tqdm

from chem_utils.molecular_fingerprints import compute_fingerprints
from chemprop.utils import load_checkpoint

from train_model import build_chemprop_data_loader, chemprop_predict


class Args(Tap):
    data_path: Path  # Path to a CSV file containing SMILES.
    model_path: Path  # Path to a directory of model checkpoints or to a specific PKL or PT file containing a trained model.
    save_path: Optional[Path] = None  # Path to a JSON file where a dictionary mapping fragments to model scores will be saved. If None, defaults to data_path.
    smiles_column: str = 'smiles'  # Name of the column containing SMILES.
    model_type: Literal['rf', 'mlp', 'chemprop']  # Type of model to train. 'rf' = random forest. 'mlp' = multilayer perceptron.
    fingerprint_type: Optional[Literal['morgan', 'rdkit']] = None  # Type of fingerprints to use as input features.
    average_preds: bool = False  # Whether to average predictions across models for an ensemble model.

    def process_args(self) -> None:
        if self.save_path is None:
            self.save_path = self.data_path

        self.save_path.parent.mkdir(parents=True, exist_ok=True)


def predict_sklearn(fingerprints: np.ndarray,
                    model_path: Path,
                    model_type: str) -> np.ndarray:
    """Make predictions with an sklearn model."""
    # Load model
    with open(model_path, 'rb') as f:
        if model_type == 'rf':
            model: RandomForestClassifier = pickle.load(f)
        elif model_type == 'mlp':
            model: MLPClassifier = pickle.load(f)
        else:
            raise ValueError(f'Model type "{model_type}" is not supported.')

    # Make predictions
    preds = model.predict_proba(fingerprints)[:, 1]

    return preds


def predict_chemprop(smiles: list[str],
                     fingerprints: Optional[np.ndarray],
                     model_path: Path) -> np.ndarray:
    """Make predictions with a chemprop model."""
    # Ensure reproducibility
    torch.manual_seed(0)
    torch.use_deterministic_algorithms(True)

    # Build data loader
    data_loader = build_chemprop_data_loader(
        smiles=smiles,
        fingerprints=fingerprints
    )

    # Load model
    model = load_checkpoint(path=model_path)

    # Make predictions
    preds = chemprop_predict(
        model=model,
        data_loader=data_loader
    )

    return preds


def predict_model(smiles: list[str],
                  fingerprints: Optional[np.ndarray],
                  model_type: str,
                  model_path: Path) -> np.ndarray:
    """Make predictions with a model."""
    # Map fragments to model scores
    if model_type == 'chemprop':
        preds = predict_chemprop(
            smiles=smiles,
            fingerprints=fingerprints,
            model_path=model_path
        )
    else:
        preds = predict_sklearn(
            fingerprints=fingerprints,
            model_path=model_path,
            model_type=model_type
        )

    return preds


def predict_ensemble(smiles: list[str],
                     fingerprint_type: Optional[str],
                     model_type: str,
                     model_path: Path,
                     average_preds: bool = True) -> np.ndarray:
    """Make predictions with an ensemble of models."""
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
        model_paths = list(model_path.glob('*.pt' if model_type == 'chemprop' else '*.pkl'))

        if len(model_paths) == 0:
            raise ValueError(f'Could not find any models in directory {model_path}.')
    else:
        model_paths = [model_path]

    preds = np.array([
        predict_model(
            smiles=smiles,
            fingerprints=fingerprints,
            model_type=model_type,
            model_path=model_path
        ) for model_path in tqdm(model_paths, desc='models')
    ])

    if average_preds:
        preds = np.mean(preds, axis=0)

    return preds


def make_predictions(args: Args) -> None:
    """Make predictions with a model and save them to a file."""
    # Load SMILES
    data = pd.read_csv(args.data_path)
    smiles = list(data[args.smiles_column])

    # Make predictions
    all_preds = predict_ensemble(
        smiles=smiles,
        fingerprint_type=args.fingerprint_type,
        model_type=args.model_type,
        model_path=args.model_path,
        average_preds=args.average_preds
    )

    # Define model string
    model_string = f'{args.model_type}{f"_{args.fingerprint_type}" if args.fingerprint_type is not None else ""}'

    if args.average_preds:
        data[f'{model_string}_ensemble_preds'] = all_preds
    else:
        for model_num, preds in enumerate(all_preds):
            data[f'{model_string}_model_{model_num}_preds'] = preds

    # Save predictions
    data.to_csv(args.save_path, index=False)


if __name__ == '__main__':
    make_predictions(Args().parse_args())
