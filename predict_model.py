"""Make predictions with a model."""
import pickle
from pathlib import Path
from time import time
from typing import Literal, Optional

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.neural_network import MLPClassifier
from tap import Tap
from tqdm import tqdm

from chem_utils.molecular_fingerprints import compute_fingerprints
from chemprop.utils import load_checkpoint

from train_model import build_chemprop_data_loader, chemprop_predict


class Args(Tap):
    data_path: Path  # Path to a CSV file containing SMILES.
    model_dir: Optional[Path] = None  # Path to a directory containing PKL or PT files with trained models. (Do not use with model_path.)
    model_path: Optional[Path] = None  # Path to a PKL or PT fil containing a trained model. (Do not use with model_dir.)
    save_path: Optional[Path] = None  # Path to a JSON file where a dictionary mapping fragments to model scores will be saved. If None, defaults to data_path.
    smiles_column: str = 'smiles'  # Name of the column containing SMILES.
    model_type: Literal['rf', 'mlp', 'chemprop']  # Type of model to train. 'rf' = random forest. 'mlp' = multilayer perceptron.
    fingerprint_type: Literal['morgan', 'rdkit']  # Type of fingerprints to use as input features.

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
                     fingerprints: np.ndarray,
                     model_path: Path) -> np.ndarray:
    """Make predictions with a chemprop model."""
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
                  fingerprints: np.ndarray,
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


def make_predictions(args: Args) -> None:
    """Make predictions with a model and save them to a file."""
    # Load SMILES
    data = pd.read_csv(args.data_path)
    smiles = list(data[args.smiles_column])

    # Compute fingerprints
    fingerprints = compute_fingerprints(smiles, fingerprint_type=args.fingerprint_type)

    # Error handling on model dir/path
    if (args.model_dir is None) == (args.model_path is None):
        raise ValueError('Must provide exactly one of model_dir and model_path.')

    # Get model paths
    if args.model_path is not None:
        model_paths = [args.model_path]
    else:
        model_paths = list(args.model_dir.glob('*.pt' if args.model_type == 'chemprop' else '*.pkl'))

    # Make predictions
    start_time = time()

    for model_num, model_path in enumerate(tqdm(model_paths, desc='models')):
        preds = predict_model(
            smiles=smiles,
            fingerprints=fingerprints,
            model_type=args.model_type,
            model_path=model_path
        )

        # Add predictions to data
        data[f'{args.model_type}_model_{model_num}_preds'] = preds

    s = 's' if len(model_paths) > 1 else ''
    print(f'Total time to make predictions using {len(model_paths)} {args.model_type} model{s} '
          f'with {args.fingerprint_type} fingerprints = {time() - start_time:.2f}')

    # Save predictions
    data.to_csv(args.save_path, index=False)


if __name__ == '__main__':
    make_predictions(Args().parse_args())
