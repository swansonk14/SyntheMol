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

from chem_utils.molecular_fingerprints import compute_fingerprints
from chemprop.data import MoleculeDataLoader, MoleculeDatapoint, MoleculeDataset
from chemprop.train import predict
from chemprop.utils import load_checkpoint


class Args(Tap):
    data_path: Path  # Path to a CSV file containing SMILES.
    model_path: Path  # Path to a PKL or PT file containing a trained model.
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
                    model_type: str) -> list[float]:
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
    preds = model.predict_proba(fingerprints)[:, 1].tolist()

    return preds


def predict_chemprop(smiles: list[str],
                     fingerprints: np.ndarray,
                     model_path: Path) -> list[float]:
    """Make predictions with a chemprop model."""
    # Build data loader
    data_loader = MoleculeDataLoader(
        dataset=MoleculeDataset([
            MoleculeDatapoint(
                smiles=[fragment],
                features=fingerprint,
            ) for fragment, fingerprint in zip(smiles, fingerprints)
        ]),
        num_workers=0,
        shuffle=False
    )

    # Load model
    model = load_checkpoint(path=model_path)

    # Make predictions
    preds = predict(
        model=model,
        data_loader=data_loader
    )
    preds = [pred[0] for pred in preds]

    return preds


def predict_with_model(smiles: list[str],
                       fingerprint_type: str,
                       model_type: str,
                       model_path: Path) -> list[float]:
    """Make predictions with a model."""
    # Compute fingerprints
    fingerprints = compute_fingerprints(smiles, fingerprint_type=fingerprint_type)

    # Map fragments to model scores
    start_time = time()
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
    print(f'Total time to make predictions using {model_type} model '
          f'with {fingerprint_type} fingerprints = {time() - start_time:.2f}')

    return preds


def make_predictions(args: Args) -> None:
    """Make predictions with a model and save them to a file."""
    # Load SMILES
    data = pd.read_csv(args.data_path)
    smiles = list(data[args.smiles_column])

    # Make predictions
    preds = predict_with_model(
        smiles=smiles,
        fingerprint_type=args.fingerprint_type,
        model_type=args.model_type,
        model_path=args.model_path
    )

    # Add predictions to data
    data[f'{args.model_type}_preds'] = preds

    # Save predictions
    data.to_csv(args.save_path, index=False)


if __name__ == '__main__':
    make_predictions(Args().parse_args())
