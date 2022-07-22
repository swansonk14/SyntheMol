"""Map fragments to prediction scores."""
import json
import pickle
from pathlib import Path
from typing import Literal

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
    fragment_path: Path  # Path to a CSV file containing fragments.
    model_path: Path  # Path to a PKL or PT file containing a trained model.
    save_path: Path  # Path to a JSON file where a dictionary mapping fragments to scores will be saved.
    smiles_column: str = 'smiles'  # Name of the column containing SMILES.
    model_type: Literal['rf', 'mlp', 'chemprop']  # Type of model to train. 'rf' = random forest. 'mlp' = multilayer perceptron.
    fingerprint_type: Literal['morgan', 'rdkit']  # Type of fingerprints to use as input features.

    def process_args(self) -> None:
        self.save_path.parent.mkdir(parents=True, exist_ok=True)


def map_fragments_to_scores_sklearn(fragments: list[str],
                                    fingerprints: np.ndarray,
                                    model_path: Path,
                                    model_type: str) -> dict[str, float]:
    """Map fragments to prediction scores using a scikit-learn model."""
    # Load model
    with open(model_path, 'rb') as f:
        if model_type == 'rf':
            model: RandomForestClassifier = pickle.load(f)
        elif model_type == 'mlp':
            model: MLPClassifier = pickle.load(f)
        else:
            raise ValueError(f'Model type "{model_type}" is not supported.')

    # Make predictions
    scores = model.predict_proba(fingerprints)[:, 1]

    # Map fragments to predictions
    fragment_to_score = dict(zip(fragments, scores))

    return fragment_to_score


def map_fragments_to_scores_chemprop(fragments: list[str],
                                     fingerprints: np.ndarray,
                                     model_path: Path) -> dict[str, float]:
    """Map fragments to prediction scores using a chemprop model."""
    # Build data loader
    data_loader = MoleculeDataLoader(
        dataset=MoleculeDataset([
            MoleculeDatapoint(
                smiles=[fragment],
                features=fingerprint,
            ) for fragment, fingerprint in zip(fragments, fingerprints)
        ]),
        num_workers=0,
        shuffle=False
    )

    # Load model
    model = load_checkpoint(path=model_path)

    # Make predictions
    scores = predict(
        model=model,
        data_loader=data_loader
    )
    breakpoint()
    scores = [score[0] for score in scores]

    # Map fragments to predictions
    fragment_to_score = dict(zip(fragments, scores))

    return fragment_to_score


def map_fragments_to_scores(args: Args) -> None:
    """Map fragments to prediction scores."""
    # Load fragments
    fragments = sorted(set(pd.read_csv(args.fragment_path)[args.smiles_column]))

    # Compute fingerprints
    fingerprints = compute_fingerprints(fragments, fingerprint_type=args.fingerprint_type)

    # Map fragments to predictions
    if args.model_type == 'chemprop':
        fragment_to_score = map_fragments_to_scores_chemprop(
            fragments=fragments,
            fingerprints=fingerprints,
            model_path=args.model_path
        )
    else:
        fragment_to_score = map_fragments_to_scores_sklearn(
            fragments=fragments,
            fingerprints=fingerprints,
            model_path=args.model_path,
            model_type=args.model_type
        )

    # Save data
    with open(args.save_path, 'w') as f:
        json.dump(fragment_to_score, f, indent=4, sort_keys=True)


if __name__ == '__main__':
    map_fragments_to_scores(Args().parse_args())
