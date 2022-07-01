"""Map fragments to prediction scores."""
import json
import pickle
from pathlib import Path
from typing import Literal

import pandas as pd
from sklearnex import patch_sklearn
patch_sklearn()
from sklearn.ensemble import RandomForestClassifier
from tap import Tap

from chem_utils.molecular_fingerprints import compute_fingerprints


class Args(Tap):
    fragment_path: Path  # Path to a CSV file containing fragments.
    model_path: Path  # Path to a PKL file containing a trained RandomForestClassifier model.
    save_path: Path  # Path to a JSON file where a dictionary mapping fragments to scores will be saved.
    smiles_column: str = 'smiles'  # Name of the column containing SMILES.
    fingerprint_type: Literal['morgan', 'rdkit']  # Type of fingerprints to use as input features.

    def process_args(self) -> None:
        self.save_path.parent.mkdir(parents=True, exist_ok=True)


def map_fragments_to_scores(args: Args) -> None:
    """Map fragments to prediction scores."""
    # Load fragments
    fragments = sorted(set(pd.read_csv(args.fragment_path)[args.smiles_column]))

    # Load model
    with open(args.model_path, 'rb') as f:
        model: RandomForestClassifier = pickle.load(f)

    # Compute Morgan fingerprints
    fingerprints = compute_fingerprints(fragments, fingerprint_type=args.fingerprint_type)

    # Make predictions
    scores = model.predict_proba(fingerprints)[:, 1]

    # Map fragments to predictions
    fragment_to_score = dict(zip(fragments, scores))

    # Save data
    with open(args.save_path, 'w') as f:
        json.dump(fragment_to_score, f, indent=4, sort_keys=True)


if __name__ == '__main__':
    map_fragments_to_scores(Args().parse_args())
