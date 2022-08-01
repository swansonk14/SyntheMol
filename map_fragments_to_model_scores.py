"""Map fragments to model prediction scores."""
import json
from pathlib import Path
from typing import Literal

import pandas as pd
from tap import Tap

from predict_with_model import predict_with_model


class Args(Tap):
    fragment_path: Path  # Path to a CSV file containing fragments.
    model_path: Path  # Path to a PKL or PT file containing a trained model.
    save_path: Path  # Path to a JSON file where a dictionary mapping fragments to model scores will be saved.
    smiles_column: str = 'smiles'  # Name of the column containing SMILES.
    model_type: Literal['rf', 'mlp', 'chemprop']  # Type of model to train. 'rf' = random forest. 'mlp' = multilayer perceptron.
    fingerprint_type: Literal['morgan', 'rdkit']  # Type of fingerprints to use as input features.

    def process_args(self) -> None:
        self.save_path.parent.mkdir(parents=True, exist_ok=True)


def map_fragments_to_scores(args: Args) -> None:
    """Map fragments to prediction scores."""
    # Load fragments
    fragments = sorted(set(pd.read_csv(args.fragment_path)[args.smiles_column]))

    # Make predictions
    scores = predict_with_model(
        smiles=fragments,
        fingerprint_type=args.fingerprint_type,
        model_type=args.model_type,
        model_path=args.model_path
    )

    # Map fragments to predictions
    fragment_to_model_score = dict(zip(fragments, scores))

    # Save fragment to model score mapping
    with open(args.save_path, 'w') as f:
        json.dump(fragment_to_model_score, f, indent=4, sort_keys=True)


if __name__ == '__main__':
    map_fragments_to_scores(Args().parse_args())
