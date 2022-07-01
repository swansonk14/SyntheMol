"""Trains a machine learning classifier model."""
import pickle
from pathlib import Path
from typing import Literal

import pandas as pd
from sklearnex import patch_sklearn
patch_sklearn()
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.model_selection import train_test_split
from sklearn.neural_network import MLPClassifier
from tap import Tap

from chem_utils.molecular_fingerprints import compute_fingerprints


class Args(Tap):
    data_path: Path  # Path to CSV file containing data.
    save_path: Path  # Path to a PKL file where the trained model will be saved.
    model_type: Literal['rf', 'mlp']  # Type of model to train. 'rf' = random forest. 'mlp' = multilayer perceptron.
    fingerprint_type: Literal['morgan', 'rdkit']  # Type of fingerprints to use as input features.
    smiles_column: str = 'smiles'  # The name of the column containing SMILES.
    activity_column: str = 'activity'  # The name of the column containing binary activity values.

    def process_args(self) -> None:
        self.save_path.parent.mkdir(parents=True, exist_ok=True)


def train_model(args: Args) -> None:
    """Trains a machine learning classifier model."""
    # Load data
    data = pd.read_csv(args.data_path)
    print(f'Data size = {len(data):,}')

    # Get Morgan fingerprints
    fingerprints = compute_fingerprints(data[args.smiles_column], fingerprint_type=args.fingerprint_type)

    # Get activity
    activities = data[args.activity_column]

    # Split into train and test
    train_fingerprints, test_fingerprints, train_activities, test_activities = train_test_split(
        fingerprints,
        activities,
        test_size=0.1,
        random_state=0
    )
    print(f'Train size = {len(train_fingerprints):,}')
    print(f'Test size = {len(test_fingerprints):,}')

    # Build model
    if args.model_type == 'rf':
        model = RandomForestClassifier(n_jobs=-1, random_state=0)
    elif args.model_type == 'mlp':
        model = MLPClassifier(hidden_layer_sizes=(100, 100, 100), random_state=0)
    else:
        raise ValueError(f'Model type "{args.model_type}" is not supported.')

    print(model)

    # Train model
    model.fit(train_fingerprints, train_activities)

    # Evaluate model
    test_probs = model.predict_proba(test_fingerprints)[:, 1]
    print(f'Test ROC-AUC = {roc_auc_score(test_activities, test_probs):.3f}')
    print(f'Test PRC-AUC = {average_precision_score(test_activities, test_probs):.3f}')

    # Save model
    with open(args.save_path, 'wb') as f:
        pickle.dump(model, f)


if __name__ == '__main__':
    train_model(Args().parse_args())
