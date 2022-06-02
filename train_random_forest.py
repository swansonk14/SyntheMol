"""Trains a random forest classifier model."""
import pickle
from pathlib import Path

import pandas as pd
# TODO: sklearn intelex
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.model_selection import train_test_split
from tap import Tap

from morgan_fingerprint import compute_morgan_fingerprints


class Args(Tap):
    data_path: Path  # Path to CSV file containing data.
    save_path: Path  # Path to a PKL file where the trained model will be saved.
    smiles_column: str = 'smiles'  # The name of the column containing SMILES.
    activity_column: str = 'activity'  # The name of the column containing binary activity values.

    def process_args(self) -> None:
        self.save_path.parent.mkdir(parents=True, exist_ok=True)


def train_random_forest(args: Args) -> None:
    """Trains a random forest classifier model."""
    # Load data
    data = pd.read_csv(args.data_path)
    print(f'Data size = {len(data):,}')

    # Get Morgan fingerprints
    fingerprints = compute_morgan_fingerprints(data[args.smiles_column])

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
    model = RandomForestClassifier(n_jobs=-1, random_state=0)

    # Train model
    model.fit(train_fingerprints, train_activities)

    # Evaluate model
    test_probs = model.predict_proba(test_fingerprints)[:, 1]
    print(f'Test ROC-AUC = {roc_auc_score(test_activities, test_probs):.2f}')
    print(f'Test PRC-AUC = {average_precision_score(test_activities, test_probs):.2f}')

    # Save model
    with open(args.save_path, 'wb') as f:
        pickle.dump(model, f)


if __name__ == '__main__':
    train_random_forest(Args().parse_args())
