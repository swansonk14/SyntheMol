"""Processes antibiotics data from potentially multiple files."""
from collections import defaultdict
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
from rdkit import Chem
from tap import Tap


class Args(Tap):
    data_paths: list[Path]  # A list of paths to CSV files containing data.
    save_path: Path  # A path to a CSV file where the processed data will be saved.
    save_hits_path: Optional[Path] = None  # A path to a CSV file where only the hits (actives) will be saved.
    smiles_column: str = 'SMILES'  # The name of the column containing the SMILES.
    mean_column: str = 'Mean'  # The name of the column containing the activity.
    num_std: int = 2  # The number of standard deviations to use when binarizing the data.

    def process_args(self) -> None:
        self.save_path.parent.mkdir(parents=True, exist_ok=True)

        if self.save_hits_path is not None:
            self.save_hits_path.parent.mkdir(parents=True, exist_ok=True)


ACTIVITY_COLUMN = 'activity'
SMILES_COLUMN = 'smiles'


def process_data(args: Args) -> None:
    """Processes antibiotics data from potentially multiple files."""
    # Load and process each data file
    all_smiles, all_activities = [], []

    for data_path in args.data_paths:
        # Load data
        data = pd.read_csv(data_path)
        print(data_path.stem)
        print(f'Data size = {len(data):,}')

        # Canonicalize SMILES
        all_smiles += [
            Chem.MolToSmiles(Chem.MolFromSmiles(smiles))
            for smiles in data[args.smiles_column]
        ]

        # Binarize data using mean - k * std
        activities = data[args.mean_column]
        mean_activity = np.mean(activities)
        std_activity = np.std(activities)
        threshold = mean_activity - args.num_std * std_activity
        print(f'Mean activity = {mean_activity}')
        print(f'Std activity = {std_activity}')
        print(f'Activity threshold of mean - {args.num_std} std = {threshold}')

        binary_activities = (activities < threshold).astype(int)
        print(f'Number of hits = {sum(binary_activities == 1):,}')
        print(f'Number of non-hits = {sum(binary_activities == 0):,}')
        print()

        all_activities += binary_activities.tolist()

    # Combine the data
    data = pd.DataFrame({
        SMILES_COLUMN: all_smiles,
        ACTIVITY_COLUMN: all_activities
    })
    print(f'Full data size = {len(data):,}')

    # Drop duplicates with non-conflicting values
    data.drop_duplicates(inplace=True)
    print(f'Data size after dropping non-conflicting duplicates = {len(data):,}')

    # Drop duplicates with conflicting values
    smiles_to_values = defaultdict(set)
    for smiles, value in zip(data[SMILES_COLUMN], data[ACTIVITY_COLUMN]):
        smiles_to_values[smiles].add(value)

    conflicting_smiles = {smiles for smiles, values in smiles_to_values.items() if len(values) > 1}

    data = data[~data[SMILES_COLUMN].isin(conflicting_smiles)]
    print(f'Data size after dropping conflicting duplicates = {len(data):,}')
    print()

    # Save data
    data.to_csv(args.save_path, index=False)

    print(f'Final data size = {len(data):,}')
    print(f'Number of hits = {sum(data[ACTIVITY_COLUMN] == 1):,}')
    print(f'Number of non-hits = {sum(data[ACTIVITY_COLUMN] == 0):,}')

    # Save hits
    if args.save_hits_path is not None:
        hits = data[data[ACTIVITY_COLUMN] == 1]
        hits.to_csv(args.save_hits_path, index=False)


if __name__ == '__main__':
    process_data(Args().parse_args())
