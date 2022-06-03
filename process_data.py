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
    # Load data files
    all_data = [pd.read_csv(data_path) for data_path in args.data_paths]
    all_smiles = []
    all_activities = []

    # Process each data file separately first
    for single_data in all_data:
        # Canonicalize SMILES
        all_smiles += [
            Chem.MolToSmiles(Chem.MolFromSmiles(smiles))
            for smiles in single_data[args.smiles_column]
        ]

        # Binarize data
        activities = single_data[args.mean_column]
        threshold = np.mean(activities) - args.num_std * np.std(activities)
        all_activities += (activities < threshold).astype(int).tolist()

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

    # Save data
    data.to_csv(args.save_path, index=False)

    # Save hits
    if args.save_hits_path is not None:
        hits = data[data[ACTIVITY_COLUMN] == 1]
        print(f'Number of hits = {len(hits):,}')
        hits.to_csv(args.save_hits_path, index=False)


if __name__ == '__main__':
    process_data(Args().parse_args())
