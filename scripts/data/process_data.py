"""Process regression data from potentially multiple files, including binarization and deduplication."""
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem

from synthemol.constants import ACTIVITY_COL, SMILES_COL


def process_data(
        data_paths: list[Path],
        save_path: Path,
        save_hits_path: Path | None = None,
        smiles_column: str = SMILES_COL,
        mean_column: str = 'mean',
        activity_column: str = ACTIVITY_COL,
        num_std: int = 2
) -> None:
    """Process regression data from potentially multiple files, including binarization and deduplication.

    :param data_paths: A list of paths to CSV files containing data.
    :param save_path: A path to a CSV file where the processed data will be saved.
    :param save_hits_path: A path to a CSV file where only the hits (actives) will be saved.
    :param smiles_column: The name of the column containing the SMILES.
    :param mean_column: The name of the column containing the mean activity.
    :param activity_column: The name of the column where the binary activity will be stored.
    :param num_std: The number of standard deviations to use when binarizing the data.
    """
    # Load and process each data file
    all_smiles, all_activities = [], []

    for data_path in data_paths:
        # Load data
        data = pd.read_csv(data_path)
        print(data_path.stem)
        print(f'Data size = {len(data):,}')

        # Canonicalize SMILES
        all_smiles += [
            Chem.MolToSmiles(Chem.MolFromSmiles(smiles))
            for smiles in data[smiles_column]
        ]

        # Binarize data using mean - k * std
        activities = data[mean_column]
        mean_activity = np.mean(activities)
        std_activity = np.std(activities)
        threshold = mean_activity - num_std * std_activity
        print(f'Mean activity = {mean_activity}')
        print(f'Std activity = {std_activity}')
        print(f'Activity threshold of mean - {num_std} std = {threshold}')

        binary_activities = (activities < threshold).astype(int)
        print(f'Number of hits = {sum(binary_activities == 1):,}')
        print(f'Number of non-hits = {sum(binary_activities == 0):,}')
        print()

        all_activities += binary_activities.tolist()

    # Combine the data
    data = pd.DataFrame({
        smiles_column: all_smiles,
        activity_column: all_activities
    })
    print(f'Full data size = {len(data):,}')

    # Drop duplicates with non-conflicting values
    data.drop_duplicates(inplace=True)
    print(f'Data size after dropping non-conflicting duplicates = {len(data):,}')

    # Drop duplicates with conflicting values
    smiles_to_values = defaultdict(set)
    for smiles, value in zip(data[smiles_column], data[activity_column]):
        smiles_to_values[smiles].add(value)

    conflicting_smiles = {smiles for smiles, values in smiles_to_values.items() if len(values) > 1}

    data = data[~data[smiles_column].isin(conflicting_smiles)]
    print(f'Data size after dropping conflicting duplicates = {len(data):,}')
    print()

    # Save data
    save_path.parent.mkdir(parents=True, exist_ok=True)
    data.to_csv(save_path, index=False)

    print(f'Final data size = {len(data):,}')
    print(f'Number of hits = {sum(data[activity_column] == 1):,}')
    print(f'Number of non-hits = {sum(data[activity_column] == 0):,}')

    # Save hits
    if save_hits_path is not None:
        save_hits_path.parent.mkdir(parents=True, exist_ok=True)
        hits = data[data[activity_column] == 1]
        hits.to_csv(save_hits_path, index=False)


if __name__ == '__main__':
    from tap import tapify

    tapify(process_data)
