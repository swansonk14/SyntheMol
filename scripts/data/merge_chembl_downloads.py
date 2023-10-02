"""Merges and deduplicates CSV files downloaded from ChEMBL using various search terms."""
from collections import defaultdict
from pathlib import Path

import pandas as pd
from rdkit import Chem

from synthemol.constants import CHEMBL_SMILES_COL


def merge_chembl_downloads(
        data_paths: list[Path],
        labels: list[str],
        save_path: Path,
        smiles_column: str = CHEMBL_SMILES_COL
) -> None:
    """Merges CSV files downloaded from ChEMBL using various search terms.

    :param data_paths: Paths to CSV files containing downloaded ChEMBL data.
    :param labels: Names of each file in data_paths for the purpose of labeling where the molecules came from.
    :param save_path: Path to CSV file where combined ChEMBL data will be saved.
    :param smiles_column: Name of column containing SMILES in the files in data_paths.
    """
    if len(data_paths) != len(labels):
        raise ValueError('Must have an equal number of data paths and labels.')

    # Load data
    datasets = []
    for data_path, label in zip(data_paths, labels):
        dataset = pd.read_csv(data_path, sep=';')
        dataset['label'] = label
        print(f'{data_path.name} size = {len(dataset):,}')
        datasets.append(dataset)

    # Merge datasets
    data = pd.concat(datasets)

    # Only keep compounds with SMILES
    data.dropna(subset=[smiles_column], inplace=True)

    # Compute canonical smiles
    data['smiles'] = [Chem.MolToSmiles(Chem.MolFromSmiles(smiles)) for smiles in data[smiles_column]]
    del data[smiles_column]

    # Map SMILES to labels
    smiles_to_labels = defaultdict(set)
    for smiles, label in zip(data['smiles'], data['label']):
        smiles_to_labels[smiles].add(label)

    # Deduplicate by canonical SMILES
    data.drop_duplicates(subset='smiles', inplace=True)

    # Add labels to rows
    data['labels'] = [';'.join(sorted(smiles_to_labels[smiles])) for smiles in data['smiles']]
    del data['label']

    # Sort by canonical SMILES
    data.sort_values(by='smiles', inplace=True)

    # Change order of columns
    new_columns = ['smiles', 'labels']
    data = data[new_columns + list(data.columns)[:-len(new_columns)]]

    print(f'Merged data size = {len(data):,}')

    # Save data
    save_path.parent.mkdir(parents=True, exist_ok=True)
    data.to_csv(save_path, index=False)


if __name__ == '__main__':
    from tap import tapify

    tapify(merge_chembl_downloads)
