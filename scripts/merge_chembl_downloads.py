""" Merges and deduplicates CSV files downloaded from ChEMBL using various search terms."""
from collections import defaultdict
from pathlib import Path

import pandas as pd
from rdkit import Chem
from tap import Tap


class Args(Tap):
    data_paths: list[Path]  # Paths to CSV files containing downloaded ChEMBL data.
    labels: list[str]  # Names of each file in data_paths for the purpose of labeling where the molecules came from.
    save_path: Path  # Path to CSV file where combined ChEMBL data will be saved.
    id_column: str = 'ChEMBL ID'  # Name of the column containing the compound ID in the files in data_paths.
    smiles_column: str = 'Smiles'  # Name of column containing SMILES in the files in data_paths.


def merge_chembl_downloads(args: Args) -> None:
    """Merges CSV files downloaded from ChEMBL using various search terms."""
    if len(args.data_paths) != len(args.labels):
        raise ValueError('Must have an equal number of data paths and labels.')

    # Load data
    datasets = []
    for data_path, label in zip(args.data_paths, args.labels):
        dataset = pd.read_csv(data_path, sep=';')
        dataset['label'] = label
        print(f'{data_path.name} size = {len(dataset):,}')
        datasets.append(dataset)

    # Merge datasets
    data = pd.concat(datasets)

    # Only keep compounds with SMILES
    data.dropna(subset=[args.smiles_column], inplace=True)

    # Compute canonical smiles
    data['smiles'] = [Chem.MolToSmiles(Chem.MolFromSmiles(smiles)) for smiles in data[args.smiles_column]]
    del data[args.smiles_column]

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
    new_columns = ['smiles', 'chembl_smiles', 'labels']
    data = data[new_columns + list(data.columns)[:-len(new_columns)]]

    print(f'Merged data size = {len(data):,}')

    # Save data
    args.save_path.parent.mkdir(parents=True, exist_ok=True)
    data.to_csv(args.save_path, index=False)


if __name__ == '__main__':
    merge_chembl_downloads(Args().parse_args())
