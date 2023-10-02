"""Merges CSV files containing WuXi GalaXi building blocks into a single file."""
from pathlib import Path

import pandas as pd
from rdkit import Chem
from tqdm import tqdm

from synthemol.constants import (
    WUXI_BUILDING_BLOCK_ID_COL,
    WUXI_BUILDING_BLOCK_SUBSET_COL,
    WUXI_BUILDING_BLOCK_SMILES_COL
)


def merge_wuxi_building_blocks(
        data_dir: Path,
        save_path: Path,
        smiles_column: str = WUXI_BUILDING_BLOCK_SMILES_COL,
        id_column: str = WUXI_BUILDING_BLOCK_ID_COL,
        subset_column: str = WUXI_BUILDING_BLOCK_SUBSET_COL
) -> None:
    """Merges CSV files containing WuXi GalaXi building blocks into a single file.

    :param data_dir: Path to directory containing CSV files.
    :param save_path: Path to CSV file where merged data will be saved.
    :param smiles_column: Name of column containing SMILES in the files in data_dir.
    :param id_column: Name of column containing IDs in the files in data_dir.
    :param subset_column: Name of column that will contain the name of the building block subset (from file name).
    """
    # Get data paths
    data_paths = sorted(data_dir.glob('*.csv'))
    print(f'Number of building blocks files = {len(data_paths):,}')

    # Load each data file
    datasets = []
    for data_path in tqdm(data_paths):
        data_name = data_path.stem.split('_')[0]
        dataset = pd.read_csv(data_path)
        dataset[subset_column] = data_name
        datasets.append(dataset)

    # Combine data
    data = pd.concat(datasets)

    print(f'Data size = {len(data):,}')

    # Remove rows with missing IDs or SMILES
    data.dropna(subset=[smiles_column], how='any', inplace=True)
    print(f'Data size after removing rows with missing SMILES = {len(data):,}')

    # Print dataset statistics
    print(f'Number of unique building block IDs = {data[id_column].nunique():,}')
    print(f'Number of unique building block SMILES = {data[smiles_column].nunique():,}')
    print(f'Number of unique building block subsets = {data[subset_column].nunique():,}')

    # Save data
    save_path.parent.mkdir(parents=True, exist_ok=True)
    data.to_csv(save_path, index=False)


if __name__ == '__main__':
    from tap import tapify

    tapify(merge_wuxi_building_blocks)
