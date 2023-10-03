"""Merges CSV files containing WuXi GalaXi building blocks into a single file."""
from pathlib import Path

import pandas as pd
from rdkit import Chem
from tqdm import tqdm

from synthemol.constants import (
    ID_COL,
    WUXI_BUILDING_BLOCK_ID_COL,
    WUXI_BUILDING_BLOCK_SUBSET_COL,
    WUXI_BUILDING_BLOCK_SMILES_COL,
    WUXI_CORES_CORE_COL,
    WUXI_CORES_ID_COL,
    WUXI_CORES_SMILES_COL
)


def merge_wuxi_building_blocks(
        building_blocks_dir: Path,
        phase_2_cores_path: Path,
        phase_3_cores_path: Path,
        save_path: Path,
        smiles_column: str = WUXI_BUILDING_BLOCK_SMILES_COL,
        wuxi_id_column: str = WUXI_BUILDING_BLOCK_ID_COL,
        subset_column: str = WUXI_BUILDING_BLOCK_SUBSET_COL,
        id_column: str = ID_COL
) -> None:
    """Merges CSV files containing WuXi GalaXi building blocks into a single file.

    :param building_blocks_dir: Path to directory containing CSV files.
    :param phase_2_cores_path: Path to XLSX file containing phase 2 cores.
    :param phase_3_cores_path: Path to XLSX file containing phase 3 cores.
    :param save_path: Path to CSV file where merged data will be saved.
    :param smiles_column: Name of column containing SMILES in the files in data_dir.
    :param wuxi_id_column: Name of column containing WuXi IDs in the files in data_dir.
    :param subset_column: Name of column that will contain the name of the building block subset (from file name).
    :param id_column: Name of the column that will contain the building block IDs (WuXi ID + subset).
    """
    # Get building blocks paths
    data_paths = sorted(building_blocks_dir.glob('*.csv'))
    print(f'Number of building blocks paths = {len(data_paths):,}')

    # Load each building blocks file
    datasets = []
    for data_path in tqdm(data_paths):
        data_name = data_path.stem.split('_')[0]
        dataset = pd.read_csv(data_path)
        dataset[subset_column] = data_name
        datasets.append(dataset)

    # Load each cores file
    phase_2_cores = pd.read_excel(phase_2_cores_path)
    phase_3_cores = pd.read_excel(phase_3_cores_path)

    # Process each cores file
    phase_2 = pd.DataFrame({
        smiles_column: phase_2_cores[WUXI_CORES_SMILES_COL],
        wuxi_id_column: phase_2_cores[WUXI_CORES_ID_COL],
        subset_column: [f'phase_2_{core}' for core in phase_2_cores[WUXI_CORES_CORE_COL]]
    })
    phase_3 = pd.DataFrame({
        smiles_column: phase_3_cores[WUXI_CORES_SMILES_COL],
        wuxi_id_column: phase_3_cores[WUXI_CORES_ID_COL],
        subset_column: [f'phase_3_{core}' for core in phase_3_cores[WUXI_CORES_CORE_COL]]
    })

    # Combine data
    data = pd.concat(datasets + [phase_2, phase_3])

    print(f'Data size = {len(data):,}')

    # Remove rows with missing IDs or SMILES
    data.dropna(subset=[wuxi_id_column], how='any', inplace=True)
    print(f'Data size after removing rows with missing IDs = {len(data):,}')

    data.dropna(subset=[smiles_column], how='any', inplace=True)
    print(f'Data size after removing rows with missing SMILES = {len(data):,}')

    # Remove rows with invalid SMILES
    data = data[data[smiles_column].apply(Chem.MolFromSmiles).notna()]
    print(f'Data size after removing rows with invalid SMILES = {len(data):,}')

    # Create ID column with both wxid and subset
    data[id_column] = data[[wuxi_id_column, subset_column]].apply(
        lambda row: f'{row[wuxi_id_column]}_{row[subset_column]}',
        axis=1
    )

    # Print dataset statistics
    print(f'Number of unique building block IDs = {data[id_column].nunique():,}')
    print(f'Number of unique building block WuXi IDs = {data[wuxi_id_column].nunique():,}')
    print(f'Number of unique building block subsets = {data[subset_column].nunique():,}')
    print(f'Number of unique building block SMILES = {data[smiles_column].nunique():,}')

    # Save data
    save_path.parent.mkdir(parents=True, exist_ok=True)
    data.to_csv(save_path, index=False)


if __name__ == '__main__':
    from tap import tapify

    tapify(merge_wuxi_building_blocks)
