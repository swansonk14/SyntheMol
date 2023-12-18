"""Sample molecules uniformly at random from WuXi GalaXi."""
from pathlib import Path

import numpy as np
import pandas as pd
from tqdm import tqdm

from synthemol.constants import WUXI_GALAXI_SIZE, WUXI_ID_COL, WUXI_SMILES_COL


def sample_wuxi_galaxi(data_dir: Path, save_path: Path, num_molecules: int,) -> None:
    """Sample molecules uniformly at random from WuXi GalaXi.

    :param data_dir: Path to directory with CSV files containing the WuXi GalaXi.
    :param save_path: Path to CSV file where sampled molecules will be saved.
    :param num_molecules: Number of molecules to sample.
    """
    # Get paths to data files
    data_paths = sorted(data_dir.rglob("*.csv"))
    print(f"Number of files = {len(data_paths):,}")

    # Determine sample proportion and slightly increase it to ensure we get enough molecules
    sample_proportion = num_molecules / WUXI_GALAXI_SIZE * 1.1
    print(f"Sample proportion = {sample_proportion:.3e}")

    # Loop through all WuXi GalaXi files
    data = []
    count = 0
    for data_path in tqdm(data_paths):
        # Load WuXi GalaXi file
        file_data = pd.read_csv(data_path, header=None, names=[WUXI_ID_COL, WUXI_SMILES_COL])
        count += len(file_data)

        # Sample rows
        rng = np.random.default_rng(seed=abs(hash(data_path.stem)))
        probs = rng.uniform(size=len(file_data))
        sampled_data = file_data.iloc[probs < sample_proportion]

        data.append(sampled_data)

    # Print WuXi GalaXi size
    print(f"WuXi GalaXi size = {count:,}")

    # Merge sampled data
    data = pd.concat(data, ignore_index=True)
    print(f"Number of sampled molecules = {len(data):,}")

    # Extract desired number of sampled molecules
    data = data.sample(n=num_molecules, random_state=0, replace=False)
    print(f"Number of molecules after further sampling = {len(data):,}")

    # Save sampled molecules
    save_path.parent.mkdir(parents=True, exist_ok=True)
    data.to_csv(save_path, index=False)


if __name__ == "__main__":
    from tap import tapify

    tapify(sample_wuxi_galaxi)
