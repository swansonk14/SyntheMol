"""Sample molecules uniformly at random from REAL space."""
from functools import partial
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd
from tqdm import tqdm

from synthemol.constants import (
    REAL_BUILDING_BLOCK_COLS,
    REAL_ID_COL,
    REAL_REACTION_COL,
    REAL_SMILES_COL,
    REAL_SPACE_SIZE,
)


USE_COLS_MULTIPLE_ID = (
    [REAL_REACTION_COL] + REAL_BUILDING_BLOCK_COLS + [REAL_SMILES_COL]
)
USE_COLS_SINGLE_ID = [REAL_ID_COL, REAL_SMILES_COL]


def sample_real_space_for_file(
    path: Path, sample_proportion: float, use_cols: list[str]
) -> tuple[pd.DataFrame, int]:
    """Sample molecules uniformly at random from a single REAL space file proportional to the file size.

    :param path: Path to a REAL space file.
    :param sample_proportion: Proportion of molecules to sample from the file.
    :param use_cols: Columns to use from the REAL space file.
    :return: A tuple containing the sampled molecules and the number of molecules in the file.
    """
    # Load REAL data file
    data = pd.read_csv(path, sep="\t", usecols=use_cols)

    # Sample rows
    rng = np.random.default_rng(seed=abs(hash(path.stem)))
    probs = rng.uniform(size=len(data))
    sampled_data = data.iloc[probs < sample_proportion]

    return sampled_data, len(data)


def sample_real_space(
    data_dir: Path, save_path: Path, num_molecules: int, single_id_column: bool = False
) -> None:
    """Sample molecules uniformly at random from REAL space.

    :param data_dir: Path to directory with CXSMILES files containing the REAL database.
    :param save_path: Path to CSV file where sampled molecules will be saved.
    :param num_molecules: Number of molecules to sample.
    :param single_id_column: Whether the reaction and building blocks are in a single ID column (newer versions of REAL)
                             or in separate columns (older versions of REAL).
    """
    # Get paths to data files
    data_paths = sorted(data_dir.rglob("*.cxsmiles.bz2"))
    print(f"Number of files = {len(data_paths):,}")

    # Determine sample proportion and double it to ensure we get enough molecules
    sample_proportion = num_molecules / REAL_SPACE_SIZE * 2
    print(f"Sample proportion = {sample_proportion:.3e}")

    sample_real_space_for_file_partial = partial(
        sample_real_space_for_file,
        sample_proportion=sample_proportion,
        use_cols=USE_COLS_SINGLE_ID if single_id_column else USE_COLS_MULTIPLE_ID,
    )

    # Loop through all REAL space files
    data = []
    with Pool() as pool:
        with tqdm(total=REAL_SPACE_SIZE) as progress_bar:
            for sampled_data, data_size in pool.imap(
                sample_real_space_for_file_partial, data_paths
            ):
                # Gather sampled data
                data.append(sampled_data)

                # Update progress bar
                progress_bar.update(data_size)

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

    tapify(sample_real_space)
