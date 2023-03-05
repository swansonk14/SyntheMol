"""Sample molecules uniformly at random from REAL space."""
from functools import partial
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd
from tap import Tap
from tqdm import tqdm

from constants import REAL_REACTION_COL, REAL_REAGENT_COLS, REAL_SMILES_COL, REAL_SPACE_SIZE


class Args(Tap):
    data_dir: Path  # Path to directory with CXSMILES files containing the REAL database.
    save_path: Path  # Path to CSV file where sampled molecules will be saved.
    num_molecules: int  # Number of molecules to sample.
    parallel: bool = False  # Whether to run the script in parallel across files rather than sequentially.

    def process_args(self) -> None:
        self.save_path.parent.mkdir(parents=True, exist_ok=True)


USE_COLS = [REAL_REACTION_COL] + REAL_REAGENT_COLS + [REAL_SMILES_COL]


def sample_real_space_for_file(path: Path, sample_proportion: float) -> tuple[pd.DataFrame, int]:
    """Sample molecules uniformly at random from a single REAL space file proportional to the file size."""
    # Load REAL data file
    data = pd.read_csv(path, sep='\t', usecols=USE_COLS)

    # Sample rows
    rng = np.random.default_rng(seed=abs(hash(path.stem)))
    probs = rng.uniform(size=len(data))
    sampled_data = data.iloc[probs < sample_proportion]

    return sampled_data, len(data)


def sample_real_space(args: Args) -> None:
    """Sample molecules uniformly at random from REAL space."""
    # Get paths to data files
    data_paths = sorted(args.data_dir.rglob('*.cxsmiles.bz2'))
    print(f'Number of files = {len(data_paths):,}')

    # Determine sample proportion and double it to ensure we get enough molecules
    sample_proportion = args.num_molecules / REAL_SPACE_SIZE * 2
    print(f'Sample proportion = {sample_proportion:.3e}')

    # Set up map function
    if args.parallel:
        pool = Pool()
        map_fn = pool.imap
    else:
        pool = None
        map_fn = map

    sample_real_space_for_file_partial = partial(sample_real_space_for_file, sample_proportion=sample_proportion)

    # Loop through all REAL space files
    data = []
    with tqdm(total=REAL_SPACE_SIZE) as progress_bar:
        for sampled_data, data_size in map_fn(sample_real_space_for_file_partial, data_paths):
            # Gather sampled data
            data.append(sampled_data)

            # Update progress bar
            progress_bar.update(data_size)

    if pool is not None:
        pool.close()

    # Merge sampled data
    data = pd.concat(data, ignore_index=True)
    print(f'Number of sampled molecules = {len(data):,}')

    # Extract desired number of sampled molecules
    data = data.sample(n=args.num_molecules, random_state=0, replace=False)
    print(f'Number of molecules after further sampling = {len(data):,}')

    # Save sampled molecules
    data.to_csv(args.save_path, index=False)


if __name__ == '__main__':
    sample_real_space(Args().parse_args())
