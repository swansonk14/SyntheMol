"""Sample molecules uniformly at random from REAL space."""
from functools import partial
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd
from tap import Tap
from tqdm import tqdm

from constants import REAL_SMILES_COL, REAL_SPACE_SIZE


class Args(Tap):
    data_dir: Path  # Path to directory with CXSMILES files containing the REAL database.
    save_path: Path  # Path to CSV file where sampled molecules will be saved.
    num_molecules: int  # Number of molecules to sample.
    parallel: bool = False  # Whether to run the script in parallel across files rather than sequentially.

    def process_args(self) -> None:
        self.save_path.parent.mkdir(parents=True, exist_ok=True)


def sample_real_space_for_file(path: Path, sample_proportion: float) -> tuple[list[str], int]:
    """Sample molecules uniformly at random from a single REAL space file proportional to the file size."""
    # Load REAL data file
    data = pd.read_csv(path, sep='\t', usecols=[REAL_SMILES_COL])

    # Sample molecules
    rng = np.random.default_rng(seed=abs(hash(path.stem)))
    probs = rng.uniform(size=len(data))
    sampled_smiles = list(data[REAL_SMILES_COL][probs < sample_proportion])

    return sampled_smiles, len(data)


def sample_real_space(args: Args) -> None:
    """Sample molecules uniformly at random from REAL space."""
    # Get paths to data files
    data_paths = sorted(args.data_dir.rglob('*.cxsmiles.bz2'))
    print(f'Number of files = {len(data_paths):,}')

    # Determine sample proportion and double it to ensure we get enough molecules
    sample_proportion = args.num_molecules / REAL_SPACE_SIZE * 2
    print(f'Sample proportion = {sample_proportion:.3e}')

    # Create combined SMILES list
    combined_smiles = []

    # Set up map function
    if args.parallel:
        pool = Pool()
        map_fn = pool.imap
    else:
        pool = None
        map_fn = map

    sample_real_space_for_file_partial = partial(sample_real_space_for_file, sample_proportion=sample_proportion)

    # Loop through all REAL space files
    with tqdm(total=REAL_SPACE_SIZE) as progress_bar:
        for sampled_smiles, data_size in map_fn(sample_real_space_for_file_partial, data_paths):
            # Merge SMILES
            combined_smiles.extend(sampled_smiles)

            # Update progress bar
            progress_bar.update(data_size)

    if pool is not None:
        pool.close()

    print(f'Number of sampled molecules = {len(combined_smiles):,}')

    # Extract desired number of sampled molecules
    rng = np.random.default_rng(seed=0)
    sampled_smiles = rng.choice(combined_smiles, size=args.num_molecules, replace=False)
    print(f'Number of molecules after further sampling = {len(sampled_smiles):,}')

    # Save sampled molecules
    data = pd.DataFrame({REAL_SMILES_COL: combined_smiles})
    data.to_csv(args.save_path, index=False)


if __name__ == '__main__':
    sample_real_space(Args().parse_args())
