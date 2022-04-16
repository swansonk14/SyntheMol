"""Samples molecules from the REAL database."""
from functools import partial
from pathlib import Path

from multiprocessing import Pool
import numpy as np
from tap import Tap
from tqdm import tqdm


class Args(Tap):
    data_dir: Path  # Path to directory with CXSMILES files containing the REAL database.
    num_molecules: int  # Number of molecules to sample per CXSMILES file.
    save_dir: Path  # Path to directory where sampled molecules will be saved.

    def process_args(self) -> None:
        self.save_dir.mkdir(parents=True, exist_ok=True)


def sample_real_database_file(path: Path, num_molecules: int, save_dir: Path) -> None:
    """Sample molecules from a single REAL database file."""
    # Count number of rows in file
    with open(path) as f:
        num_rows = sum(1 for _ in tqdm(f, desc='Counting rows')) - 1  # minus 1 to account for header

    # Sample row indices
    sampled_row_indices = set(np.random.choice(num_rows, size=num_molecules, replace=False))

    # Select and save sampled data
    with open(path) as full_file, open(save_dir / f'{path.stem}_sampled.cxsmiles', 'w') as sample_file:
        # Copy header
        header = next(full_file)
        sample_file.write(header)

        # Copy selected rows
        for index, row in tqdm(enumerate(f), total=num_rows, desc='Saving rows'):
            if index in sampled_row_indices:
                sample_file.write(row)


def sample_real_database(args: Args) -> None:
    """Samples molecules from the REAL database."""
    # Get paths to data files
    data_paths = list(args.data_dir.glob('*.cxsmiles'))

    # Set up sampling function
    sample_real_database_file_fn = partial(
        sample_real_database_file,
        num_molecules=args.num_molecules,
        save_dir=args.save_dir
    )

    # Loop through all REAL database files and sample molecules
    with Pool() as pool:
        for _ in tqdm(pool.imap_unordered(sample_real_database_file_fn, data_paths), desc='Loop over files'):
            pass


if __name__ == '__main__':
    sample_real_database(Args().parse_args())
