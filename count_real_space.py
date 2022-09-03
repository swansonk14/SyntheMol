"""Counts reactions in REAL space."""
from collections import Counter
from multiprocessing import Pool
from pathlib import Path

import pandas as pd
from tap import Tap
from tqdm import tqdm

from constants import REAL_REACTION_COL, REAL_SPACE_SIZE
from count_real_database import save_counts_as_csv


class Args(Tap):
    data_dir: Path  # Path to directory with CXSMILES files containing the REAL database.
    save_path: Path  # Path to CSV file with reaction counts will be saved.
    parallel: bool = False  # Whether to run the script in parallel across files rather than sequentially.

    def process_args(self) -> None:
        self.save_path.parent.mkdir(parents=True, exist_ok=True)


def count_real_space_for_file(path: Path) -> Counter:
    """Counts reactions for a single REAL space file."""
    # Load REAL data file
    data = pd.read_csv(path, sep='\t', usecols=[REAL_REACTION_COL])

    # Count reactions
    reaction_counts = Counter(data[REAL_REACTION_COL])

    return reaction_counts


def count_real_space(args: Args) -> None:
    """Counts reactions in REAL space."""
    # Get paths to data files
    data_paths = sorted(args.data_dir.rglob('*.cxsmiles.bz2'))
    print(f'Number of files = {len(data_paths):,}')

    # Create combined dictionary
    combined_reaction_counts = Counter()

    # Set up map function
    if args.parallel:
        pool = Pool()
        map_fn = pool.imap
    else:
        pool = None
        map_fn = map

    # Loop through all REAL space files
    with tqdm(total=REAL_SPACE_SIZE) as progress_bar:
        for reaction_counts in map_fn(count_real_space_for_file, data_paths):
            # Merge counts
            combined_reaction_counts.update(reaction_counts)

            # Update progress bar
            progress_bar.update(reaction_counts.total())

    if pool is not None:
        pool.close()

    # Save counts
    save_counts_as_csv(
        counts=combined_reaction_counts,
        count_name='Reaction',
        save_path=args.save_path
    )


if __name__ == '__main__':
    count_real_space(Args().parse_args())
