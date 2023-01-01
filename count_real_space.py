"""Counts reactions and reagents in REAL space."""
from collections import Counter
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd
from tap import Tap
from tqdm import tqdm

from constants import REAL_REACTION_COL, REAL_REAGENT_COLS, REAL_SPACE_SIZE
from count_real_database import save_counts_as_csv


class Args(Tap):
    data_dir: Path  # Path to directory with CXSMILES files containing the REAL database.
    save_dir: Path  # Path to directory where reaction and reagent counts will be saved.
    parallel: bool = False  # Whether to run the script in parallel across files rather than sequentially.

    def process_args(self) -> None:
        self.save_dir.mkdir(parents=True, exist_ok=True)


USE_COLS = [REAL_REACTION_COL] + REAL_REAGENT_COLS


def count_real_space_for_file(path: Path) -> tuple[Counter, Counter]:
    """Counts reactions and reagents for a single REAL space file."""
    # Load REAL data file
    data = pd.read_csv(path, sep='\t', usecols=USE_COLS)

    # Count reactions
    reaction_counts = Counter(data[REAL_REACTION_COL])

    # Count reagents
    reagent_counts = Counter()
    for reaction, reagent_1, reagent_2, reagent_3, reagent_4 in data.itertuples(index=False):
        reagents = [reagent_1, reagent_2, reagent_3, reagent_4]
        unique_reagents = {reagent for reagent in reagents if not np.isnan(reagent)}
        reagent_counts.update(unique_reagents)

    return reaction_counts, reagent_counts


def count_real_space(args: Args) -> None:
    """Counts reactions and reagents in REAL space."""
    # Get paths to data files
    data_paths = sorted(args.data_dir.rglob('*.cxsmiles.bz2'))
    print(f'Number of files = {len(data_paths):,}')

    # Create combined counters
    combined_reaction_counts = Counter()
    combined_reagent_counts = Counter()

    # Set up map function
    if args.parallel:
        pool = Pool()
        map_fn = pool.imap
    else:
        pool = None
        map_fn = map

    # Loop through all REAL space files
    with tqdm(total=REAL_SPACE_SIZE) as progress_bar:
        for reaction_counts, reagent_counts in map_fn(count_real_space_for_file, data_paths):
            # Merge counts
            combined_reaction_counts.update(reaction_counts)
            combined_reagent_counts.update(reagent_counts)

            # Update progress bar
            progress_bar.update(reaction_counts.total())

    if pool is not None:
        pool.close()

    # Save reaction counts
    save_counts_as_csv(
        counts=combined_reaction_counts,
        count_name='reaction',
        save_path=args.save_dir / 'real_space_reaction_counts.csv'
    )

    # Create reagent counts DataFrame
    combined_regent_counts_data = pd.DataFrame([
        {
            'reagent': reagent,
            'count': count,
            'percent_of_real_molecules': count / REAL_SPACE_SIZE
        } for reagent, count in combined_reagent_counts.items()
    ])

    # Sort reagent counts by count from largest to smallest
    combined_regent_counts_data.sort_values(by='count', ascending=False, inplace=True)

    # Save reagent counts
    combined_regent_counts_data.to_csv(args.save_dir / 'real_space_reagent_counts.csv', index=False)


if __name__ == '__main__':
    count_real_space(Args().parse_args())
