"""Counts reactions and reagents in REAL space."""
from collections import Counter
from functools import partial
from multiprocessing import Pool
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
from tap import Tap
from tqdm import tqdm

from constants import REAL_REACTION_COL, REAL_REAGENT_COLS, REAL_SPACE_SIZE
from count_real_database import save_counts_as_csv
from real_reactions import REAL_REACTIONS


class Args(Tap):
    data_dir: Path  # Path to directory with CXSMILES files containing the REAL database.
    save_dir: Path  # Path to directory where reaction and reagent counts will be saved.
    fragment_path: Optional[Path] = None  # If provided, only count reactions and reagents that contain the fragments in this file.
    fragment_id_column: str = 'Reagent_ID'  # Column in fragment file that contains fragment IDs.
    only_selected_reactions: bool = False  # If True, only count reactions that are in the selected reactions in real_reactions.py.

    def process_args(self) -> None:
        self.save_dir.mkdir(parents=True, exist_ok=True)


USE_COLS = [REAL_REACTION_COL] + REAL_REAGENT_COLS
REAL_REACTION_IDS = {reaction.id for reaction in REAL_REACTIONS}


def count_real_space_for_file(
        path: Path,
        fragment_set: Optional[set] = None,
        only_selected_reactions: bool = False
) -> tuple[Counter, Counter, int, int]:
    """Counts reactions and reagents for a single REAL space file."""
    # Load REAL data file
    data = pd.read_csv(path, sep='\t', usecols=USE_COLS)

    # Get number of molecules in file
    num_molecules_in_file = len(data)

    # Optionally filter by fragments
    if fragment_set is not None:
        data = data[data[REAL_REAGENT_COLS].isin(fragment_set).all(axis=1)]

    # Optionally filter by selected reactions
    if only_selected_reactions:
        data = data[data[REAL_REACTION_COL].isin(REAL_REACTION_IDS)]

    # Get number of molecules that will be counted
    num_molecules_counted = len(data)

    # Count reactions
    reaction_counts = Counter(data[REAL_REACTION_COL])

    # Count reagents
    reagent_counts = Counter()
    for reagents in data[REAL_REAGENT_COLS].itertuples(index=False):
        unique_reagents = {reagent for reagent in reagents if not np.isnan(reagent)}
        reagent_counts.update(unique_reagents)

    return reaction_counts, reagent_counts, num_molecules_in_file, num_molecules_counted


def count_real_space(args: Args) -> None:
    """Counts reactions and reagents in REAL space."""
    # Get paths to data files
    data_paths = sorted(args.data_dir.rglob('*.cxsmiles.bz2'))
    print(f'Number of files = {len(data_paths):,}')

    # Optionally get set of fragments to filter by
    if args.fragment_path is not None:
        fragment_set = set(pd.read_csv(args.fragment_path)[args.fragment_id_column])
        print(f'Number of fragments = {len(fragment_set):,}')
        fragment_set |= {np.nan}
    else:
        fragment_set = None

    # Optionally only count selected reactions
    if args.only_selected_reactions:
        print(f'Number of selected reactions = {len(REAL_REACTION_IDS):,}')

    # Set up function to count reactions and reagents for a single file
    count_real_space_for_file_fn = partial(
        count_real_space_for_file,
        fragment_set=fragment_set,
        only_selected_reactions=args.only_selected_reactions
    )

    # Create combined counters
    combined_reaction_counts = Counter()
    combined_reagent_counts = Counter()
    total_num_molecules = 0
    total_num_molecules_counted = 0

    # Loop through all REAL space files
    with Pool() as pool:
        with tqdm(total=REAL_SPACE_SIZE) as progress_bar:
            for reaction_counts, reagent_counts, num_molecules_in_file, num_molecules_counted in pool.imap(
                    count_real_space_for_file_fn,
                    data_paths
            ):
                # Merge counts
                combined_reaction_counts.update(reaction_counts)
                combined_reagent_counts.update(reagent_counts)
                total_num_molecules += num_molecules_in_file
                total_num_molecules_counted += num_molecules_counted

                # Update progress bar
                progress_bar.update(num_molecules_in_file)

    print(f'Total number of molecules = {total_num_molecules:,}')
    print(f'Total number of molecules with selected fragments/reactions = {total_num_molecules_counted:,}')

    # Create reaction counts DataFrame
    combined_reaction_counts_data = pd.DataFrame(data=[
        {
            'reaction': reaction,
            'count': count,
            'percent': 100 * count / REAL_SPACE_SIZE,
        } for reaction, count in combined_reaction_counts.items()
    ])

    # Sort data by count from largest to smallest
    combined_reaction_counts_data.sort_values(by='count', ascending=False, inplace=True)

    # Add cumulative sum and cumulative percent
    combined_reaction_counts_data['cumulative_count'] = np.cumsum(combined_reaction_counts_data['count'])
    combined_reaction_counts_data['cumulative_percent'] = 100 * combined_reaction_counts_data['cumulative_count'] / REAL_SPACE_SIZE

    # Save reaction counts
    combined_reaction_counts_data.to_csv(args.save_dir / 'real_space_reaction_counts.csv', index=False)

    # Create reagent counts DataFrame
    combined_reagent_counts_data = pd.DataFrame([
        {
            'building_block': reagent,
            'count': count,
            'percent': 100 * count / REAL_SPACE_SIZE
        } for reagent, count in combined_reagent_counts.items()
    ])

    # Sort reagent counts by count from largest to smallest
    combined_reagent_counts_data.sort_values(by='count', ascending=False, inplace=True)

    # Save reagent counts
    combined_reagent_counts_data.to_csv(args.save_dir / 'real_space_reagent_counts.csv', index=False)


if __name__ == '__main__':
    count_real_space(Args().parse_args())
