"""Counts reactions and building blocks in REAL space."""
from collections import Counter
from functools import partial
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd
from tqdm import tqdm

from synthemol.constants import (
    REAL_BUILDING_BLOCK_COLS,
    REAL_REACTION_COL,
    REAL_BUILDING_BLOCK_ID_COL,
    REAL_SPACE_SIZE
)
from synthemol.reactions import REACTIONS


USE_COLS = [REAL_REACTION_COL] + REAL_BUILDING_BLOCK_COLS
REAL_REACTION_IDS = {reaction.id for reaction in REACTIONS}


def count_real_space_for_file(
        path: Path,
        building_block_set: set | None = None,
        only_selected_reactions: bool = False
) -> tuple[Counter, Counter, int, int]:
    """Counts reactions and building blocks for a single REAL space file.

    :param path: Path to a REAL space file.
    :param building_block_set: Set of building blocks to filter by.
    :param only_selected_reactions: Whether to only count reactions in REAL_REACTION_IDS.
    :return: A tuple containing the reaction counts, building block counts,
             number of molecules in the file, and number of molecules counted.
    """
    # Load REAL data file
    data = pd.read_csv(path, sep='\t', usecols=USE_COLS)

    # Get number of molecules in file
    num_molecules_in_file = len(data)

    # Optionally filter by building blocks
    if building_block_set is not None:
        data = data[data[REAL_BUILDING_BLOCK_COLS].isin(building_block_set).all(axis=1)]

    # Optionally filter by selected reactions
    if only_selected_reactions:
        data = data[data[REAL_REACTION_COL].isin(REAL_REACTION_IDS)]

    # Get number of molecules that will be counted
    num_molecules_counted = len(data)

    # Count reactions
    reaction_counts = Counter(data[REAL_REACTION_COL])

    # Count building blocks
    building_block_counts = Counter()
    for building_blocks in data[REAL_BUILDING_BLOCK_COLS].itertuples(index=False):
        unique_building_blocks = {building_block for building_block in building_blocks if not np.isnan(building_block)}
        building_block_counts.update(unique_building_blocks)

    return reaction_counts, building_block_counts, num_molecules_in_file, num_molecules_counted


def count_real_space(
        data_dir: Path,
        save_dir: Path,
        building_blocks_path: Path | None = None,
        building_block_id_column: str = REAL_BUILDING_BLOCK_ID_COL,
        only_selected_reactions: bool = False
) -> None:
    """Counts reactions and building blocks in REAL space.

    :param data_dir: Path to directory with CXSMILES files containing the REAL database.
    :param save_dir: Path to directory where reaction and building block counts will be saved.
    :param building_blocks_path: If provided, only count reactions and building blocks that contain the building blocks in this file.
    :param building_block_id_column: Column in building block file that contains building block IDs.
    :param only_selected_reactions: If True, only count reactions that are in the selected reactions in real_reactions.py.
    """
    # Get paths to data files
    data_paths = sorted(data_dir.rglob('*.cxsmiles.bz2'))
    print(f'Number of files = {len(data_paths):,}')

    # Optionally get set of building blocks to filter by
    if building_blocks_path is not None:
        building_block_set = set(pd.read_csv(building_blocks_path)[building_block_id_column])
        print(f'Number of building blocks = {len(building_block_set):,}')
        building_block_set |= {np.nan}
    else:
        building_block_set = None

    # Optionally only count selected reactions
    if only_selected_reactions:
        print(f'Number of selected reactions = {len(REAL_REACTION_IDS):,}')

    # Set up function to count reactions and building blocks for a single file
    count_real_space_for_file_fn = partial(
        count_real_space_for_file,
        building_block_set=building_block_set,
        only_selected_reactions=only_selected_reactions
    )

    # Create combined counters
    combined_reaction_counts = Counter()
    combined_building_block_counts = Counter()
    total_num_molecules = 0
    total_num_molecules_counted = 0

    # Loop through all REAL space files
    with Pool() as pool:
        with tqdm(total=REAL_SPACE_SIZE) as progress_bar:
            for reaction_counts, building_block_counts, num_molecules_in_file, num_molecules_counted in pool.imap(
                    count_real_space_for_file_fn,
                    data_paths
            ):
                # Merge counts
                combined_reaction_counts.update(reaction_counts)
                combined_building_block_counts.update(building_block_counts)
                total_num_molecules += num_molecules_in_file
                total_num_molecules_counted += num_molecules_counted

                # Update progress bar
                progress_bar.update(num_molecules_in_file)

    print(f'Total number of molecules = {total_num_molecules:,}')
    print(f'Total number of molecules with selected building blockss/reactions = {total_num_molecules_counted:,}')

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
    combined_reaction_counts_data.to_csv(
        save_dir / f'reaction_counts{"_selected" if only_selected_reactions else ""}.csv',
        index=False
    )

    # Create building block counts DataFrame
    combined_building_block_counts_data = pd.DataFrame([
        {
            'building_block': building_block,
            'count': count,
            'percent': 100 * count / REAL_SPACE_SIZE
        } for building_block, count in combined_building_block_counts.items()
    ])

    # Sort building block counts by count from largest to smallest
    combined_building_block_counts_data.sort_values(by='count', ascending=False, inplace=True)

    # Save building block counts
    combined_building_block_counts_data.to_csv(
        save_dir / f'building_block_counts{"_selected" if only_selected_reactions else ""}.csv',
        index=False
    )


if __name__ == '__main__':
    from tap import tapify

    tapify(count_real_space)
