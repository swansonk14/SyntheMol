"""Determines which REAL building blocks can be used in which REAL reactions."""
import pickle
from collections import defaultdict
from functools import partial
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd
from tqdm import tqdm

from count_real_space import expand_single_id_column
from synthemol.constants import (
    REAL_BUILDING_BLOCK_COLS,
    REAL_ID_COL,
    REAL_REACTION_COL,
    REAL_SPACE_SIZE
)


USE_COLS_MULTIPLE_ID = [REAL_REACTION_COL] + REAL_BUILDING_BLOCK_COLS
USE_COLS_SINGLE_ID = [REAL_ID_COL]


def map_reactions_for_file(
        path: Path,
        single_id_column: bool = False
) -> tuple[str, int, dict[int, dict[int, set[int]]]]:
    """Computes a mapping from reactions to building blocks for a single REAL file.

    :param path: Path to a REAL file.
    :param single_id_column: Whether the reaction and building blocks are in a single ID column (newer versions of REAL)
                             or in separate columns (older versions of REAL).
    :return: A tuple containing the name of the file, the number of molecules in the file,
             and a mapping from reaction ID to reactant index to reaction type to valid building block IDs.
    """
    # Create mapping from reaction ID to reactant index to valid building block IDs
    reaction_to_reactants_to_building_blocks: dict[int, dict[int, set[int]]] = defaultdict(lambda: defaultdict(set))

    # Load REAL data file (ensures cols are in the order of usecols for itertuples below)
    use_cols = USE_COLS_SINGLE_ID if single_id_column else USE_COLS_MULTIPLE_ID
    data = pd.read_csv(path, sep='\t', usecols=use_cols)[use_cols]

    # If single ID column, expand to multiple ID columns and ensure only those columns remain in order
    if single_id_column:
        data = expand_single_id_column(data)[USE_COLS_MULTIPLE_ID]

    # Update mapping
    for row_tuple in data.itertuples(index=False):
        reaction, building_blocks = row_tuple[0], row_tuple[1:]

        for reactant, building_block in enumerate(building_blocks):
            if not np.isnan(building_block):
                reaction_to_reactants_to_building_blocks[reaction][reactant].add(int(building_block))

    # Convert to regular dict for compatibility with multiprocessing
    reaction_to_reactants_to_building_blocks = {
        reaction: {
            reactant: building_block_ids
            for reactant, building_block_ids in reactant_to_building_blocks.items()
        }
        for reaction, reactant_to_building_blocks in reaction_to_reactants_to_building_blocks.items()
    }

    # Get data name
    name = path.stem.split('.')[0]

    # Get number of molecules
    num_molecules = len(data)

    return name, num_molecules, reaction_to_reactants_to_building_blocks


def map_real_reactions_to_building_blocks(
        data_dir: Path,
        save_path: Path,
        single_id_column: bool = False
) -> None:
    """Determines which REAL building blocks can be used in which REAL reactions.

    :param data_dir: Path to directory with CXSMILES files containing the REAL database.
    :param save_path: Path to PKL file where mapping will be saved.
    :param single_id_column: Whether the reaction and building blocks are in a single ID column (newer versions of REAL)
                             or in separate columns (older versions of REAL).
    """
    # Get paths to data files
    data_paths = sorted(data_dir.rglob('*.cxsmiles.bz2'))
    print(f'Number of files = {len(data_paths):,}')

    # Create combined dictionary
    combined_reaction_to_reactants_to_building_blocks = defaultdict(lambda: defaultdict(set))

    # Set up mapping function
    map_reactions_for_file_fn = partial(
        map_reactions_for_file,
        single_id_column=single_id_column
    )

    # Loop through all REAL space files
    num_files = total_num_molecules = 0
    with Pool() as pool:
        with tqdm(total=REAL_SPACE_SIZE) as progress_bar:
            for name, num_molecules, reaction_to_reactants_to_building_blocks in pool.imap(map_reactions_for_file_fn, data_paths):
                num_files += 1
                total_num_molecules += num_molecules
                print(f'{name}: file num = {num_files:,} / {len(data_paths):,} | '
                      f'num mols = {num_molecules:,} | cumulative num mols = {total_num_molecules:,}\n')

                # Merge dictionary with combined dictionary
                for reaction, reactant_to_building_blocks in reaction_to_reactants_to_building_blocks.items():
                    for reactant, building_blocks in reactant_to_building_blocks.items():
                        combined_reaction_to_reactants_to_building_blocks[reaction][reactant] |= building_blocks

                # Update progress bar
                progress_bar.update(num_molecules)

    print(f'Total number of molecules = {total_num_molecules:,}')

    # Convert sets to sorted lists
    combined_reaction_to_reactants_to_building_blocks = {
        reaction: {
            reactant: sorted(building_blocks)
            for reactant, building_blocks in reactant_to_building_blocks.items()
        }
        for reaction, reactant_to_building_blocks in combined_reaction_to_reactants_to_building_blocks.items()
    }

    # Save mapping
    save_path.parent.mkdir(parents=True, exist_ok=True)

    with open(save_path, 'wb') as f:
        pickle.dump(combined_reaction_to_reactants_to_building_blocks, f)


if __name__ == '__main__':
    from tap import tapify

    tapify(map_real_reactions_to_building_blocks)
