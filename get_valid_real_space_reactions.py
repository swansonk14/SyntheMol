"""Gets the valid reactant and product combinations for reactions across the whole REAL space."""
import pickle
from functools import partial
from multiprocessing import Pool
from pathlib import Path

import pandas as pd
from tap import Tap
from tqdm import tqdm

from constants import (
    REAL_REACTION_COL,
    REAL_REAGENT_COLS,
    REAL_SMILES_COL,
    REAL_SPACE_SIZE
)
from real_reactions import REAL_REACTIONS


REACTIONS = [reaction.id for reaction in REAL_REACTIONS]


class Args(Tap):
    data_dir: Path  # Path to directory with CXSMILES files containing the REAL database.
    save_dir: Path  # Path to directory where the valid reactant combinations will be saved.
    reactions: list[int] = REACTIONS  # Reaction IDs to get valid reactant combinations for.
    parallel: bool = False  # Whether to run the script in parallel across files rather than sequentially.

    def process_args(self) -> None:
        self.save_dir.mkdir(parents=True, exist_ok=True)


def get_valid_real_space_reactant_combinations_for_file(
        path: Path, reactions: list[int]
) -> tuple[dict[int, set[tuple[int, ..., str]]], int]:
    """Gets the valid reactant and product combinations for reactions for a single REAL space file."""
    # Load REAL data file
    data = pd.read_csv(path, sep='\t', usecols=[REAL_SMILES_COL, REAL_REACTION_COL, *REAL_REAGENT_COLS])
    data_size = len(data)

    # Create dictionary of valid reactants
    valid_reactant_product_combos = {}

    # Loop through reactions
    for reaction in reactions:
        # Limit data to reaction
        reaction_data = data[data[REAL_REACTION_COL] == reaction]

        # Determine non-na reactant columns
        reactant_cols = [col for col in REAL_REAGENT_COLS if reaction_data[col].notna().any()]

        # Set up all cols to extract
        cols = [*reactant_cols, REAL_SMILES_COL]

        # Get valid reactant combinations
        valid_reactant_product_combos[reaction] = set(reaction_data[cols].itertuples(index=False, name=None))

    return valid_reactant_product_combos, data_size


def get_valid_real_space_reactant_combinations(args: Args) -> None:
    """Gets the valid reactant and product combinations for reactions across the whole REAL space."""
    # Get paths to data files
    data_paths = sorted(args.data_dir.rglob('*.cxsmiles.bz2'))
    print(f'Number of files = {len(data_paths):,}')

    # Create combined valid reactant/product dictionary
    all_valid_reactant_product_combos = {}

    # Set up map function
    if args.parallel:
        pool = Pool()
        map_fn = pool.imap
    else:
        pool = None
        map_fn = map

    # Set up reactant/product for file function
    reactants_for_file_fn = partial(get_valid_real_space_reactant_combinations_for_file, reactions=args.reactions)

    # Loop through all REAL space files
    with tqdm(total=REAL_SPACE_SIZE) as progress_bar:
        for valid_reactant_product_combos, data_size in map_fn(reactants_for_file_fn, data_paths):
            # Merge valid reactant/product combos
            for reaction, reactant_product_combo in valid_reactant_product_combos.items():
                if reaction not in all_valid_reactant_product_combos:
                    all_valid_reactant_product_combos[reaction] = reactant_product_combo
                else:
                    all_valid_reactant_product_combos[reaction].update(reactant_product_combo)

            # Update progress bar
            progress_bar.update(data_size)

    if pool is not None:
        pool.close()

    # Print number of valid reactants for each reaction
    for reaction, reactant_product_combo in all_valid_reactant_product_combos.items():
        print(f'Number of valid reactant/product combinations for reaction {reaction} = {len(reactant_product_combo):,}')

    # Save valid reactants with one file per reaction
    for reaction, reactant_product_combo in all_valid_reactant_product_combos.items():
        with open(args.save_dir / f'reaction_{reaction}.pkl', 'wb') as f:
            pickle.dump(reactant_product_combo, f)


if __name__ == '__main__':
    get_valid_real_space_reactant_combinations(Args().parse_args())
