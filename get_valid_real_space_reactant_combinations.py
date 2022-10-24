"""Gets the valid reactant combinations for reactions across the whole REAL space."""
import pickle
from functools import partial
from multiprocessing import Pool
from pathlib import Path

import pandas as pd
from tap import Tap
from tqdm import tqdm

from constants import REAL_REACTION_COL, REAL_REAGENT_COLS, REAL_SPACE_SIZE
from real_reactions import REAL_REACTIONS


REACTIONS = [reaction.id for reaction in REAL_REACTIONS]


class Args(Tap):
    data_dir: Path  # Path to directory with CXSMILES files containing the REAL database.
    save_dir: Path  # Path to directory where the valid reactant combinations will be saved.
    reactions: list[int] = REACTIONS  # Reaction IDs to get valid reactant combinations for.
    parallel: bool = False  # Whether to run the script in parallel across files rather than sequentially.

    def process_args(self) -> None:
        self.save_dir.mkdir(parents=True, exist_ok=True)


def get_valid_real_space_reactant_combinations_for_file(path: Path, reactions: list[int]) -> tuple[dict[int, set[tuple[int, ...]]], int]:
    """Gets the valid reactant combinations for reactions for a single REAL space file."""
    # Load REAL data file
    data = pd.read_csv(path, sep='\t', usecols=[REAL_REACTION_COL, *REAL_REAGENT_COLS])
    data_size = len(data)

    # Create dictionary of valid reactants
    valid_reactants = {}

    # Loop through reactions
    for reaction in reactions:
        # Limit data to reaction
        reaction_data = data[data[REAL_REACTION_COL] == reaction]

        # Determine non-na reactant columns
        reactant_cols = [col for col in REAL_REAGENT_COLS if reaction_data[col].notna().any()]

        # Get valid reactant combinations
        valid_reactants[reaction] = set(reaction_data[reactant_cols].itertuples(index=False, name=None))

    return valid_reactants, data_size


def get_valid_real_space_reactant_combinations(args: Args) -> None:
    """Gets the valid reactant combinations for reactions across the whole REAL space."""
    # Get paths to data files
    data_paths = sorted(args.data_dir.rglob('*.cxsmiles.bz2'))
    print(f'Number of files = {len(data_paths):,}')

    # Create combined valid reactants dictionary
    all_valid_reactants = {}

    # Set up map function
    if args.parallel:
        pool = Pool()
        map_fn = pool.imap
    else:
        pool = None
        map_fn = map

    # Set up reactants for file function
    reactants_for_file_fn = partial(get_valid_real_space_reactant_combinations_for_file, reactions=args.reactions)

    # Loop through all REAL space files
    with tqdm(total=REAL_SPACE_SIZE) as progress_bar:
        for valid_reactants, data_size in map_fn(reactants_for_file_fn, data_paths):
            # Merge valid reactants
            for reaction, reactants in valid_reactants.items():
                if reaction not in all_valid_reactants:
                    all_valid_reactants[reaction] = reactants
                else:
                    all_valid_reactants[reaction].update(reactants)

            # Update progress bar
            progress_bar.update(data_size)

    if pool is not None:
        pool.close()

    # Print number of valid reactants for each reaction
    for reaction, reactants in all_valid_reactants.items():
        print(f'Number of valid reactants for reaction {reaction} = {len(reactants):,}')

    # Save valid reactants with one file per reaction
    for reaction, reactants in all_valid_reactants.items():
        with open(args.save_dir / f'reaction_{reaction}.pkl', 'wb') as f:
            pickle.dump(reactants, f)


if __name__ == '__main__':
    get_valid_real_space_reactant_combinations(Args().parse_args())
