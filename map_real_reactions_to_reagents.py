"""Determines which REAL reagents can be used in which REAL reactions."""
import json
from multiprocessing import Pool
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from tap import Tap
from tqdm import tqdm


class Args(Tap):
    data_dir: Path  # Path to directory with CXSMILES files containing the REAL database.
    save_path: Path  # Path to JSON file where mapping will be saved.

    def process_args(self) -> None:
        self.save_path.parent.mkdir(parents=True, exist_ok=True)


REACTION_COL = 'reaction'
REAGENT_COLS = ['reagent1', 'reagent2', 'reagent3', 'reagent4']
TYPE_COL = 'Type'
USECOLS = [REACTION_COL] + REAGENT_COLS + [TYPE_COL]


def map_reactions_for_file(path: Path) -> tuple[int, dict[int, dict[int, dict[str, set[int]]]]]:
    """Computes the mapping for a single file."""
    # Create mapping from reaction ID to reagent number to reaction type to valid reagent IDs
    reaction_to_reagents: dict[int, dict[int, dict[str, set[int]]]] = \
        defaultdict(lambda: defaultdict(lambda: defaultdict(set)))

    # Load REAL data file (ensures cols are in the order of USECOLS for itertuples below)
    data = pd.read_csv(path, sep='\t', usecols=USECOLS)[USECOLS]

    # Update mapping
    for reaction, reagent_1, reagent_2, reagent_3, reagent_4, reaction_type in data.itertuples(index=False):
        for reagent_index, reagent in enumerate([reagent_1, reagent_2, reagent_3, reagent_4]):
            if not np.isnan(reagent):
                reaction_to_reagents[reaction][reagent_index][reaction_type].add(int(reagent))

    # Convert to regular dict for compatibility with multiprocessing
    reaction_to_reagents = {
        reaction: {
            reagent_index: {
                reaction_type: reagent_ids
                for reaction_type, reagent_ids in reaction_type_to_reagent_ids.items()
            }
            for reagent_index, reaction_type_to_reagent_ids in reagent_mapping.items()
        }
        for reaction, reagent_mapping in reaction_to_reagents.items()
    }

    return len(data), reaction_to_reagents


def map_real_reactions_to_reagents(args: Args) -> None:
    """Determines which REAL reagents can be used in which REAL reactions."""
    # Create list of all mappings
    all_reaction_to_reagents = []

    # Get paths to data files
    data_paths = sorted(args.data_dir.rglob('*.cxsmiles.bz2'))

    # Loop through all REAL database files
    with tqdm(total=31000000000, desc='Mapping') as pbar:
        with Pool() as pool:
            for data_size, reaction_to_reagents in pool.imap(map_reactions_for_file, data_paths):
                all_reaction_to_reagents.append(reaction_to_reagents)
                pbar.update(data_size)

    # Merge dictionaries
    combined_reaction_to_reagents = defaultdict(lambda: defaultdict(lambda: defaultdict(set)))
    for reaction_to_reagents in tqdm(all_reaction_to_reagents, desc='Merging'):
        for reaction, reagent_mapping in reaction_to_reagents.items():
            for reagent_index, reaction_type_to_reagent_ids in reagent_mapping.items():
                for reaction_type, reagent_ids in reaction_type_to_reagent_ids.items():
                    combined_reaction_to_reagents[reaction][reagent_index][reaction_type] |= reagent_ids

    # Convert sets to sorted lists for JSON serializability
    combined_reaction_to_reagents = {
        reaction: {
            reagent_index: {
                reaction_type: sorted(reagent_ids)
                for reaction_type, reagent_ids in reaction_type_to_reagent_ids.items()
            }
            for reagent_index, reaction_type_to_reagent_ids in reagent_mapping.items()
        }
        for reaction, reagent_mapping in combined_reaction_to_reagents.items()
    }

    # Save mapping
    with open(args.save_path, 'w') as f:
        json.dump(combined_reaction_to_reagents, f, indent=4, sort_keys=True)


if __name__ == '__main__':
    map_real_reactions_to_reagents(Args().parse_args())
