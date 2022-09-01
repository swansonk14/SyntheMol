"""Determines which REAL reagents can be used in which REAL reactions."""
import json
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


def map_real_reactions_to_reagents(args: Args) -> None:
    """Determines which REAL reagents can be used in which REAL reactions."""
    # Create mapping from reaction ID to reagent number to reaction type to valid reagent IDs
    reaction_to_reagent_to_ids: dict[int, dict[int, dict[str, set[int]]]] = \
        defaultdict(lambda: defaultdict(lambda: defaultdict(set)))

    # Loop through all REAL database files
    for path in tqdm(list(args.data_dir.rglob('*.cxsmiles.bz2'))):
        # Load REAL data file (ensures cols are in the order of USECOLS for itertuples below)
        data = pd.read_csv(path, sep='\t', usecols=USECOLS)[USECOLS]

        # Update mapping
        file_name = path.stem.split('.')[0]
        for reaction, reagent_1, reagent_2, reagent_3, reagent_4, reaction_type in tqdm(data.itertuples(index=False),
                                                                                        total=len(data), leave=False,
                                                                                        desc=file_name):
            for reagent_index, reagent in enumerate([reagent_1, reagent_2, reagent_3, reagent_4]):
                if not np.isnan(reagent):
                    reaction_to_reagent_to_ids[reaction][reagent_index][reaction_type].add(int(reagent))

    # Convert to JSON serializable
    reaction_to_reagent_to_ids = {
        reaction: {
            reagent_index: {
                reaction_type: sorted(reagent_ids)
                for reaction_type, reagent_ids in reaction_type_to_reagent_ids.items()
            }
            for reagent_index, reaction_type_to_reagent_ids in reagent_mapping.items()
        }
        for reaction, reagent_mapping in reaction_to_reagent_to_ids.items()
    }

    # Save mapping
    with open(args.save_path, 'w') as f:
        json.dump(reaction_to_reagent_to_ids, f, indent=4, sort_keys=True)


if __name__ == '__main__':
    map_real_reactions_to_reagents(Args().parse_args())
