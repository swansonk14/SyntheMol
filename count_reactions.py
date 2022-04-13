"""Count reactions in the REAL database."""
import json
from collections import Counter
from pathlib import Path

import pandas as pd
from tap import Tap
from tqdm import tqdm


class Args(Tap):
    data_dir: Path  # Path to directory with CXSMILES files containing the REAL database.
    save_path: Path  # Path to JSON file where results will be saved.

    def process_args(self) -> None:
        self.save_path.parent.mkdir(parents=True, exist_ok=True)


def count_reactions(args: Args) -> None:
    """Count reactions in the REAL database."""
    # Get files
    paths = list(args.data_dir.glob('*.cxsmiles'))

    # Count reactions
    reaction_counts = Counter()

    for path in tqdm(paths):
        reactions = pd.read_csv(path, sep='\t', usecols=['reaction'], dtype={'reaction': int})
        reaction_counts.update(reactions['reaction'])

    # Count total molecules
    num_molecules = sum(reaction_counts.values())
    print(f'Number of molecules = {num_molecules:,}')
    print(f'Number of reactions = {len(reaction_counts):,}')

    # Save reaction counts
    with open(args.save_path, 'w') as f:
        json.dump(dict(reaction_counts), f, indent=4, sort_keys=True)


if __name__ == '__main__':
    count_reactions(Args().parse_args())
