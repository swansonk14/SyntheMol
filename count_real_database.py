"""Count reactions and reagents in the REAL database."""
import json
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from tap import Tap
from tqdm import tqdm


class Args(Tap):
    data_dir: Path  # Path to directory with CXSMILES files containing the REAL database.
    save_dir: Path  # Path to dreictoy where results will be saved.

    def process_args(self) -> None:
        self.save_dir.mkdir(parents=True, exist_ok=True)


REACTION_COL = 'reaction'
REAGENT_COLS = ['reagent1', 'reagent2', 'reagent3', 'reagent4']
USECOLS = [REACTION_COL] + REAGENT_COLS


def save_counts_as_json(counts: dict[int, int], count_name: str, save_dir: Path) -> None:
    """Save counts as JSON."""
    with open(save_dir / f'{count_name.lower()}_counts.json', 'w') as f:
        json.dump(counts, f, indent=4, sort_keys=True)


def save_counts_as_csv(counts: dict[int, int], count_name: str, save_dir: Path) -> pd.DataFrame:
    """Save counts as CSV with additional statistics."""
    # Build DataFrame with counts
    counts_data = pd.DataFrame(data=[
        {
            count_name: item,
            'count': count

        } for item, count in counts.items()
    ])

    # Sort data by count from largest to smallest
    counts_data.sort_values(by='count', ascending=False, inplace=True)

    # Add percent, cumulative sum, and cumulative percent
    num_molecules = counts_data['count'].sum()

    counts_data['percent'] = counts_data['count'] / num_molecules
    counts_data['cumulative_count'] = np.cumsum(counts_data['count'])
    counts_data['cumulative_percent'] = counts_data['cumulative_count'] / num_molecules

    # Save counts
    counts_data.to_csv(save_dir / f'{count_name.lower()}_counts.csv', index=False)

    return counts_data


def plot_counts(counts_data: pd.DataFrame, count_name: str, save_dir: Path) -> None:
    """Plot counts and save to PDF."""
    # Generate plot of counts
    plt.clf()
    plt.scatter(np.arange(len(counts_data)), counts_data['cumulative_count'])
    plt.gca().set_ylim(bottom=0)
    plt.xlabel(f'{count_name} Index')
    plt.ylabel('Cumulative Sum of Molecules')
    plt.title(f'{count_name} Counts')

    # Save plot
    plt.savefig(save_dir / f'{count_name.lower()}_counts.pdf')


def save_counts(counts: dict[int, int], count_name: str, save_dir: Path) -> None:
    """Saves counts as JSON, CSV, and PDF plot."""
    # Save counts as JSON
    save_counts_as_json(counts=counts, count_name=count_name, save_dir=save_dir)

    # Save counts as CSV and get back DataFrame
    counts_data = save_counts_as_csv(counts=counts, count_name=count_name, save_dir=save_dir)

    # Plot counts and save as PDF
    plot_counts(counts_data=counts_data, count_name=count_name, save_dir=save_dir)


def count_real_database(args: Args) -> None:
    """Count reactions and reagents in the REAL database."""
    # Count reactions and reagents
    reaction_counts = Counter()
    reagent_counts = Counter()
    reagent_counts_by_reaction = defaultdict(Counter)

    # Loop through all REAL database files
    for path in tqdm(list(args.data_dir.glob('*.cxsmiles'))):
        # Load REAL data file
        data = pd.read_csv(path, sep='\t', usecols=USECOLS)

        # Update reaction counts
        reaction_counts.update(data['reaction'])

        # Update reagent counts
        for reagent_col in REAGENT_COLS:
            reagent_counts.update(data[reagent_col].dropna())

        # Update reagent counts by reaction
        for reaction, reagent_1, reagent_2, reagent_3, reagent_4 in data.itertuples(index=False):
            reagents = [reagent_1, reagent_2, reagent_3, reagent_4]
            reagents = [reagent for reagent in reagents if not np.isnan(reagent)]
            reagent_counts_by_reaction[reaction].update(reagents)

    # Save reaction counts
    save_counts(counts=dict(reaction_counts), count_name='Reaction', save_dir=args.save_dir)

    # Save reagent counts
    save_counts(counts=dict(reagent_counts), count_name='Reagent', save_dir=args.save_dir)

    # Save reagent counts by reaction
    for reaction, reagent_counts_for_reaction in reagent_counts_by_reaction.items():
        # Create save directory for this reaction
        save_dir = args.save_dir / 'reagent_counts_by_reaction' / str(reaction)
        save_dir.mkdir(parents=True, exist_ok=True)

        # Save reagent counts for this reaction
        save_counts(counts=dict(reagent_counts_for_reaction), count_name='Reagent', save_dir=save_dir)


if __name__ == '__main__':
    count_real_database(Args().parse_args())
