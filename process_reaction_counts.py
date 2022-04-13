"""Process reactions counts from the REAL database."""
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from tap import Tap


class Args(Tap):
    counts_path: Path  # Path to JSON file with reaction counts output by count_reactions.py.
    save_table_path: Path  # Path to CSV file where the reaction counts tabular data will be saved.
    save_plot_path: Path  # Path to PDF/PNG file where the reaction counts plot will be saved.

    def process_args(self) -> None:
        self.save_table_path.parent.mkdir(parents=True, exist_ok=True)
        self.save_plot_path.parent.mkdir(parents=True, exist_ok=True)


def process_reaction_counts(args: Args) -> None:
    """Process reactions counts from the REAL database."""
    # Load reaction counts data
    with open(args.counts_path) as f:
        reaction_counts = json.load(f)

    # Build DataFrame with reactions and counts
    data = pd.DataFrame(data=[
        {
            'reaction': reaction,
            'count': count

        } for reaction, count in reaction_counts.items()
    ])

    # Sort data by count from largest to smallest
    data.sort_values(by='count', reversed=True)

    # Add cumulative sum and percent
    num_molecules = data['count'].sum()
    data['cumulative_count'] = np.cumsum(reaction_counts)
    data['cumulative_percent'] = data['cumulative_count'] / num_molecules

    # Save data
    data.to_csv(args.save_table_path, index=False)

    # Generate and save plot
    plt.scatter(np.arange(len(data)), data['cumulative_count'])
    plt.xlabel('Reaction Index')
    plt.ylabel('Cumulative Sum of Molecules')
    plt.title('Reaction Counts')
    plt.savefig(args.save_plot_path)


if __name__ == '__main__':
    process_reaction_counts(Args().parse_args())
