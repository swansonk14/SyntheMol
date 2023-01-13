"""Plot score distribution."""
from pathlib import Path

import pandas as pd
from tap import Tap

from assess_generated_molecules import plot_scores


class Args(Tap):
    data_path: Path  # Path to CSV file containing scores.
    score_column: str  # Name of the column containing scores.
    save_dir: Path  # Path to directory where the plot will be saved.

    def process_args(self) -> None:
        self.save_dir.mkdir(parents=True, exist_ok=True)


def plot_score_distribution(args: Args) -> None:
    """Plot score distribution"""
    # Load data
    data = pd.read_csv(args.data_path)

    # Plot score distribution
    plot_scores(scores=data[args.score_column], save_dir=args.save_dir)


if __name__ == '__main__':
    plot_score_distribution(Args().parse_args())
