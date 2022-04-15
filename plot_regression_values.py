"""Plot the raw regression values of the antibiotic inhibition data."""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from tap import Tap


class Args(Tap):
    data_path: Path  # Path to a CSV file containing antibiotic inhibition data.
    rep1_column: str  # Name of the column containing the raw regression values from the first replicate.
    rep2_column: str  # Name of the column containing the raw regression values from the second replicate.
    remove_outliers: bool = False  # Whether to remove outliers (5 points with largest mean inhibition).
    save_dir: Path  # Path to a directory where the plots will be saved.

    def process_args(self) -> None:
        self.save_dir.mkdir(parents=True, exist_ok=True)


def plot_regression_values(args: Args) -> None:
    """Plot the raw regression values of the antibiotic inhibition data."""
    # Load data
    data = pd.read_csv(args.data_path)

    # Optionally remove outliers (5 points with largest mean inhibition)
    if args.remove_outliers:
        data['mean'] = (data[args.rep1_column] + data[args.rep2_column]) / 2
        data.sort_values(by='mean', inplace=True)
        data = data.iloc[:-5]

    # Get regression values
    rep1 = data[args.rep1_column]
    rep2 = data[args.rep2_column]
    index = np.arange(len(data))

    # Plot r1 and r2
    for rep_num, rep in [(1, rep1), (2, rep2)]:
        plt.clf()
        plt.scatter(index, sorted(rep), s=5)
        plt.xlabel('Molecule Index')
        plt.ylabel('Inhibition')
        plt.title(f'Replicate {rep_num} Inhibition')
        plt.savefig(args.save_dir / f'replicate_{rep_num}.pdf')

    # Plot r1 vs r2
    plt.clf()
    plt.scatter(rep1, rep2, s=5)
    plt.xlabel('Replicate 1 Inhibition')
    plt.ylabel('Replicate 2 Inhibition')
    plt.title('Replicate 1 vs 2 Inhibition')
    plt.savefig(args.save_dir / 'replicate_1_vs_2.pdf')


if __name__ == '__main__':
    plot_regression_values(Args().parse_args())
