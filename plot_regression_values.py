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

    # Compute mean value
    data['mean'] = (data[args.rep1_column] + data[args.rep2_column]) / 2
    data.sort_values(by='mean', inplace=True)

    # Optionally remove outliers (5 points with largest mean inhibition)
    if args.remove_outliers:
        data = data.iloc[:-5]

    # Compute binary activity threshold
    mean = data['mean'].mean()
    std = data['mean'].std()
    threshold = mean - 2 * std
    threshold_color = ['red']
    threshold_name = 'mean - 2 * std'

    # Get regression values
    rep1 = data[args.rep1_column]
    rep2 = data[args.rep2_column]
    index = np.arange(len(data))

    # Get data name
    data_name = args.data_path.stem

    # Plot r1 and r2
    for rep_num, rep in [(1, rep1), (2, rep2)]:
        plt.clf()
        plt.scatter(index, sorted(rep), s=5)

        plt.hlines(threshold, 0, len(rep), colors=threshold_color, label=threshold_name)

        plt.xlabel('Molecule Index')
        plt.ylabel('Inhibition')
        plt.title(f'{data_name} Replicate {rep_num} Inhibition')
        plt.legend()
        plt.savefig(args.save_dir / f'replicate_{rep_num}.pdf')

    # Plot mean (sorted)
    plt.clf()
    plt.scatter(index, sorted(data['mean']), s=5)

    plt.hlines(threshold, 0, len(data), colors=threshold_color, label=threshold_name)

    plt.xlabel('Molecule Index')
    plt.ylabel('Mean Inhibition')
    plt.title(f'{data_name} Mean Inhibition')
    plt.legend()
    plt.savefig(args.save_dir / 'mean_sorted.pdf')

    # Plot mean (unsorted)
    mean_values = np.array(sorted(data['mean']))
    np.random.seed(0)
    np.random.shuffle(mean_values)

    plt.clf()
    plt.scatter(index, mean_values, s=5)

    plt.hlines(threshold, 0, len(data), colors=threshold_color, label=threshold_name)

    plt.xlabel('Molecule Index')
    plt.ylabel('Mean Inhibition')
    plt.title(f'{data_name} Mean Inhibition')
    plt.legend()
    plt.savefig(args.save_dir / 'mean_unsorted.pdf')

    # Plot r1 vs r2
    plt.clf()
    plt.scatter(rep1, rep2, s=5)

    plt.hlines(threshold, 0, threshold, colors=threshold_color, label=threshold_name)
    plt.vlines(threshold, 0, threshold, colors=threshold_color)

    plt.xlabel('Replicate 1 Inhibition')
    plt.ylabel('Replicate 2 Inhibition')
    plt.title(f'{data_name} Replicate 1 vs 2 Inhibition')
    plt.legend()
    plt.savefig(args.save_dir / 'replicate_1_vs_2.pdf')


if __name__ == '__main__':
    plot_regression_values(Args().parse_args())
