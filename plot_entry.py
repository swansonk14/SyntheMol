"""Plot eNTRY scores."""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from tap import Tap


class Args(Tap):
    data_paths: list[Path]  # Paths to files containing eNTRy scores.
    names: list[str]  # Names of the datasets
    save_dir: Path  # Path to directory where plots will be saved.

    def process_args(self) -> None:
        self.save_dir.mkdir(parents=True, exist_ok=True)


ENTRY_SCORES = np.array([0, 1, 2, 3])
ZERO_ONE = np.array([0, 1])
ENTRY_PROPS_TO_THRESHOLDS = {
    'primary_amine': 0.5,  # technically only 1 or 0
    'globularity': 0.25,
    'rotatable_bonds': 5
}


def plot_entry(args: Args) -> None:
    """Plot eNTRY scores."""
    if len(args.data_paths) != len(args.names):
        raise ValueError('Number of data paths must match number of names.')

    # Load data
    datas = [pd.read_csv(path) for path in args.data_paths]

    # Set bar width
    width = 0.8 / len(datas)

    # Plot eNTRy scores
    plt.clf()
    for i, (data, name) in enumerate(zip(datas, args.names)):
        score_counts = dict(data['entry_score'].value_counts())
        plt.bar(
            ENTRY_SCORES + i * width,
            [score_counts.get(score, 0) / len(data) for score in ENTRY_SCORES],
            width=width,
            label=name,
        )

    plt.xticks(ENTRY_SCORES + width * (len(datas) - 1) / 2, ENTRY_SCORES)
    plt.xlabel('eNTRy Score')
    plt.ylabel('Fraction of Molecules')
    plt.title('eNTRy Score Proportions')
    plt.legend()
    plt.savefig(args.save_dir / 'entry_scores.pdf', bbox_inches='tight')

    # TODO: Save eNTRy scores

    # Plot eNTRy properties
    for entry_prop, threshold in ENTRY_PROPS_TO_THRESHOLDS.items():
        plt.clf()
        for i, (data, name) in enumerate(zip(datas, args.names)):
            has_prop_percent = len(data[data[entry_prop] <= threshold]) / len(data)

            if entry_prop == 'primary_amine':
                has_prop_percent = 1 - has_prop_percent

            plt.bar(
                ZERO_ONE + i * width,
                [1 - has_prop_percent, has_prop_percent],
                width=width,
                label=name,
            )

        entry_prop_name = entry_prop.replace('_', ' ').title()
        plt.xticks(ZERO_ONE + width * (len(datas) - 1) / 2, ['Fail', 'Pass'])
        plt.xlabel(f'{entry_prop_name} Status')
        plt.ylabel('Fraction of Molecules')
        plt.title(f'{entry_prop_name} Proportions')
        plt.legend()
        plt.savefig(args.save_dir / f'entry_{entry_prop}.pdf', bbox_inches='tight')

    # TODO: Save property scores


if __name__ == '__main__':
    plot_entry(Args().parse_args())
