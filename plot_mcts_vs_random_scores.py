"""Plot the scores of MCTS vs randomly generated molecules."""
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from tap import Tap


class Args(Tap):
    random_path: Path  # Path to CSV file containing randomly generated molecules with model scores.
    mcts_path: Path  # Path to CSV file containing MCTS molecules with model scores.
    mcts_name: str  # Name of MCTS model.
    save_path: Path  # Path to PDF/PNG file where the plot will be saved.

    def process_args(self) -> None:
        self.save_path.parent.mkdir(parents=True, exist_ok=True)


def plot_mcts_vs_random_scores(args: Args) -> None:
    """Plot the scores of MCTS vs randomly generated molecules."""
    # Load data
    random = pd.read_csv(args.random_path)
    mcts = pd.read_csv(args.mcts_path)

    # Plot histogram of scores of random and mcts molecules
    plt.clf()
    plt.hist(random[f'{args.mcts_name}_ensemble_preds'], bins=100, alpha=0.5, label='random')
    plt.hist(mcts['score'], bins=100, alpha=0.5, label=f'{args.mcts_name}')
    plt.xlabel('Score')
    plt.ylabel('Count')
    plt.title(f'Random vs MCTS with {args.mcts_name} Model')
    plt.legend()
    plt.savefig(args.save_path, bbox_inches='tight')


if __name__ == '__main__':
    plot_mcts_vs_random_scores(Args().parse_args())
