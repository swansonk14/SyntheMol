"""Plots the MCTS results over time."""
from pathlib import Path
from typing import Literal, Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from tap import Tap


def plot_mcts_over_time(
        data_path: Path,
        save_dir: Path,
        model_name: str,
        plot_type: Literal['histogram', 'line', 'violin'] = 'histogram',
        increment: int = 50000,
        min_score: Optional[float] = None,
) -> None:
    """Plots the MCTS results over time.

    :param data_path: Path to CSV file containing MCTS generated molecules.
    :param save_dir: Path to directory where the plot will be saved.
    :param model_name: The name of the model used during MCTS.
    :param plot_type: The type of plot to generate.
    :param increment: The number of rollouts between each plot.
    :param min_score: If provided, only molecules with scores >= this threshold are plotted.
    """
    # Load data
    data = pd.read_csv(data_path)

    # Optionally threshold by score
    if min_score is not None:
        data = data[data['score'] >= min_score]

    # Plot MCTS results over time
    plt.clf()
    num_rollouts = data['rollout_num'].max()
    alpha = increment / num_rollouts

    if plot_type == 'histogram':
        # Histograms of scores by rollout bin
        for rollout_num in range(0, num_rollouts, increment):
            rollout_data = data[(data['rollout_num'] > rollout_num) & (data['rollout_num'] <= rollout_num + increment)]
            plt.hist(rollout_data['score'], bins=100, alpha=alpha,
                     label=f'Rollouts {rollout_num + 1:,} - {rollout_num + increment:,}')

        plt.legend()
        plt.xlabel('Score')
        plt.ylabel('Count')
    elif plot_type == 'line':
        # Line plot of mean and standard deviation score by rollout bin
        rollouts, means, stds = [], [], []

        for rollout_num in range(0, num_rollouts, increment):
            rollout_data = data[(data['rollout_num'] > rollout_num) & (data['rollout_num'] <= rollout_num + increment)]
            rollouts.append(rollout_num)
            means.append(np.mean(rollout_data['score']))
            stds.append(np.std(rollout_data['score']))

        plt.errorbar(rollouts, means, yerr=stds, fmt='-o', capsize=5)
        plt.xlabel('Rollout')
        plt.ylabel('Score')
    elif plot_type == 'violin':
        # Violin plot of scores by rollout bin
        rollouts, scores = [], []

        for rollout_num in range(0, num_rollouts, increment):
            rollout_data = data[(data['rollout_num'] > rollout_num) & (data['rollout_num'] <= rollout_num + increment)]
            rollouts.append(rollout_num)
            scores.append(rollout_data['score'])

        plt.violinplot(scores, rollouts, widths=0.95 * increment, showmedians=True)
        plt.xlabel('Rollout')
        plt.ylabel('Score')
    else:
        raise ValueError(f'Invalid plot type: {plot_type}')

    # Save plot
    plt.title(f'{model_name} MCTS Score Over Rollouts')
    save_dir.mkdir(parents=True, exist_ok=True)
    plt.savefig(save_dir / 'mcts_rollout_scores.pdf', bbox_inches='tight')

    # Save data
    fig_data = data[['rollout_num', 'score']].sort_values('rollout_num')
    fig_data.to_csv(save_dir / 'mcts_rollout_scores.csv', index=False)


if __name__ == '__main__':
    class Args(Tap):
        data_path: Path  # Path to CSV file containing MCTS generated molecules.
        save_dir: Path  # Path to directory where the plot will be saved.
        model_name: str  # The name of the model used during MCTS.
        plot_type: Literal['histogram', 'line', 'violin'] = 'histogram'  # The type of plot to generate.
        increment: int = 50000  # The number of rollouts between each plot.
        min_score: Optional[float] = None  # If provided, only molecules with scores >= this threshold are plotted.

    plot_mcts_over_time(**Args().parse_args().as_dict())
