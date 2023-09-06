"""Plots the MCTS results over time."""
from pathlib import Path
from typing import Literal

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from synthemol.constants import ROLLOUT_COL, SCORE_COL


def plot_mcts_over_time(
        data_path: Path,
        save_dir: Path,
        model_name: str,
        plot_type: Literal['histogram', 'line', 'violin'] = 'violin',
        increment: int = 50000,
        min_score: float | None = None,
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
        data = data[data[SCORE_COL] >= min_score]

    # Bin rollouts
    num_rollouts = data[ROLLOUT_COL].max()
    rollout_bins, rollout_bin_scores = [], []

    for rollout_num in range(0, num_rollouts, increment):
        rollout_data = data[(data[ROLLOUT_COL] > rollout_num) & (data[ROLLOUT_COL] <= rollout_num + increment)]
        rollout_bins.append((rollout_num, rollout_num + increment))
        rollout_bin_scores.append(rollout_data[SCORE_COL])

    rollout_bin_labels = [
        f'{rollout_bin_start + 1:,} - {rollout_bin_end:,}' for
        rollout_bin_start, rollout_bin_end in rollout_bins
    ]

    # Plot MCTS results over time
    plt.clf()
    alpha = increment / num_rollouts

    if plot_type == 'histogram':
        # Histograms of scores by rollout bin
        for (rollout_bin_start, rollout_bin_end), rollout_bin_score in zip(rollout_bins, rollout_bin_scores):
            plt.hist(rollout_bin_scores, bins=100, alpha=alpha,
                     label=f'Rollouts {rollout_bin_start + 1:,} - {rollout_bin_end:,}')

        plt.legend()
        plt.xlabel('Score')
        plt.ylabel('Count')
    elif plot_type == 'line':
        # Line plot of mean and standard deviation score by rollout bin
        means, stds = [], []

        for rollout_bin_score in rollout_bin_scores:
            means.append(np.mean(rollout_bin_score))
            stds.append(np.std(rollout_bin_score))

        xticks = np.arange(len(means))
        plt.errorbar(xticks, means, yerr=stds, fmt='-o', capsize=5)
        plt.xticks(xticks, rollout_bin_labels)
        plt.xlabel('Rollout')
        plt.ylabel('Score')
    elif plot_type == 'violin':
        # Violin plot of scores by rollout bin
        xticks = np.arange(len(rollout_bin_scores))
        plt.violinplot(rollout_bin_scores, xticks, widths=0.95, showmedians=True)
        plt.xticks(xticks, rollout_bin_labels, rotation=45)
        plt.xlabel('Rollout')
        plt.ylabel('Score')
    else:
        raise ValueError(f'Invalid plot type: {plot_type}')

    # Save plot
    plt.title(f'{model_name} MCTS Score Over Rollouts')
    save_dir.mkdir(parents=True, exist_ok=True)
    plt.savefig(save_dir / 'mcts_rollout_scores.pdf', bbox_inches='tight')

    # Save data
    max_len = max(len(rollout_bin_score) for rollout_bin_score in rollout_bin_scores)
    fig_data = pd.DataFrame({
        f'rollouts_{rollout_bin_start + 1}-{rollout_bin_end}': np.pad(
            rollout_bin_score,
            (0, max_len - len(rollout_bin_score)),
            constant_values=np.nan
        )
        for (rollout_bin_start, rollout_bin_end), rollout_bin_score in zip(rollout_bins, rollout_bin_scores)
    })
    fig_data.to_csv(save_dir / 'mcts_rollout_scores.csv', index=False)


if __name__ == '__main__':
    from tap import tapify

    tapify(plot_mcts_over_time)
