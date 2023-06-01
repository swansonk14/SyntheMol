"""Plots a heatmap of the MCTS node scores over rollouts."""
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from tqdm import trange

from synthemol.constants import REAL_BUILDING_BLOCK_ID_COL, SCORE_COL


def plot_heatmap(
        data_path: Path,
        building_blocks_path: Path,
        save_dir: Path,
        num_reactions: int,
        building_blocks_id_column: str = REAL_BUILDING_BLOCK_ID_COL,
        building_blocks_score_column: str = SCORE_COL,
) -> None:
    """Plots a heatmap of the MCTS node scores over rollouts.

    :param data_path: Path to a CSV file containing MCTS results.
    :param building_blocks_path: Path to a CSV file containing building block IDs, SMILES, and model scores.
    :param save_dir: Path to directory where the plot will be saved.
    :param num_reactions: Number of reactions in the final molecule for which to plot the heatmap.
    :param building_blocks_id_column: Name of the column in the building blocks file containing building block IDs.
    :param building_blocks_score_column: Name of the column in the building blocks file containing building block scores.
    """
    # TODO: change MCTS saving to save a CSV file with the node scores for each rollout

    # Load data
    data = pd.read_csv(data_path)
    building_blocks = pd.read_csv(building_blocks_path)

    # Map from building block ID to score
    building_block_id_to_score = dict(zip(
        building_blocks[building_blocks_id_column],
        building_blocks[building_blocks_score_column]
    ))

    # Map rollout to row indices
    rollout_to_row_index = defaultdict(list)
    for row_index, rollout in enumerate(data['rollout_num']):
        rollout_to_row_index[rollout].append(row_index)

    # Build heatmap with rollouts as columns and nodes as rows
    building_block_id_columns = sorted(
        column
        for column in data.columns
        if column.startswith('building_block_') and column.endswith('_id') and int(column.split('_')[2]) <= num_reactions
    )
    max_rollout = max(data['rollout_num'])

    row_size = 100
    heatmap = np.zeros((row_size * (len(building_block_id_columns) + 1), len(data)))

    molecule_index = 0
    for rollout in trange(1, max_rollout + 1):
        for row_index in rollout_to_row_index[rollout]:
            # Skip if wrong number of reactions
            if data.loc[row_index, 'num_reactions'] != num_reactions:
                continue

            # Get building block scores
            for i, building_block_id_column in enumerate(building_block_id_columns):
                heatmap[i * row_size:(i + 1) * row_size, molecule_index] = building_block_id_to_score.get(
                    data.loc[row_index, building_block_id_column], 0
                )

            # Get full molecule score
            heatmap[-row_size:, molecule_index] = data.loc[row_index, 'score']
            molecule_index += 1

    heatmap = heatmap[:, :molecule_index]

    # Plot and save heatmap
    save_dir.mkdir(parents=True, exist_ok=True)

    num_im_rows = 10
    fig, axes = plt.subplots(num_im_rows, 1)

    # Split the rows of heatmap into equal parts
    heatmap_split = np.array_split(heatmap, num_im_rows, axis=1)

    for heatmap_section, ax in zip(heatmap_split, axes):
        ax.imshow(heatmap_section, cmap='hot', interpolation='none')
        ax.axis('off')

    # fig.colorbar(heatmap)
    # plt.xlabel('Molecule')
    # plt.ylabel('Node')
    # plt.title('MCTS Node Scores')

    plt.savefig(save_dir / f'mcts_node_scores_{num_reactions}_reactions.pdf', bbox_inches='tight')


if __name__ == '__main__':
    from tap import tapify

    tapify(plot_heatmap)
