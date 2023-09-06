"""Plot score vs similarity for generated molecules."""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def plot_score_vs_similarity(
    data_dir: Path,
    save_dir: Path,
) -> None:
    """Plot score vs similarity for generated molecules.

    :param data_dir: Path to a directory containing directories with CSV files containing score vs similarity data.
    :param save_dir: Path to a directory where the plots will be saved.
    """
    # Load results
    mcts_paths = list(data_dir.glob("mcts_*/score_vs_similarity.csv"))
    rl_rdkit_paths = list(data_dir.glob("rl_rdkit_*/score_vs_similarity.csv"))
    rl_chemprop_paths = list(data_dir.glob("rl_chemprop_*/score_vs_similarity.csv"))

    # Get score thresholds
    data = pd.read_csv(mcts_paths[0])
    score_thresholds = data["score_threshold"].tolist()

    # Plot score vs similarity for each score threshold
    for score_threshold in sorted(score_thresholds):
        plt.clf()
        plt.subplots()

        # Plot each method
        for paths, marker, label, cmap in [
            (mcts_paths, "o", "mcts", "spring"),
            (rl_rdkit_paths, "x", "rl_rdkit", "winter"),
            (rl_chemprop_paths, "*", "rl_chemprop", "winter"),
        ]:
            # Get relevant statistics
            num_hits_novel = []
            similarity_novel = []
            explore = []
            for path in paths:
                # Load data
                data = pd.read_csv(path)

                # Select results matching score threshold
                data = data[data["score_threshold"] == score_threshold]

                # Collect statistics
                num_hits_novel.append(data["num_hits_novel"].iloc[0])
                similarity_novel.append(data["max_independent_set_novel"].iloc[0])
                explore.append(float(path.parent.name.split("_")[-1]))

            # Plot data
            plt.scatter(
                num_hits_novel,
                similarity_novel,
                marker=marker,
                c=np.log(explore),
                cmap=cmap,
                label=label,
            )

            # Set up colorbar
            if label == "mcts":
                cbar = plt.colorbar()
                cbar.set_label("log(MCTS Explore)")
            elif label == "rl_rdkit":
                cbar = plt.colorbar()
                cbar.set_label("log(RL Temperature)")

        # Label plot
        plt.legend()
        plt.xlabel("Number of Novel Hits")
        plt.ylabel("Max Independent Set of Novel Hits")
        plt.title(f"Hits vs Similarity for Score Threshold {score_threshold}")

        # Save plot
        save_dir.mkdir(parents=True, exist_ok=True)
        plt.savefig(
            save_dir / f"score_threshold_{score_threshold}_novel_max_independent.pdf",
            bbox_inches="tight",
        )


if __name__ == '__main__':
    from tap import tapify

    tapify(plot_score_vs_similarity)
