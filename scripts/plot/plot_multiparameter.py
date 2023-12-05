"""Plot the distribution of two scores for a set of molecules."""
from collections import Counter
from pathlib import Path

import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt


def plot_multiparameter(
    data_paths: list[Path],
    data_names: list[str],
    save_path: Path,
    property_column_x: str,
    property_threshold_x: float,
    property_column_y: str,
    property_threshold_y: float,
) -> None:
    """Plot the distribution of two scores for a set of molecules.

    :param data_paths: Path to CSV files containing scores.
    :param data_names: Names of the datasets.
    :param save_path: Path to a PDF file where the plot will be saved.
    :param property_column_x: Name of the column containing the x-axis property.
    :param property_threshold_x: Threshold for the x-axis property.
    :param property_column_y: Name of the column containing the y-axis property.
    :param property_threshold_y: Threshold for the y-axis property.
    """
    # Check lengths
    if len(data_paths) != len(data_names):
        raise ValueError(
            f"Number of data paths ({len(data_paths):,}) must equal number of data names ({len(data_names):,})"
        )

    # Load data
    data_dfs = [pd.read_csv(data_path) for data_path in data_paths]

    # Initialize mapping from name to quadrant counts
    name_to_quadrant_counts = {}

    # Get seaborn default colors
    colors = sns.color_palette()

    # Plot scores
    for data, name, color in zip(data_dfs, data_names, colors):
        # Plot
        sns.scatterplot(
            x=data[property_column_x],
            y=data[property_column_y],
            label=name,
            s=1,
            color=color,
            edgecolor=None,
            alpha=1 / len(data_paths),
        )

        # Count points in each quadrant
        quadrants = []
        for x, y in zip(data[property_column_x], data[property_column_y]):
            if x >= property_threshold_x and y >= property_threshold_y:
                quadrants.append("Q1")
            elif x < property_threshold_x and y >= property_threshold_y:
                quadrants.append("Q2")
            elif x < property_threshold_x and y < property_threshold_y:
                quadrants.append("Q3")
            else:
                quadrants.append("Q4")

        # Add quadrant counts
        name_to_quadrant_counts[name] = Counter(quadrants)

    # Plot quadrant counts
    quadrant_to_location = {
        "Q1": (0.95, 0.95),
        "Q2": (0.15, 0.95),
        "Q3": (0.15, 0.15),
        "Q4": (0.95, 0.15),
    }
    for quadrant, (x, y) in quadrant_to_location.items():
        for i, (name, counts) in enumerate(name_to_quadrant_counts.items()):
            plt.text(
                x=x,
                y=y - i * 0.04,
                s=f"{name}: {counts[quadrant]:,}",
                transform=plt.gca().transAxes,
                color=colors[i],
                ha="right" if x > 0.5 else "left",
                va="top" if y > 0.5 else "bottom",
            )

    # Plot thresholds
    plt.axvline(property_threshold_x, color="red", linestyle="--")
    plt.axhline(property_threshold_y, color="red", linestyle="--")

    # Turn off legend
    plt.legend().remove()

    # Save plot
    save_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(save_path, bbox_inches="tight")


if __name__ == "__main__":
    from tap import tapify

    tapify(plot_multiparameter)
