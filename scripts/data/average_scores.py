"""Computes a weighted average of the scores of multiple models in a CSV file."""
from pathlib import Path

import pandas as pd


def average_scores(
    data_path: Path,
    score_columns: list[str],
    score_weights: list[float],
    average_column: str = "average_score",
    save_path: Path | None = None,
) -> None:
    """Computes a weighted average of the scores of multiple models in a CSV file.

    :param data_path: Path to a CSV file containing scores.
    :param score_columns: A list of columns containing scores.
    :param score_weights: A list of weights for each score column.
    :param average_column: Name of the column containing the average score.
    :param save_path: Path to a CSV file where the average scores will be saved. If None, uses data_path.
    """
    # Load data
    data = pd.read_csv(data_path)

    # Compute weighted average of scores
    data[average_column] = sum(
        weight * data[score_column]
        for score_column, weight in zip(score_columns, score_weights)
    )

    # Save data
    if save_path is None:
        save_path = data_path

    save_path.parent.mkdir(parents=True, exist_ok=True)
    data.to_csv(save_path, index=False)


if __name__ == "__main__":
    from tap import tapify

    tapify(average_scores)
