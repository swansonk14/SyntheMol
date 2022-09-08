"""Average scores across multiple columns of a CSV file."""
from pathlib import Path
from typing import Optional

import pandas as pd
from tap import Tap


class Args(Tap):
    data_path: Path  # Path to a CSV file containing multiple score columns.
    score_columns: list[str]  # Names of the columns containing scores to average.
    average_column: str = 'average_score'  # Name of the column that will contain the average score.
    save_path: Optional[Path] = None  # Path to a CSV file where the data with averaged scores will be saved. If None, defaults to data_path.

    def process_args(self) -> None:
        if self.save_path is None:
            self.save_path = self.data_path

        self.save_path.parent.mkdir(parents=True, exist_ok=True)


def average_scores(args: Args) -> None:
    """Average scores across multiple columns of a CSV file."""
    # Load data
    data = pd.read_csv(args.data_path)

    # Average scores
    data[args.average_column] = data[args.score_columns].mean(axis=1)

    # Save data
    data.to_csv(args.save_path, index=False)


if __name__ == '__main__':
    average_scores(Args().parse_args())
