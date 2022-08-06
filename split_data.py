"""Splits data into train and test splits."""
from pathlib import Path

import pandas as pd
from tap import Tap


class Args(Tap):
    data_path: Path  # Path to CSV file containing data.
    train_save_path: Path  # Path to CSV file where train data will be saved.
    test_save_path: Path  # Path to CSV file where test data will be saved.
    train_size: float = 0.8  # Size of the train set.

    def process_args(self) -> None:
        self.train_save_path.parent.mkdir(parents=True, exist_ok=True)
        self.test_save_path.parent.mkdir(parents=True, exist_ok=True)


def split_data(args: Args) -> None:
    """Splits data into train and test splits."""
    # Load data
    data = pd.read_csv(args.data_path)
    print(f'Data size = {len(data)}')

    # Split data
    train_data = data.sample(frac=args.train_size, random_state=0)
    test_data = data.drop(train_data.index)

    print(f'Train size = {len(train_data)}')
    print(f'Test size = {len(test_data)}')

    # Save data
    train_data.to_csv(args.train_save_path, index=False)
    test_data.to_csv(args.test_save_path, index=False)


if __name__ == '__main__':
    split_data(Args().parse_args())
