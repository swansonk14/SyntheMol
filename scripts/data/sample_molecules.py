"""Sample molecules uniformly at random from a CSV file."""
from pathlib import Path

import pandas as pd


def sample_molecules(data_path: Path, save_path: Path, num_molecules: int) -> None:
    """Sample molecules uniformly at random from a CSV file.

    :param data_path: Path to CSV file containing molecules.
    :param save_path: Path to CSV file where sampled molecules will be saved.
    :param num_molecules: Number of molecules to sample.
    """
    # Load data
    data = pd.read_csv(data_path)

    # Sample molecules
    sampled_data = data.sample(n=num_molecules, random_state=0, replace=False)

    # Save sampled data
    save_path.parent.mkdir(parents=True, exist_ok=True)
    sampled_data.to_csv(save_path, index=False)


if __name__ == "__main__":
    from tap import tapify

    tapify(sample_molecules)
