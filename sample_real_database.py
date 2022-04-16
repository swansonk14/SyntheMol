"""Samples molecules from the REAL database."""
from pathlib import Path

import pandas as pd
from tap import Tap
from tqdm import tqdm


class Args(Tap):
    data_dir: Path  # Path to directory with CXSMILES files containing the REAL database.
    num_molecules: int  # Number of molecules to sample per CXSMILES file.
    save_dir: Path  # Path to directory where sampled molecules will be saved.

    def process_args(self) -> None:
        self.save_dir.mkdir(parents=True, exist_ok=True)


def sample_real_database(args: Args) -> None:
    """Samples molecules from the REAL database."""
    # Loop through all REAL database files
    for path in tqdm(list(args.data_dir.glob('*.cxsmiles')), desc='Sampling'):
        # Load REAL data file
        data = pd.read_csv(path, sep='\t')

        # Sample REAL data
        data = data.sample(n=args.num_molecules, random_state=0)

        # Save sampled data
        data.to_csv(args.save_dir / f'{path.stem}_sampled.cxsmiles', index=False)


if __name__ == '__main__':
    sample_real_database(Args().parse_args())
