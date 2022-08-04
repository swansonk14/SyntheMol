"""Map fragments to train similarities."""
import json
from pathlib import Path
from time import time
from typing import Literal

import pandas as pd
from tap import Tap

from chem_utils.molecular_similarities import compute_top_similarities


class Args(Tap):
    fragment_path: Path  # Path to a CSV file containing fragments.
    train_path: Path  # Path to CSV file containing SMILES for the active molecules in the training set.
    save_path: Path  # Path to a JSON file where a dictionary mapping fragments to train similarities will be saved.
    similarity_type: Literal['tanimoto', 'tversky']  # Type of similarity to compute.
    top_k: int = 1  # Computes the similarity of the top_k most similarity training molecule.
    smiles_column: str = 'smiles'  # Name of the column containing SMILES.

    def process_args(self) -> None:
        self.save_path.parent.mkdir(parents=True, exist_ok=True)


def map_fragments_to_scores(args: Args) -> None:
    """Map fragments to prediction scores."""
    # Load fragments
    fragments = sorted(set(pd.read_csv(args.fragment_path)[args.smiles_column]))

    # Load train hits
    train = pd.read_csv(args.train_path)

    # Map fragments to train similarity
    start_time = time()
    train_similarities = compute_top_similarities(
        similarity_type=args.similarity_type,
        mols=fragments,
        reference_mols=train[args.smiles_column],
        top_k=args.top_k
    )
    fragment_to_train_similarity = dict(zip(fragments, train_similarities.tolist()))
    print(f'Total time to map fragments to train similarities = {time() - start_time:.2f}')

    # Save fragment to train similarity mapping
    with open(args.save_path, 'w') as f:
        json.dump(fragment_to_train_similarity, f, indent=4, sort_keys=True)


if __name__ == '__main__':
    map_fragments_to_scores(Args().parse_args())
