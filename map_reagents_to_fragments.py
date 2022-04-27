"""Maps REAL reagents to fragments that match those reagents."""
import json
from pathlib import Path

import pandas as pd
from rdkit import Chem
from tap import Tap
from tqdm import tqdm

from real_reactions import REAL_REACTIONS


class Args(Tap):
    fragment_path: Path  # Path to CSV file containing REAL fragments.
    smiles_column: str = 'smiles'  # Name of the column containing SMILES.
    save_path: Path  # Path to JSON file where a dictionary mapping reagents to fragments will be saved.

    def process_args(self) -> None:
        self.save_path.parent.mkdir(parents=True, exist_ok=True)


def map_reagents_to_fragments(args: Args) -> None:
    """Maps REAL reagents to fragments that match those reagents."""
    # Load data
    fragment_data = pd.read_csv(args.fragment_path)

    # Get fragment SMILES
    fragments = sorted(set(fragment_data[args.smiles_column]))

    # Convert fragment SMILES to mols
    fragment_mols = [Chem.AddHs(Chem.MolFromSmiles(fragment)) for fragment in tqdm(fragments, desc='SMILES to mol')]

    # Map reagents to fragments
    reagent_to_fragments = {}

    for reaction in tqdm(REAL_REACTIONS, desc='Mapping reagents to fragments'):
        for reagent in reaction.reagents:
            if reagent.smarts in reagent_to_fragments:
                continue

            reagent_to_fragments[reagent.smarts] = [
                fragment
                for fragment, fragment_mol in tqdm(zip(fragments, fragment_mols), total=len(fragments), leave=False)
                if reagent.has_substruct_match(fragment_mol)
            ]

    # Save data
    with open(args.save_path, 'w') as f:
        json.dump(reagent_to_fragments, f, indent=4, sort_keys=True)


if __name__ == '__main__':
    map_reagents_to_fragments(Args().parse_args())
