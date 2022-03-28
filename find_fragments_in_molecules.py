"""Finds fragments in molecules."""
from pathlib import Path

import pandas as pd
from rdkit import Chem
from rdkit.Chem.FragmentMatcher import FragmentMatcher
from tap import Tap
from tqdm import tqdm


class Args(Tap):
    fragment_path: Path  # Path to CSV file containing molecular fragments.
    molecule_path: Path  # Path to CSV file containing full molecules.
    smiles_column: str = 'smiles'  # Name of the column in molecule_path containing SMILES.
    save_path: Path  # Path to CSV file where fragment indicators for each molecule will be saved.

    def process_args(self) -> None:
        self.save_path.parent.mkdir(parents=True, exist_ok=True)


def find_fragments_in_molecules(args: Args) -> None:
    """Finds fragments in molecules."""
    # Load data
    fragments = pd.read_csv(args.fragment_path)
    molecules = pd.read_csv(args.molecule_path)

    # Convert SMILES to RDKit molecules
    molecules['mol'] = [Chem.MolFromSmiles(smiles) for smiles in tqdm(molecules[args.smiles_column])]

    # Find fragments in molecules
    for index, fragment in tqdm(enumerate(fragments)):
        frag_matcher = FragmentMatcher()
        frag_matcher.Init(fragment)
        molecules[f'fragment_{index}'] = [frag_matcher.HasMatch(mol) for mol in tqdm(molecules['mol'], leave=False)]

    # Delete RDKit molecules
    del molecules['mol']

    # Save data
    molecules.to_csv(args.save_path, index=False)


if __name__ == '__main__':
    find_fragments_in_molecules(Args().parse_args())
