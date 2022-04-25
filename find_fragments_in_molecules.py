"""Finds fragments in molecules and plots a histogram of number of fragments that match to each molecule."""
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from rdkit import Chem
from rdkit.Chem.FragmentMatcher import FragmentMatcher
from tap import Tap
from tqdm import tqdm


class Args(Tap):
    fragment_path: Path  # Path to CSV file containing molecular fragments.
    fragment_smiles_path: str = 'smiles'  # Name of the column in fragment_path containing SMILES.
    molecule_path: Path  # Path to CSV file containing full molecules.
    molecule_smiles_column: str = 'smiles'  # Name of the column in molecule_path containing SMILES.
    save_path: Path  # Path to PDF file where fragment count histogram will be saved.

    def process_args(self) -> None:
        self.save_path.parent.mkdir(parents=True, exist_ok=True)


def find_fragments_in_molecules(args: Args) -> None:
    """Finds fragments in molecules and plots a histogram of number of fragments that match to each molecule."""
    # Load fragments
    fragments = pd.read_csv(args.fragment_path)
    print(f'Number of fragments = {len(fragments):,}')

    fragments = set(fragments[args.fragment_smiles_path])
    print(f'Number of unique fragments = {len(fragments):,}')

    # Map from fragment SMILES to fragment matcher
    frag_matchers = {}
    for fragment in tqdm(fragments, desc='Initialize fragments'):
        frag_matcher = FragmentMatcher()
        frag_matcher.Init(fragment)
        frag_matchers[fragment] = frag_matcher

    # Load molecules
    molecules = pd.read_csv(args.molecule_path)
    print(f'Number of molecules = {len(molecules):,}')

    # Convert SMILES to RDKit molecules
    mols = [Chem.MolFromSmiles(smiles) for smiles in tqdm(molecules[args.molecule_smiles_column], desc='SMILES to mol')]

    # Filter out bad molecules
    mols = [mol for mol in mols if mol is not None]
    print(f'Number of valid molecules = {len(molecules):,}')

    # Count fragments in molecules
    num_frags_per_molecule = [
        sum(frag_matchers[fragment].HasMatch(mol) for fragment in fragments)
        for mol in tqdm(mols, desc='Matching fragments')
    ]

    # Histogram of fragment counts per molecule
    plt.hist(num_frags_per_molecule, bins=100)
    plt.xlabel('Number of fragments per molecule')
    plt.ylabel('Count')
    plt.title('Fragment counts per molecule')
    plt.savefig(args.save_path)


if __name__ == '__main__':
    find_fragments_in_molecules(Args().parse_args())
