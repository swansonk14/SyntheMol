"""Counts the number of feasible molecules using reaction SMARTS from the top REAL reactions."""
from pathlib import Path

import pandas as pd
from tap import Tap
from tqdm import tqdm

from reactions import convert_to_mol
from real_reactions import REAL_REACTIONS


class Args(Tap):
    fragment_path: Path  # Path to CSV file containing building blocks.
    smiles_column: str = 'smiles'  # Name of the column in molecule_path containing SMILES.


def count_feasible_molecules_one_reaction(args: Args) -> None:
    """Counts the number of feasible molecules using reaction SMARTS from the top REAL reactions."""
    # Load data
    fragment_data = pd.read_csv(args.fragment_path)

    # Remove duplicate fragments
    fragments = set(fragment_data[args.smiles_column])

    # Convert SMILES to mols
    fragment_mols = [
        convert_to_mol(smiles, add_hs=True)
        for smiles in tqdm(fragments, desc='SMILES to mol')
    ]

    # Count feasible molecules for each reaction
    total_num_molecules = 0
    for reaction in tqdm(REAL_REACTIONS, desc='Looping reactions'):
        num_molecules_per_reagent = [
            sum(
                reagent.has_substruct_match(fragment_mol)
                for fragment_mol in tqdm(fragment_mols, desc='Matching fragments', leave=False)
            )
            for reagent in tqdm(reaction.reagents, desc='Looping reagents', leave=False)
        ]

        num_molecules = reaction.count_feasible_products(*num_molecules_per_reagent)
        print(f'Reaction {reaction.id} can produce {num_molecules:,} molecules')
        total_num_molecules += num_molecules

    print(f'All reactions can produce {total_num_molecules:,} molecules')


if __name__ == '__main__':
    count_feasible_molecules_one_reaction(Args().parse_args())
