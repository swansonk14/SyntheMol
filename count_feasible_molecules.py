"""Counts the number of feasible molecules using reaction SMARTS from the top REAL reactions."""
from pathlib import Path

import pandas as pd
from rdkit import Chem
from tap import Tap
from tqdm import tqdm

from real_reactions import REAL_REACTIONS


class Args(Tap):
    data_path: Path  # Path to CSV file containing building blocks.
    smiles_column: str = 'smiles'  # Name of the column in molecule_path containing SMILES.


def count_feasible_reactions(args: Args) -> None:
    """Counts the number of feasible molecules using reaction SMARTS from the top REAL reactions."""
    # Load data
    building_block_data = pd.read_csv(args.data_path)

    # Convert SMILES to mols
    building_block_mols = [
        Chem.AddHs(Chem.MolFromSmiles(smiles))
        for smiles in tqdm(building_block_data[args.smiles_column], desc='SMILES to mol')
    ]

    # Process data
    total_num_molecules = 0
    for index, reaction in enumerate(tqdm(REAL_REACTIONS, desc='Looping reactions')):
        reagents = []
        for reagent_index, reagent in enumerate(tqdm(reaction['reagents'], desc='Looping reagents', leave=False)):
            params = Chem.SubstructMatchParameters()

            if 'checkers' in reaction and reagent_index in reaction['checkers']:
                params.setExtraFinalCheck(reaction['checkers'][reagent_index])

            reagents.append([
                building_block_mol
                for building_block_mol in tqdm(building_block_mols, desc='Matching reagents', leave=False)
                if building_block_mol.HasSubstructMatch(reagent, params)
            ])

        num_molecules = reaction['counting_fn'](*[len(reagent) for reagent in reagents])
        print(f'Reaction {index + 1} can produce {num_molecules:,} molecules')
        total_num_molecules += num_molecules

    print(f'All reactions can produce {total_num_molecules:,} molecules')


if __name__ == '__main__':
    count_feasible_reactions(Args().parse_args())
