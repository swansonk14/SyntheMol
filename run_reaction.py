"""Runs a chemical reaction (ID 22) using building blocks."""
from copy import deepcopy
from itertools import product
from multiprocessing import Pool
from pathlib import Path

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from tap import Tap
from tqdm import tqdm


# Reaction ID 22
# TODO: how to avoid matching R1 and R2 that are part of the same connected component (e.g., ring) or is it okay as is?
REAGENT_1 = '[*:1][NH1:2][*:3]'
PATTERN_1 = Chem.MolFromSmarts(REAGENT_1)
REAGENT_2 = '[OH1][C:4]([*:5])=[O:6]'
PATTERN_2 = Chem.MolFromSmarts(REAGENT_2)
PRODUCT = '[*:5][C:4](=[O:6])[N:2]([*:1])[*:3]'
REACTION_SMARTS = f'{REAGENT_1}.{REAGENT_2}>>{PRODUCT}'
REACTION = AllChem.ReactionFromSmarts(REACTION_SMARTS)


class Args(Tap):
    fragment_path: Path  # Path to CSV file containing molecular fragments.
    smiles_column: str = 'smiles'  # Path to column in fragment_path containing SMILES.
    save_path: Path  # Path to CSV file where the products of the reaction will be saved.

    def process_args(self) -> None:
        self.save_path.parent.mkdir(parents=True, exist_ok=True)


def run(mol_pair):
    mol_pair = (Chem.MolFromSmiles(mol_pair[0]), Chem.MolFromSmiles(mol_pair[1]))
    reaction = AllChem.ReactionFromSmarts(REACTION_SMARTS)
    return reaction.RunReactants(mol_pair)


def run_reaction(args: Args) -> None:
    """Runs a chemical reaction (ID 22) using building blocks."""
    # Load fragments
    fragments = pd.read_csv(args.fragment_path)

    # Convert SMILES to mols
    smiles = list(fragments[args.smiles_column])
    mols = [Chem.MolFromSmiles(smiles) for smiles in tqdm(fragments[args.smiles_column], desc='Mol to SMILES')]

    # Match reactants
    # mols_1 = [mol for mol in tqdm(mols, desc='Matching reagent 1') if mol.HasSubstructMatch(PATTERN_1)]
    mols_1 = [smile for mol, smile in tqdm(zip(mols, smiles), desc='Matching reagent 1') if mol.HasSubstructMatch(PATTERN_1)]
    print(f'Number of first reactant = {len(mols_1):,}')

    # mols_2 = [mol for mol in tqdm(mols, desc='Matching reagent 2') if mol.HasSubstructMatch(PATTERN_2)]
    mols_2 = [smile for mol, smile in tqdm(zip(mols, smiles), desc='Matching reagent 2') if mol.HasSubstructMatch(PATTERN_2)]
    print(f'Number of second reactant = {len(mols_2):,}')

    # Determine number of combinations of the two reactants
    num_combinations = len(mols_1) * len(mols_2)
    print(f'Number of combinations = {num_combinations:,}')

    pairs = list(product(mols_1[:10], mols_2[:100000]))
    # pairs = [(deepcopy(mol_1), deepcopy(mol_2)) for mol_1, mol_2 in tqdm(pairs, desc='copy')]
    # Run reaction with all combinations
    with Pool() as pool:
        for products in tqdm(map(run, pairs), total=len(pairs)):
            if any(len(p) != 1 for p in products):
                breakpoint()
    #
    # for mol_pair in tqdm(product(mols_1, mols_2), total=num_combinations, desc='Running reaction'):
    #     # print(Chem.MolToSmiles(mol_pair[0]))
    #     # print(Chem.MolToSmiles(mol_pair[1]))
    #     products = REACTION.RunReactants(mol_pair)
    #     # if len(products) != 2:
    #     #     breakpoint()
    #
    #     if any(len(p) != 1 for p in products):
    #         breakpoint()


        # print(len(products))
        # print(products)
        # print(Chem.MolToSmiles(products[0][0]))
        # print(Chem.MolToSmiles(products[1][0]))
        # breakpoint()


    # from rdkit.Chem.Draw import rdMolDraw2D
    # patt = Chem.MolFromSmarts(REAGENT_2)
    # for i, mol in enumerate(tqdm(mols)):
    #     match = mol.GetSubstructMatch(patt)
    #     if len(match) > 0:
    #         d = rdMolDraw2D.MolDraw2DSVG(500, 500)
    #         rdMolDraw2D.PrepareAndDrawMolecule(d, mol, highlightAtoms=match)
    #         d.FinishDrawing()
    #         with open(f'images/{i}.svg', 'w') as f:
    #             f.write(d.GetDrawingText())




if __name__ == '__main__':
    run_reaction(Args().parse_args())
