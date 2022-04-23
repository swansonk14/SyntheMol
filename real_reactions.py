"""SMARTS representations of the REAL reactions."""
"""
from collections import deque

from rdkit import Chem

class CarbonChainChecker:
    def __init__(self, start_atoms: set[int], end_atoms: set[int]) -> None:
        self.start_atoms = start_atoms  # Atoms that start side chains in one fragment
        self.end_atoms = end_atoms  # Atoms that end side chains in the other fragment

    def __call__(self, mol: Chem.Mol, matched_atoms: list[int]) -> bool:
        visited_atoms = set(matched_atoms)
        start_atom_indices = {matched_atoms[start_atom] for start_atom in self.start_atoms}
        end_atom_indices = {matched_atoms[end_atom] for end_atom in self.end_atoms}

        # Loop over the atoms that begin side chains
        for start_atom_index in start_atom_indices:
            atom = mol.GetAtomWithIdx(start_atom_index)

            # Do a breadth-first search from the start atom to try to find a path with only carbons to an end atom
            queue = deque([atom])
            while queue:
                atom = queue.pop()
                visited_atoms.add(atom.GetIdx())

                for neighbor_atom in atom.GetNeighbors():
                    # Check if we've reached an end index
                    if neighbor_atom.GetIdx() in end_atom_indices:
                        return True

                    # Add neighbor atom if it is not visited and is a non-aromatic carbon
                    if (neighbor_atom.GetIdx() not in visited_atoms
                            and neighbor_atom.GetAtomicNum() == 6
                            and not neighbor_atom.GetIsAromatic()):
                        queue.append(neighbor_atom)

        return False  # There is no path from any start index to any end index that is all carbon
"""

# TODO: requires Chem.AddHs()
REAL_REACTIONS = [
    {
        'reagents': [
            # reagent_with_carbon_chain := 'CC(C)(C)OC(=O)[N:1]([*:2])[*:3].[*:4][N:5]([H])[*:6]',
            'CC(C)(C)OC(=O)[N:1]([*:2])[C:3].[C:4][N:5]([H])[*:6]',
            '[OH1][C:7]([*:8])=[O:9]',
            '[OH1][C:10]([*:11])=[O:12]'
        ],
        # 'product': '[*:4][N:5]([*:6])[C:7](=[O:9])[*:8].[*:3][N:1]([*:2])[C:10](=[O:12])[*:11]',
        'product': '[C:4][N:5]([*:6])[C:7](=[O:9])[*:8].[C:3][N:1]([*:2])[C:10](=[O:12])[*:11]',
        # 'checkers': {  # TODO: maybe check explicitly for carbon chains?
        #     reagent_with_carbon_chain: CarbonChainChecker(start_atoms={8, 9}, end_atoms={10, 13})
        # },
        'real_ids': {275592}
    },
    {
        'reagents': [
            '[*:1][NH1:2][*:3]',
            '[OH1][C:4]([*:5])=[O:6]'
        ],
        'product': '[*:5][C:4](=[O:6])[N:2]([*:1])[*:3]',
        'real_ids': {11, 22, 527, 188690, 240690, 269946, 270006, 270112, 272270, 272430, 273450, 273652, 280130}
    },
    {
        'reagents': [
            '[*:1][NH1:2][*:3]',
            '[*:4][NH1:5][*:6]'
        ],
        'product': 'O=C([N:2]([*:1])[*:3])[N:5]([*:4])[*:6]',
        'real_ids': {512, 2430, 2554, 2708, 272164, 273390, 273392, 273452, 273454, 273458, 273464, 273574}
    },
    {
        'reagents': [
            '[*:1][NH1:2][*:3]',
            '[F,Cl,Br,I:4][*:5]'
        ],
        'product': '[*:1][N:2]([*:2])[*:3]',
        'real_ids': {7, 27, 34, 38, 44, 61, 2230, 269956, 269982, 270122, 270166, 270344, 272692, 272710, 273654}
    },
    {
        'reagents': [
            '[*:1][O,S:2]',
            '[F,Cl,Br,I:3][*:4]'
        ],
        'product': '[*:1][O,S:2][*:4]',
        'real_ids': {7, 34, 272692, 272710, 273654}
    },
    {
        'reagents': [
            '[*:1][NH1:2][*:3]',
            '[*:4][NH1:5][*:6]',
        ],
        'product': 'O=C(C(=O)[N:2]([*:1])[*:3])[N:5]([*:4])[*:6]',
        'real_ids': {2718, 271948, 271949, 273460, 273462, 273492, 273494, 273496, 273498}
    },
    {
        'reagents': [
            '[*:1][NH1:2][*:3]',
            '[O:4]=[S:5](=[O:6])([F,Cl,Br,I:7])[*:8]'
        ],
        'product': '[O:4]=[S:5](=[O:6])([*:8])[N:2]([*:1])[*:3]',
        'real_ids': {20, 40, 196680, 232682, 270084, 270188, 271082, 273578, 274078}
    },
    {
        'reagents': [
            '[O:1]=[C:2]([O:3])[*:4]'
            '[F,Cl,Br,I:5][*:6]'
        ],
        'product': '[O:1]=[C:2]([*:4])[O:3][*:6]',
        'real_ids': {1458}
    }
]

# Create reaction SMARTS from reagent and product SMARTS
for reaction in REAL_REACTIONS:
    # TODO: do we need parentheses for the reagents to force them to be separate molecules?
    # reaction = reagent_1.reagent_2...reagent_n>>product
    reaction['reaction'] = f'{".".join(reaction["reagent"])}>>{reaction["product"]}'
