"""SMARTS representations of the REAL reactions."""
from rdkit import Chem
from rdkit.Chem import AllChem
from scipy.special import comb

from collections import deque


class CarbonChainChecker:
    def __init__(self, start_atoms: set[int], end_atoms: set[int]) -> None:
        self.start_atoms = start_atoms  # Atoms that start side chains in one fragment
        self.end_atoms = end_atoms  # Atoms that end side chains in the other fragment

    def __call__(self, mol: Chem.Mol, matched_atoms: list[int]) -> bool:
        visited_atoms = set(matched_atoms)

        # Get start and end indices on the two side chains which are carbons
        start_atom_indices = {
            atom_index for start_atom in self.start_atoms
            if mol.GetAtomWithIdx(atom_index := matched_atoms[start_atom]).GetAtomicNum() == 6
        }
        end_atom_indices = {
            atom_index for end_atom in self.end_atoms
            if mol.GetAtomWithIdx(atom_index := matched_atoms[end_atom]).GetAtomicNum() == 6
        }

        if len(start_atom_indices) == 0 or len(end_atom_indices) == 0:
            return False

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


def count_two_different_reagents(num_r1: int, num_r2: int) -> int:
    """Counts the number of feasible molecules created from two different reagents."""
    return num_r1 * num_r2


def count_two_same_reagents(num_r1: int, num_r2: int) -> int:
    """Counts the number of feasible molecules created from two of the same reagent."""
    assert num_r1 == num_r2

    return comb(num_r1, 2, exact=True, repetition=True)


def count_three_reagents_with_two_same(num_r1: int, num_r2: int, num_r3: int) -> int:
    """Counts the number of feasible molecules created from three reagents
       with the last two the same and the first different."""
    assert num_r2 == num_r3

    return num_r1 * comb(num_r2, 2, exact=True, repetition=True)


# TODO: requires Chem.AddHs()
REAL_REACTIONS = [
    {
        'reagents': [
            'CC(C)(C)OC(=O)[N:1]([*:2])[*:3].[*:4][N:5]([H])[*:6]',
            # 'CC(C)(C)OC(=O)[N:1]([*:2])[C:3].[C:4][N:5]([H])[*:6]',
            '[OH1][C:7]([*:8])=[O:9]',
            '[OH1][C:10]([*:11])=[O:12]'
        ],
        'product': '[*:4][N:5]([*:6])[C:7](=[O:9])[*:8].[*:3][N:1]([*:2])[C:10](=[O:12])[*:11]',
        # 'product': '[C:4][N:5]([*:6])[C:7](=[O:9])[*:8].[C:3][N:1]([*:2])[C:10](=[O:12])[*:11]',
        'checkers': {  # TODO: maybe check explicitly for carbon chains?
            0: CarbonChainChecker(start_atoms={8, 9}, end_atoms={10, 13})
        },
        'real_ids': {275592},
        'counting_fn': count_three_reagents_with_two_same
    },
    {
        'reagents': [
            '[*:1][OH1,SH1:2]',
            '[F,Cl,Br,I:3][*:4]'
        ],
        'product': '[*:1][O,S:2][*:4]',
        'real_ids': {7, 34, 272692, 272710, 273654},
        'counting_fn': count_two_different_reagents
    },
    {
        'reagents': [
            '[*:1][NH1:2][*:3]',
            '[F,Cl,Br,I:4][*:5]'
        ],
        'product': '[*:1][N:2]([*:2])[*:3]',
        'real_ids': {7, 27, 34, 38, 44, 61, 2230, 269956, 269982, 270122, 270166, 270344, 272692, 272710, 273654},
        'counting_fn': count_two_different_reagents
    },
    {
        'reagents': [
            '[O:1]=[C:2]([OH1:3])[*:4]',
            '[F,Cl,Br,I:5][*:6]'
        ],
        'product': '[O:1]=[C:2]([*:4])[O:3][*:6]',
        'real_ids': {1458},
        'counting_fn': count_two_different_reagents
    },
    {
        'reagents': [
            '[*:1][NH1:2][*:3]',
            '[OH1][C:4]([*:5])=[O:6]'
        ],
        'product': '[*:5][C:4](=[O:6])[N:2]([*:1])[*:3]',
        'real_ids': {11, 22, 527, 188690, 240690, 269946, 270006, 270112, 272270, 272430, 273450, 273652, 280130},
        'counting_fn': count_two_different_reagents
    },
    {
        'reagents': [
            '[*:1][NH1:2][*:3]',
            '[*:4][NH1:5][*:6]'
        ],
        'product': 'O=C([N:2]([*:1])[*:3])[N:5]([*:4])[*:6]',
        'real_ids': {512, 2430, 2554, 2708, 272164, 273390, 273392, 273452, 273454, 273458, 273464, 273574},
        'counting_fn': count_two_same_reagents
    },
    {
        'reagents': [
            '[*:1][NH1:2][*:3]',
            '[*:4][NH1:5][*:6]',
        ],
        'product': 'O=C(C(=O)[N:2]([*:1])[*:3])[N:5]([*:4])[*:6]',
        'real_ids': {2718, 271948, 271949, 273460, 273462, 273492, 273494, 273496, 273498},
        'counting_fn': count_two_same_reagents
    },
    {
        'reagents': [
            '[*:1][NH1:2][*:3]',
            '[O:4]=[S:5](=[O:6])([F,Cl,Br,I:7])[*:8]'
        ],
        'product': '[O:4]=[S:5](=[O:6])([*:8])[N:2]([*:1])[*:3]',
        'real_ids': {20, 40, 196680, 232682, 270084, 270188, 271082, 273578, 274078},
        'counting_fn': count_two_different_reagents
    }
]

# Create reaction SMARTS from reagent and product SMARTS
for reaction in REAL_REACTIONS:
    # TODO: do we need parentheses for the reagents to force them to be separate molecules?
    # reaction = reagent_1.reagent_2...reagent_n>>product
    reaction['reaction'] = AllChem.ReactionFromSmarts(f'{".".join(reaction["reagents"])}>>{reaction["product"]}')

    # Convert strings to SMARTS objects
    reaction['reagents'] = [Chem.MolFromSmarts(smarts) for smarts in reaction['reagents']]
    reaction['product'] = Chem.MolFromSmarts(reaction['product'])
