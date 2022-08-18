"""SMARTS representations of the **generalized** REAL reactions.
Reference: https://docs.google.com/document/d/1jEFqvKT4bFe4Il90Pg53AHCtQQf5eA46RST4ZWQJhyI/edit?usp=sharing
"""
from collections import deque
from typing import Any

from rdkit import Chem

from reactions import (
    count_three_reagents_with_two_same,
    count_two_different_reagents,
    count_two_same_reagents,
    QueryMol,
    Reaction
)


class CarbonPathChecker:
    """Checks whether a SMARTS match with two fragments contains a path
     from one fragment to the other with only non-aromatic carbon atoms."""

    def __init__(self, smarts: str) -> None:
        """Initializes the carbon chain checker.

        Note: This is dependent on RDKit assigning atom indices in the order of the atoms in the SMARTS.

        :param smarts: A SMARTS string containing the query. Must contain precisely two fragments.
        """
        if sum(char == '.' for char in smarts) != 1:
            raise ValueError('SMARTS must contain precisely two fragments (separate by ".").')

        self.smarts = smarts

        first_smarts, second_smarts = smarts.split('.')
        first_mol, second_mol = Chem.MolFromSmarts(first_smarts), Chem.MolFromSmarts(second_smarts)
        first_num_atoms = first_mol.GetNumAtoms()

        # Get the indices of all the "*" atoms, i.e., the beginning of side chains
        self.start_atoms = {
            atom.GetIdx()
            for atom in first_mol.GetAtoms()
            if atom.GetAtomicNum() == 0
        }
        self.end_atoms = {
            atom.GetIdx() + first_num_atoms
            for atom in second_mol.GetAtoms()
            if atom.GetAtomicNum() == 0
        }

    def __call__(self, mol: Chem.Mol, matched_atoms: list[int]) -> bool:
        """Checks whether a molecule that has matched to the SMARTS query satisfies the carbon chain criterion.

        :param mol: The Mol object of the molecule that matches the SMARTS query.
        :param matched_atoms: A list of indices of atoms in the molecule that match the SMARTS query.
                              The ith element of this list is the index of the atom in the molecule
                              that matches the atom at index i in the SMARTS query.
                              (Note: This is technically an RDKit vect object, but it can be treated like a list.)
        :return: Whether the matched molecule satisfies the carbon chain criterion.
        """
        # Initialize a set of visited atoms
        visited_atoms = set(matched_atoms)

        # Get start and end indices in the molecule of the side chains that being with carbon atoms
        start_atom_indices = {
            atom_index for start_atom in self.start_atoms
            if mol.GetAtomWithIdx(atom_index := matched_atoms[start_atom]).GetAtomicNum() == 6
        }
        end_atom_indices = {
            atom_index for end_atom in self.end_atoms
            if mol.GetAtomWithIdx(atom_index := matched_atoms[end_atom]).GetAtomicNum() == 6
        }

        # If none of the side chains of a fragment begin with carbon, return False since there is no path of carbons
        if len(start_atom_indices) == 0 or len(end_atom_indices) == 0:
            return False

        # Loop over the atoms that begin side chains
        for start_atom_index in start_atom_indices:
            # Get the starting atom from its index
            atom = mol.GetAtomWithIdx(start_atom_index)

            # Do a breadth-first search from the start atom to try to find a path with only carbons to an end atom
            queue = deque([atom])
            while queue:
                # Get the next atom in the queue
                atom = queue.pop()

                # Add the atom to the visited set to avoid visiting it again
                visited_atoms.add(atom.GetIdx())

                # Loop through neighboring atoms
                for neighbor_atom in atom.GetNeighbors():
                    # Check if we've reached an end atom and return True if so since we've found a path of carbon atoms
                    if neighbor_atom.GetIdx() in end_atom_indices:
                        return True

                    # Add neighbor atom to the queue if it is not visited and is a non-aromatic carbon
                    if (neighbor_atom.GetIdx() not in visited_atoms
                            and neighbor_atom.GetAtomicNum() == 6
                            and not neighbor_atom.GetIsAromatic()):
                        queue.append(neighbor_atom)

        # If we get here, then there is no path that is all non-aromatic carbon atoms so return False
        return False

    def __hash__(self) -> int:
        return hash(self.smarts)

    def __eq__(self, other: Any) -> bool:
        return isinstance(other, CarbonPathChecker) and self.smarts == other.smarts


REAL_REACTIONS = [
    Reaction(
        reagents=[
            QueryMol('CC(C)(C)OC(=O)[N:1]([*:2])[*:3].[*:4][N:5]([H])[*:6]', checker=CarbonPathChecker),
            QueryMol('[OH1][C:7]([*:8])=[O:9]'),
            QueryMol('[OH1][C:10]([*:11])=[O:12]')
        ],
        product=QueryMol('[*:4][N:5]([*:6])[C:7](=[O:9])[*:8].[*:3][N:1]([*:2])[C:10](=[O:12])[*:11]',
                         checker=CarbonPathChecker),
        reaction_id=1,
        real_ids={275592},
        counting_fn=count_three_reagents_with_two_same
    ),
    Reaction(
        reagents=[
            QueryMol('[*:1][N:2]([H])[*:3]'),
            QueryMol('[F,Cl,Br,I][*:4]')
        ],
        product=QueryMol('[*:1][N:2]([*:3])[*:4]'),
        reaction_id=2,
        real_ids={7, 27, 34, 38, 44, 61, 2230, 269956, 269982, 270122, 270166, 270344, 272692, 272710, 273654},
        counting_fn=count_two_different_reagents
    ),
    Reaction(
        reagents=[
            QueryMol('[*:1][OH1,SH1:2]'),
            QueryMol('[F,Cl,Br,I][*:3]')
        ],
        product=QueryMol('[*:1][O,S:2][*:3]'),
        reaction_id=3,
        real_ids={7, 34, 272692, 272710, 273654},
        counting_fn=count_two_different_reagents
    ),
    Reaction(
        reagents=[
            QueryMol('[*:1][N:2]([H])[*:3]'),
            QueryMol('[*:4][N:5]([H])[*:6]')
        ],
        product=QueryMol('O=C([N:2]([*:1])[*:3])[N:5]([*:4])[*:6]'),
        reaction_id=4,
        real_ids={512, 2430, 2554, 2708, 272164, 273390, 273392, 273452, 273454, 273458, 273464, 273574},
        counting_fn=count_two_same_reagents
    ),
    Reaction(
        reagents=[
            QueryMol('[*:1][N:2]([H])[*:3]'),
            QueryMol('[*:4][N:5]([H])[*:6]')
        ],
        product=QueryMol('O=C(C(=O)[N:2]([*:1])[*:3])[N:5]([*:4])[*:6]'),
        reaction_id=5,
        real_ids={2718, 271948, 271949, 273460, 273462, 273492, 273494, 273496, 273498},
        counting_fn=count_two_same_reagents
    ),
    Reaction(
        reagents=[
            QueryMol('[*:1][N:2]([H])[*:3]'),
            QueryMol('[OH1][C:4]([*:5])=[O:6]')
        ],
        product=QueryMol('[*:5][C:4](=[O:6])[N:2]([*:1])[*:3]'),
        reaction_id=6,
        real_ids={11, 22, 527, 188690, 240690, 269946, 270006, 270112, 272270, 272430, 273450, 273652, 280130},
        counting_fn=count_two_different_reagents
    ),
    Reaction(
        reagents=[
            QueryMol('[O:1]=[C:2]([OH1:3])[*:4]'),
            QueryMol('[F,Cl,Br,I][*:5]')
        ],
        product=QueryMol('[O:1]=[C:2]([*:4])[O:3][*:5]'),
        reaction_id=7,
        real_ids={1458},
        counting_fn=count_two_different_reagents
    ),
    Reaction(
        reagents=[
            QueryMol('[*:1][N:2]([H])[*:3]'),
            QueryMol('[O:4]=[S:5](=[O:6])([F,Cl,Br,I])[*:7]')
        ],
        product=QueryMol('[O:4]=[S:5](=[O:6])([*:7])[N:2]([*:1])[*:3]'),
        reaction_id=8,
        real_ids={20, 40, 196680, 232682, 270084, 270188, 271082, 273578, 274078},
        counting_fn=count_two_different_reagents
    )
]

# Check that reaction IDs are unique and count from 1 to the number of reactions
reaction_ids = {reaction.id for reaction in REAL_REACTIONS}

if reaction_ids != set(range(1, len(REAL_REACTIONS) + 1)):
    raise ValueError('REAL reaction IDs are not unique and/or do not count from 1 to the number of reactions.')
