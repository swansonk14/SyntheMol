"""SMARTS representations of the REAL reactions."""
import re
from typing import Any, Callable, Optional, Union

from rdkit import Chem
from rdkit.Chem import AllChem
from scipy.special import comb

from collections import deque


Molecule = Union[str, Chem.Mol]  # Either a SMILES string or an RDKit Mol object


class CarbonChainChecker:
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
        return isinstance(other, CarbonChainChecker) and self.smarts == other.smarts


# TODO: document/remove diff
def count_two_different_reagents(num_r1: int, num_r2: int, diff: bool = False) -> int:
    """Counts the number of feasible molecules created from two different reagents.

    :param num_r1: The number of molecules that match the first reagent.
    :param num_r2: The number of molecules that match the second reagent.
    :return: The number of different molecules that can be constructed using this reaction and the given reagents.
    """
    return num_r1 * num_r2


# TODO: document/remove diff
def count_two_same_reagents(num_r1: int, num_r2: int, diff: bool = False) -> int:
    """Counts the number of feasible molecules created from two of the same reagent.

    :param num_r1: The number of molecules that match the first reagent.
    :param num_r2: The number of molecules that match the second reagent (this should be the same as num_r1).
    :return: The number of different molecules that can be constructed using this reaction and the given reagents.
    """
    if diff:
        return num_r1 * num_r2

    assert num_r1 == num_r2

    return comb(num_r1, 2, exact=True, repetition=True)


# TODO: document/remove diff
def count_three_reagents_with_two_same(num_r1: int, num_r2: int, num_r3: int, diff: bool = False) -> int:
    """Counts the number of feasible molecules created from three reagents
       with the last two the same and the first different.

    :param num_r1: The number of molecules that match the first reagent.
    :param num_r2: The number of molecules that match the second reagent.
    :param num_r3: The number of molecules that match the third reagent (this should be the same as num_r2).
    :return: The number of different molecules that can be constructed using this reaction and the given reagents.
    """
    if diff:
        return num_r1 * num_r2 * num_r3

    assert num_r2 == num_r3

    return num_r1 * comb(num_r2, 2, exact=True, repetition=True)


def strip_atom_mapping(smarts: str) -> str:
    """Strips the atom mapping from a SMARTS (i.e., any ":" followed by digits).

    :param smarts: A SMARTS string with atom mapping indices.
    :return: The same SMARTS string but without the atom mapping indices.
    """
    return re.sub(r':\d+', '', smarts)


def convert_to_mol(mol: Molecule, add_hs: bool = False) -> Chem.Mol:
    """Converts a SMILES to an RDKit Mol object (if not already converted) and optionally adds Hs.

    :param mol: A SMILES string or an RDKit Mol object.
    :param add_hs: Whether to add Hs.
    :return: An RDKit Mol object with Hs added optionally.
    """
    if isinstance(mol, str):
        mol = Chem.MolFromSmiles(mol)

    if add_hs:
        mol = Chem.AddHs(mol)

    return mol


class QueryMol:
    """Contains a molecule query in the form of a SMARTS string with helper functions."""

    def __init__(self,
                 smarts: str,
                 checker_class: Optional[type] = None) -> None:
        """Initializes the QueryMol.

        :param smarts: A SMARTS string representing the molecular query.
        :param checker_class: A class that is used as an extra checker for any molecules that match the SMARTS query.
        """
        self.smarts_with_atom_mapping = smarts if ':' in smarts else None
        self.smarts = strip_atom_mapping(smarts)
        self.query_mol = Chem.MolFromSmarts(self.smarts)

        self.params = Chem.SubstructMatchParameters()

        if checker_class is not None:
            self.checker = checker_class(self.smarts)
            self.params.setExtraFinalCheck(self.checker)
        else:
            self.checker = None

    def has_substruct_match(self, mol: Molecule) -> bool:
        """Determines whether the provided molecule includes this QueryMol as a substructure.

        :param mol: A molecule, which can either be a SMILES string or an RDKit Mol object.
        :return: True if the molecule includes this QueryMol as a substructure, False otherwise.
        """
        mol = convert_to_mol(mol, add_hs=True)

        return mol.HasSubstructMatch(self.query_mol, self.params)

    def __hash__(self) -> int:
        """Gets the hash of the QueryMol. Note: The hash depends on the SMARTS *without* atom mapping."""
        return hash(self.smarts)

    def __eq__(self, other: Any) -> bool:
        """Determines equality with another object. Note: The equality depends on the SMARTS *without* atom mapping."""
        return isinstance(other, QueryMol) and self.smarts == other.smarts and self.checker == other.checker


class Reaction:
    """Contains a chemical reaction including SMARTS for the reagents, product, and reaction and helper functions."""

    def __init__(self,
                 reagents: list[QueryMol],
                 product: QueryMol,
                 reaction_id: Optional[int] = None,
                 real_ids: Optional[set[int]] = None,
                 synnet_id: Optional[int] = None,
                 counting_fn: Optional[Union[Callable[[int, int], int],
                                             Callable[[int, int, int], int]]] = None) -> None:
        """Initializes the Reaction.

        :param reagents: A list of QueryMols containing the reagents of the reaction.
        :param product: A QueryMol containing the product of the reaction.
        :param reaction_id: The ID of the reaction.
        :param real_ids: A set of reaction IDs from the REAL database that this reaction corresponds to.
        :param synnet_id: The ID of the corresponding reaction in SynNet.
        :param counting_fn: A function that takes in the number of molecules that match each possible reagent
                            and outputs the number of possible product molecules.
        """
        self.reagents = reagents
        self.product = product
        self.id = reaction_id
        self.real_ids = real_ids
        self.synnet_id = synnet_id
        self.counting_fn = counting_fn

        self.reaction = AllChem.ReactionFromSmarts(
            f'{".".join(f"({reagent.smarts_with_atom_mapping})" for reagent in self.reagents)}'
            f'>>({self.product.smarts_with_atom_mapping})'
        )

    @property
    def num_reagents(self) -> int:
        return len(self.reagents)

    def run_reactants(self, reactants: list[Molecule]) -> tuple[tuple[Chem.Mol, ...], ...]:
        return self.reaction.RunReactants([convert_to_mol(reactant, add_hs=True) for reactant in reactants])

    # TODO: document/remove diff
    def count_feasible_products(self, *num_rs: tuple[int], diff: bool = False) -> int:
        """Counts the number of feasible products of this reaction given the number of each reagent.

        :param num_rs: The number of different molecules that match each reagent in the reaction.
        :return: The number of feasible product molecules this reaction could produce given the number of reagents.
        """
        if self.counting_fn is None:
            raise ValueError('Counting function is not provided for this reaction.')

        return self.counting_fn(*num_rs, diff=diff)


REAL_REACTIONS = [
    Reaction(
        reagents=[
            QueryMol('CC(C)(C)OC(=O)[N:1]([*:2])[*:3].[*:4][N:5]([H])[*:6]', checker_class=CarbonChainChecker),
            QueryMol('[OH1][C:7]([*:8])=[O:9]'),
            QueryMol('[OH1][C:10]([*:11])=[O:12]')
        ],
        product=QueryMol('[*:4][N:5]([*:6])[C:7](=[O:9])[*:8].[*:3][N:1]([*:2])[C:10](=[O:12])[*:11]',
                         checker_class=CarbonChainChecker),
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
    ),
    Reaction(
        reagents=[
            QueryMol('[cH1:1]1:[c:2](-[CH2:7]-[CH2:8]-[NH2:9]):[c:3]:[c:4]:[c:5]:[c:6]:1'),
            QueryMol('[#6:11]-[CH1;R0:10]=[OD1]')
        ],
        product=QueryMol('[c:1]12:[c:2](-[CH2:7]-[CH2:8]-[NH1:9]-[C:10]-2(-[#6:11])):[c:3]:[c:4]:[c:5]:[c:6]:1'),
        reaction_id=9,
        synnet_id=1,
        counting_fn=count_two_different_reagents
    )
]

# Check that reaction IDs are unique and count from 1 to the number of reactions
reaction_ids = {reaction.id for reaction in REAL_REACTIONS}

if reaction_ids != set(range(1, len(REAL_REACTIONS) + 1)):
    raise ValueError('Reaction IDs are not unique and/or do not count from 1 to the number of reactions.')
