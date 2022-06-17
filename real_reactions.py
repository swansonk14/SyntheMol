"""SMARTS representations of the REAL reactions."""
import re
from functools import cache
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


def count_one_reagent(num_r1: int) -> int:
    """Counts the number of feasible molecules created from one reagent.

    :param num_r1: The number of molecules that match the first reagent.
    :return: The number of different molecules that can be constructed using this reaction and the given reagents.
    """
    return num_r1


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
    return re.sub(r'\[([\w|*]+)(:\d+)]', r'[\1]', smarts)


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

    # TODO: cache this (might need to switch to only using strings)
    @cache
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
                 synnet_ids: Optional[set[int]] = None,
                 counting_fn: Optional[Union[Callable[[int, int], int],
                                             Callable[[int, int, int], int]]] = None) -> None:
        """Initializes the Reaction.

        :param reagents: A list of QueryMols containing the reagents of the reaction.
        :param product: A QueryMol containing the product of the reaction.
        :param reaction_id: The ID of the reaction.
        :param real_ids: A set of reaction IDs from the REAL database that this reaction corresponds to.
        :param synnet_ids: A set of IDs from the SynNet database that this reaction corresponds to.
        :param counting_fn: A function that takes in the number of molecules that match each possible reagent
                            and outputs the number of possible product molecules.
        """
        self.reagents = reagents
        self.product = product
        self.id = reaction_id
        self.real_ids = real_ids
        self.synnet_id = synnet_ids
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
    )
]

SYNNET_REACTIONS = [
    Reaction(
        reagents=[
            QueryMol('[c:1](-[C;$(C-c1ccccc1):2](=[OD1:3])-[OH1]):[c:4](-[NH2:5])'),
            QueryMol('[N;!H0;!$(N-N);!$(N-C=N);!$(N(-C=O)-C=O):6]-[C;H1,$(C-[#6]):7]=[OD1]')
        ],
        product=QueryMol('[c:4]2:[c:1]-[C:2](=[O:3])-[N:6]-[C:7]=[N:5]-2'),
        reaction_id=9,
        synnet_ids={8},
        counting_fn=count_two_different_reagents
    ),
    Reaction(
        reagents=[
            QueryMol('[CH0;$(C-[#6]):1]#[NH0:2]')
        ],
        product=QueryMol('[C:1]1=[N:2]-N-N=N-1'),
        reaction_id=10,
        synnet_ids={9},
        counting_fn=count_one_reagent
    ),
    Reaction(
        reagents=[
            QueryMol('[CH0;$(C-[#6]):1]#[NH0:2]'),
            QueryMol('[C;A;!$(C=O):3]-[*;#17,#35,#53]')
        ],
        product=QueryMol('[C:1]1=[N:2]-N(-[C:3])-N=N-1'),
        reaction_id=11,
        synnet_ids={10},
        counting_fn=count_two_different_reagents
    ),
    Reaction(
        reagents=[
            QueryMol('[CH0;$(C-[#6]):1]#[NH0:2]'),
            QueryMol('[C;A;!$(C=O):3]-[*;#17,#35,#53]')
        ],
        product=QueryMol('[C:1]1=[N:2]-N=N-N-1(-[C:3])'),
        reaction_id=12,
        synnet_ids={11},
        counting_fn=count_two_different_reagents
    ),
    Reaction(
        reagents=[
            QueryMol('[CH0;$(C-[#6]):1]#[CH1:2]'),
            QueryMol('[C;H1,H2;A;!$(C=O):3]-[*;#17,#35,#53,OH1]')
        ],
        product=QueryMol('[C:1]1=[C:2]-N(-[C:3])-N=N-1'),
        reaction_id=13,
        synnet_ids={12},
        counting_fn=count_two_different_reagents
    ),
    Reaction(
        reagents=[
            QueryMol('[CH0;$(C-[#6]):1]#[CH1:2]'),
            QueryMol('[C;H1,H2;A;!$(C=O):3]-[*;#17,#35,#53,OH1]')
        ],
        product=QueryMol('[C:1]1=[C:2]-N=NN(-[C:3])-1'),
        reaction_id=14,
        synnet_ids={13},
        counting_fn=count_two_different_reagents
    ),
    Reaction(
        reagents=[
            QueryMol('[CH0;$(C-[#6]):1]#[CH0;$(C-[#6]):2]'),
            QueryMol('[C;H1,H2;A;!$(C=O):3]-[*;#17,#35,#53,OH1]')
        ],
        product=QueryMol('[C:1]1=[C:2]-N=NN(-[C:3])-1'),
        reaction_id=15,
        synnet_ids={14},
        counting_fn=count_two_different_reagents
    ),
    Reaction(
        reagents=[
            QueryMol('[CH0;$(C-[#6]):1]#[NH0:2]'),
            QueryMol('[NH2:3]-[NH1:4]-[CH0;$(C-[#6]);R0:5]=[OD1]')
        ],
        product=QueryMol('[N:2]1-[C:1]=[N:3]-[N:4]-[C:5]=1'),
        reaction_id=16,
        synnet_ids={15},
        counting_fn=count_two_different_reagents
    ),
    Reaction(
        reagents=[
            QueryMol('[CH0;$(C-[#6]):1]#[NH0:2]'),
            QueryMol('[CH0;$(C-[#6]);R0:5](=[OD1])-[#8;H1,$(O-[CH3]),$(O-[CH2]-[CH3])]')
        ],
        product=QueryMol('[N:2]1-[C:1]=N-N-[C:5]=1'),
        reaction_id=17,
        synnet_ids={16},
        counting_fn=count_two_different_reagents
    ),
    Reaction(
        reagents=[
            QueryMol('[c:1](-[C;$(C-c1ccccc1):2](=[OD1:3])-[CH3:4]):[c:5](-[OH1:6])'),
            QueryMol('[C;$(C1-[CH2]-[CH2]-[N,C]-[CH2]-[CH2]-1):7](=[OD1])')
        ],
        product=QueryMol('[O:6]1-[c:5]:[c:1]-[C:2](=[OD1:3])-[C:4]-[C:7]-1'),
        reaction_id=18,
        synnet_ids={17},
        counting_fn=count_two_different_reagents
    ),
    Reaction(
        reagents=[
            QueryMol('[c;r6:1](-[C;$(C=O):6]-[OH1]):[c;r6:2]-[C;H1,$(C-C):3]=[OD1]'),
            QueryMol('[NH2:4]-[NH1;$(N-[#6]);!$(NC=[O,S,N]):5]')
        ],
        product=QueryMol('[c:1]1:[c:2]-[C:3]=[N:4]-[N:5]-[C:6]-1'),
        reaction_id=19,
        synnet_ids={18},
        counting_fn=count_two_different_reagents
    ),
    Reaction(
        reagents=[
            QueryMol('[C;$(C-c1ccccc1):1](=[OD1])-[C;D3;$(C-c1ccccc1):2]~[O;D1,H1]'),
            QueryMol('[CH1;$(C-c):3]=[OD1]')
        ],
        product=QueryMol('[C:1]1-N=[C:3]-[NH1]-[C:2]=1'),
        reaction_id=20,
        synnet_ids={19},
        counting_fn=count_two_different_reagents
    ),
    Reaction(
        reagents=[
            QueryMol('[NH1;$(N-c1ccccc1):1](-[NH2])-[c:5]:[cH1:4]'),
            QueryMol('[C;$(C([#6])[#6]):2](=[OD1])-[CH2;$(C([#6])[#6]);!$(C(C=O)C=O):3]')
        ],
        product=QueryMol('[C:5]1-[N:1]-[C:2]=[C:3]-[C:4]:1'),
        reaction_id=21,
        synnet_ids={20},
        counting_fn=count_two_different_reagents
    ),
    Reaction(
        reagents=[
            QueryMol('[*;Br,I;$(*c1ccccc1)]-[c:1]:[c:2]-[OH1:3]'),
            QueryMol('[CH1:5]#[C;$(C-[#6]):4]')
        ],
        product=QueryMol('[c:1]1:[c:2]-[O:3]-[C:4]=[C:5]-1'),
        reaction_id=22,
        synnet_ids={22},
        counting_fn=count_two_different_reagents
    ),
    Reaction(
        reagents=[
            QueryMol('[*;Br,I;$(*c1ccccc1)]-[c:1]:[c:2]-[SD2:3]-[CH3]'),
            QueryMol('[CH1:5]#[C;$(C-[#6]):4]')
        ],
        product=QueryMol('[c:1]1:[c:2]-[S:3]-[C:4]=[C:5]-1'),
        reaction_id=23,
        synnet_ids={23},
        counting_fn=count_two_different_reagents
    ),
    Reaction(
        reagents=[
            QueryMol('[*;Br,I;$(*c1ccccc1)]-[c:1]:[c:2]-[NH2:3]'),
            QueryMol('[CH1:5]#[C;$(C-[#6]):4]')
        ],
        product=QueryMol('[c:1]1:[c:2]-[N:3]-[C:4]=[C:5]-1'),
        reaction_id=24,
        synnet_ids={24},
        counting_fn=count_two_different_reagents
    ),
    Reaction(
        reagents=[
            QueryMol('[#6:6][C:5]#[#7;D1:4]'),
            QueryMol('[#6:1][C:2](=[OD1:3])[OH1]')
        ],
        product=QueryMol('[#6:6][c:5]1[n:4][o:3][c:2]([#6:1])n1'),
        reaction_id=25,
        synnet_ids={25},
        counting_fn=count_two_different_reagents
    ),
    Reaction(
        reagents=[
            QueryMol('[#6;H0;D3;$([#6](~[#6])~[#6]):1]B(O)O'),
            QueryMol('[#6;H0;D3;$([#6](~[#6])~[#6]):2][Cl,Br,I]')
        ],
        product=QueryMol('[#6:2][#6:1]'),
        reaction_id=26,
        synnet_ids={27},
        counting_fn=count_two_different_reagents
    ),
    Reaction(
        reagents=[
            QueryMol('[C;H1&$(C([#6])[#6]),H2&$(C[#6]):1][OH1]'),
            QueryMol('[NH1;$(N(C=O)C=O):2]')
        ],
        product=QueryMol('[C:1][N:2]'),
        reaction_id=27,
        synnet_ids={29},
        counting_fn=count_two_different_reagents
    ),
    Reaction(
        reagents=[
            QueryMol('[C;H1&$(C([#6])[#6]),H2&$(C[#6]):1][OH1]'),
            QueryMol('[OH1;$(Oc1ccccc1):2]')
        ],
        product=QueryMol('[C:1][O:2]'),
        reaction_id=28,
        synnet_ids={30},
        counting_fn=count_two_different_reagents
    ),
    Reaction(
        reagents=[
            QueryMol('[C;H1&$(C([#6])[#6]),H2&$(C[#6]):1][OH1]'),
            QueryMol('[NH1;$(N([#6])S(=O)=O):2]')
        ],
        product=QueryMol('[C:1][N:2]'),
        reaction_id=29,
        synnet_ids={31},
        counting_fn=count_two_different_reagents
    ),
    Reaction(
        reagents=[
            QueryMol('[#6;$(C=C-[#6]),$(c:c):1][Br,I]'),
            QueryMol('[Cl,Br,I][c:2]')
        ],
        product=QueryMol('[c:2][#6:1]'),
        reaction_id=30,
        synnet_ids={36},
        counting_fn=count_two_different_reagents
    ),
    Reaction(
        reagents=[
            QueryMol('[#6:1][C:2]#[#7;D1]'),
            QueryMol('[Cl,Br,I][#6;$([#6]~[#6]);!$([#6]([Cl,Br,I])[Cl,Br,I]);!$([#6]=O):3]')
        ],
        product=QueryMol('[#6:1][C:2](=O)[#6:3]'),
        reaction_id=31,
        synnet_ids={37},
        counting_fn=count_two_different_reagents
    ),
    Reaction(
        reagents=[
            QueryMol('[#6:1][C;H1,$([C]([#6])[#6]):2]=[OD1:3]'),
            QueryMol('[Cl,Br,I][#6;$([#6]~[#6]);!$([#6]([Cl,Br,I])[Cl,Br,I]);!$([#6]=O):4]')
        ],
        product=QueryMol('[C:1][#6:2]([OH1:3])[#6:4]'),
        reaction_id=32,
        synnet_ids={38},
        counting_fn=count_two_different_reagents
    ),
    Reaction(
        reagents=[
            QueryMol('[Cl,Br,I][c;$(c1:[c,n]:[c,n]:[c,n]:[c,n]:[c,n]:1):1]'),
            QueryMol('[N;$(NC)&!$(N=*)&!$([N-])&!$(N#*)&!$([ND3])&!$([ND4])&!$(N[c,O])&!$(N[C,S]=[S,O,N]),H2&$(Nc1:[c,n]:[c,n]:[c,n]:[c,n]:[c,n]:1):2]')
        ],
        product=QueryMol('[c:1][N:2]'),
        reaction_id=33,
        synnet_ids={42},
        counting_fn=count_two_different_reagents
    ),
    Reaction(
        reagents=[
            QueryMol('[c;$(c1[c;$(c[C,S,N](=[OD1])[*;R0;!OH1])]cccc1):1][C;$(C(=O)[O;H1])]'),
            QueryMol('[c;$(c1aaccc1):2][Cl,Br,I]')
        ],
        product=QueryMol('[c:1][c:2]'),
        reaction_id=34,
        synnet_ids={44},
        counting_fn=count_two_different_reagents
    ),
    Reaction(
        reagents=[
            QueryMol('[Cl:1][CH2:2]-[C$([CH](C)),C$(C(C)(C)):3]=[O:4]'),
            QueryMol('[OH:12]-[c:11]1[c:6][c:7][c:8][c:9][c:10]1-[CH:13]=[O:14]')
        ],
        product=QueryMol('[C:3](=[O:4])-[c:2]1[c:13][c:10]2[c:9][c:8][c:7][c:6][c:11]2[o:12]1'),
        reaction_id=35,
        synnet_ids={54},
        counting_fn=count_two_different_reagents
    ),
    Reaction(
        reagents=[
            QueryMol('[C$(C(#C)([CX4])):2]#[C$(C(#C)([CX4])):1]'),
            QueryMol('[N$(N(~N)([CX4])):5]~[N]~[N]')
        ],
        product=QueryMol('[c:2]1[c:1][n:5][n][n]1'),
        reaction_id=36,
        synnet_ids={58},
        counting_fn=count_two_different_reagents
    ),
    Reaction(
        reagents=[
            QueryMol('[C$(C(=C)([CX4])):2]=[C$(C(=C)([CX4])):1]'),
            QueryMol('[N$(N(~N)([CX4])):5]~[N]~[N]')
        ],
        product=QueryMol('[C:2]1[C:1][N:5][N]=[N]1'),
        reaction_id=37,
        synnet_ids={59},
        counting_fn=count_two_different_reagents
    ),
    Reaction(
        reagents=[
            QueryMol('[C$(C(=C)([CX4,OX2,NX3])([CX4,OX2,NX3])),C$([CH](=C)([CX4,OX2,NX3])),C$([CH2](=C)):1]=[C$(C(=C)([CX4,OX2,NX3])([CX4,OX2,NX3])),C$([CH](=C)([CX4,OX2,NX3])),C$([CH2](=C)):2]'),
            QueryMol('[C$(C(=C)([CX4,OX2,NX3])([CX4,OX2,NX3])),C$([CH](=C)([CX4,OX2,NX3])),C$([CH2](=C)):3]=[C$([C](=C)(C)([CX4,OX2,NX3])),C$([CH](=C)(C)):4]-[C$([C](=C)(C)([CX4,OX2,NX3])),C$([CH](=C)(C)):5]=[C$(C(=C)([CX4,OX2,NX3])([CX4,OX2,NX3])),C$([CH](=C)([CX4,OX2,NX3])),C$([CH2](=C)):6]')
        ],
        product=QueryMol('[C:1]1[C:2][C:3][C:4]=[C:5][C:6]1'),
        reaction_id=38,
        synnet_ids={60},
        counting_fn=count_two_different_reagents
    ),
    Reaction(
        reagents=[
            QueryMol('[C$(C(#C)([CX4,OX2,NX3])),C$([CH](#C)):1]#[C$(C(#C)([CX4,OX2,NX3])),C$([CH](#C)):2]'),
            QueryMol('[C$(C(=C)([CX4,OX2,NX3])([CX4,OX2,NX3])),C$([CH](=C)([CX4,OX2,NX3])),C$([CH2](=C)):3]=[C$([C](=C)(C)([CX4,OX2,NX3])),C$([CH](=C)(C)):4]-[C$([C](=C)(C)([CX4,OX2,NX3])),C$([CH](=C)(C)):5]=[C$(C(=C)([CX4,OX2,NX3])([CX4,OX2,NX3])),C$([CH](=C)([CX4,OX2,NX3])),C$([CH2](=C)):6]')
        ],
        product=QueryMol('[C:1]1=[C:2][C:3][C:4]=[C:5][C:6]1'),
        reaction_id=39,
        synnet_ids={61},
        counting_fn=count_two_different_reagents
    ),
    Reaction(
        reagents=[
            QueryMol('[CH:7](=[O:8])-[c:1]1[c:2][c:3][c:4][c:5][c:6]1'),
            QueryMol('[O:24]=[C:23](-[C:22](=[O:25])-[c:15]1[c:10][c:11][c:12][c:13][c:14]1)-[c:20]1[c:21][c:16][c:17][c:18][c:19]1')
        ],
        product=QueryMol('[nH:27]-1[c:7]([n:26][c:23]([c:22]-1[c:15]1[c:10][c:11][c:12][c:13][c:14]1)-[c:20]1[c:21][c:16][c:17][c:18][c:19]1)-[c:1]1[c:2][c:3][c:4][c:5][c:6]1'),
        reaction_id=40,
        synnet_ids={64},
        counting_fn=count_two_different_reagents
    ),
    Reaction(
        reagents=[
            QueryMol('[O$(O(C)([CX4])):8][C:7](=[O:9])[CH:6][C:5][C:4][C:3][C:2]([O$(O(C)([CX4])):10])=[O:1]')
        ],
        product=QueryMol('[O:8][C:7](=[O:9])[C:6]1[C:5][C:4][C:3][C:2]1=[O:1]'),
        reaction_id=41,
        synnet_ids={66},
        counting_fn=count_one_reagent
    ),
    Reaction(
        reagents=[
            QueryMol('[O$(O(C)([CX4])):8][C:7](=[O:9])[CH:6][C:5][C:11][C:4][C:3][C:2]([O$(O(C)([CX4])):10])=[O:1]')
        ],
        product=QueryMol('[O:8][C:7](=[O:9])[C:6]1[C:5][C:11][C:4][C:3][C:2]1=[O:1]'),
        reaction_id=42,
        synnet_ids={67},
        counting_fn=count_one_reagent
    ),
    Reaction(
        reagents=[
            QueryMol('[Cl:9][C:7](=[O:8])-[c:3]1[c:2][c:1][c:6][c:5][c:4]1'),
            QueryMol('[C$([CH2](C)([CX4])),C$([CH3](C)):18]-[C:16](=[O:17])-[c:14]1[c:13][c:12][c:11][c:10][c:15]1-[OH:19]')
        ],
        product=QueryMol('[O:17]=[C:16]-1-[C:18]=[C:7](-[O:8]-[c:15]2[c:10][c:11][c:12][c:13][c:14]-12)-[c:3]1[c:2][c:1][c:6][c:5][c:4]1'),
        reaction_id=43,
        synnet_ids={68},
        counting_fn=count_two_different_reagents
    ),
    Reaction(
        reagents=[
            QueryMol('[C$(C(C)(=O)([CX4,OX2&H0])),C$(C(C)(#N)),N$([N+1](C)(=O)([O-1])):1][C$([CH]([C,N])([C,N])([CX4])),C$([CH2]([C,N])([C,N])):2][C$(C(C)(=O)([CX4,OX2&H0])),C$(C(C)(#N)),N$([N+1](C)(=O)([O-1])):3]'),
            QueryMol('[C$(C(C)(#N)),C$(C(C)([CX4,OX2&H0])([CX4,OX2&H0])([OX2&H0])),C$([CH](C)([CX4,OX2&H0])([OX2&H0])),C$([CH2](C)([OX2&H0])),C$(C(C)(=O)([OX2&H0])):6][CH:5]=[C$(C(=C)([CX4])([CX4])),C$([CH](=C)([CX4])),C$([CH2](=C)):4]')
        ],
        product=QueryMol('[C:6][C:5][C:4][C:2]([C:1])[C:3]'),
        reaction_id=44,
        synnet_ids={69},
        counting_fn=count_two_different_reagents
    ),
    Reaction(
        reagents=[
            QueryMol('[C$([C](O)([CX4])([CX4])([CX4])),C$([CH](O)([CX4])([CX4])),C$([CH2](O)([CX4])):4]-[O:3]-[C$(C(=O)([CX4])),C$([CH](=O)):2]=[O:5]'),
            QueryMol('[C$([CH](C)([CX4])([CX4])),C$([CH2](C)([CX4])),C$([CH3](C)):7]-[C$(C(=O)([CX4])),C$([CH](=O)):8]=[O:9]')
        ],
        product=QueryMol('[C:7](-[C:2]=[O:5])-[C:8]=[O:9]'),
        reaction_id=45,
        synnet_ids={70},
        counting_fn=count_two_different_reagents
    ),
    Reaction(
        reagents=[
            QueryMol('[Cl,OH,O-:3][C$(C(=O)([CX4,c])),C$([CH](=O)):2]=[O:4]'),
            QueryMol('[O$([OH]([CX4,c])),O$([OH]([CX4,c])([CX4,c])),S$([SH]([CX4,c])),S$([SH]([CX4,c])([CX4,c])):6]')
        ],
        product=QueryMol('[*:6]-[C:2]=[O:4]'),
        reaction_id=46,
        synnet_ids={71},
        counting_fn=count_two_different_reagents
    ),
    Reaction(
        reagents=[
            QueryMol('[C$(C(=O)([CX4,c])([CX4,c])),C$([CH](=O)([CX4,c])):1]=[O:2]'),
            QueryMol('[N$([NH2,NH3+1]([CX4,c])),N$([NH]([CX4,c])([CX4,c])):3]')
        ],
        product=QueryMol('[N+0:3][C:1]'),
        reaction_id=47,
        synnet_ids={72},
        counting_fn=count_two_different_reagents
    ),
    Reaction(
        reagents=[
            QueryMol('[Br:1][c$(c(Br)),n$(n(Br)),o$(o(Br)),C$([CH](Br)(=C)):2]'),
            QueryMol('[C$(C(B)([CX4])([CX4])([CX4])),C$([CH](B)([CX4])([CX4])),C$([CH2](B)([CX4])),C$([CH2](B)),C$(C(B)(=C)),c$(c(B)),o$(o(B)),n$(n(B)):3][B$(B([C,c,n,o])([OH,$(OC)])([OH,$(OC)])),B$([B-1]([C,c,n,o])(N)([OH,$(OC)])([OH,$(OC)])):4]')
        ],
        product=QueryMol('[C,c,n,o:2][C,c,n,o:3]'),
        reaction_id=48,
        synnet_ids={73},
        counting_fn=count_two_different_reagents
    ),
    Reaction(
        reagents=[
            QueryMol('[Br,I:1][C$(C([Br,I])([CX4])([CX4])([CX4])),C$([CH]([Br,I])([CX4])([CX4])),C$([CH2]([Br,I])([CX4])),C$([CH3]([Br,I])),C$([C]([Br,I])(=C)([CX4])),C$([CH]([Br,I])(=C)),C$(C([Br,I])(#C)),c$(c([Br,I])):2]'),
            QueryMol('[Br,I:3][C$(C([Br,I])([CX4])([CX4])([CX4])),C$([CH]([Br,I])([CX4])([CX4])),C$([CH2]([Br,I])([CX4])),C$([CH3]([Br,I])),C$([C]([Br,I])(=C)([CX4])),C$([CH]([Br,I])(=C)),C$(C([Br,I])(#C)),c$(c([Br,I])):4]')
        ],
        product=QueryMol('[C,c:2][C,c:4]'),
        reaction_id=49,
        synnet_ids={74},
        counting_fn=count_two_different_reagents
    ),
    Reaction(
        reagents=[
            QueryMol('[OH,O-]-[C$(C(=O)(O)([CX4,c])):2]=[O:3]'),
            QueryMol('[OH:8]-[C$([CH](O)([CX4,c])([CX4,c])),C$([CH2](O)([CX4,c])),C$([CH3](O)):6]')
        ],
        product=QueryMol('[C:6][O]-[C:2]=[O:3]'),
        reaction_id=50,
        synnet_ids={75},
        counting_fn=count_two_different_reagents
    ),
    Reaction(
        reagents=[
            QueryMol('[C$([CH](=C)([CX4])),C$([CH2](=C)):2]=[C$(C(=C)([CX4])([CX4])),C$([CH](=C)([CX4])),C$([CH2](=C)):3]'),
            QueryMol('[Br,I:7][C$([CX4]([Br,I])),c$([c]([Br,I])):4]')
        ],
        product=QueryMol('[C,c:4][C:2]=[C:3]'),
        reaction_id=51,
        synnet_ids={76},
        counting_fn=count_two_different_reagents
    ),
    Reaction(
        reagents=[
            QueryMol('[Cl,OH,O-:3][C$(C(=O)([CX4,c])),C$([CH](=O)):2]=[O:4]'),
            QueryMol('[N$([NH2,NH3+1]([CX4,c])),N$([NH]([CX4,c])([CX4,c])):6]')
        ],
        product=QueryMol('[N+0:6]-[C:2]=[O:4]'),
        reaction_id=52,
        synnet_ids={77},
        counting_fn=count_two_different_reagents
    ),
    Reaction(
        reagents=[
            QueryMol('[I:1][C$(C(I)([CX4,c])([CX4,c])([CX4,c])),C$([CH](I)([CX4,c])([CX4,c])),C$([CH2](I)([CX4,c])),C$([CH3](I)):2]'),
            QueryMol('[C$(C(=O)([Cl,OH,O-])([CX4,c])),C$([CH]([Cl,OH,O-])(=O)):3](=[O:6])[Cl,OH,O-:5]')
        ],
        product=QueryMol('[C:2]-[C:3]=[O:6]'),
        reaction_id=53,
        synnet_ids={80},
        counting_fn=count_two_different_reagents
    ),
    Reaction(
        reagents=[
            QueryMol('[C$(C(C)([CX4])([CX4])([CX4])),C$([CH](C)([CX4])([CX4])),C$([CH2](C)([CX4])),C$([CH3](C)):1][C:2]#[CH:3]'),
            QueryMol('[Br,I:4][C$(C(=O)([Br,I])([CX4])),C$([CH](=O)([Br,I])):5]=[O:6]')
        ],
        product=QueryMol('[C:1][C:2]#[C:3][C:5]=[O:6]'),
        reaction_id=54,
        synnet_ids={83},
        counting_fn=count_two_different_reagents
    ),
    Reaction(
        reagents=[
            QueryMol('[OH,O-:4]-[C$(C(=O)([OH,O-])([CX4])),C$([CH](=O)([OH,O-])):2]=[O:3]')
        ],
        product=QueryMol('[Cl:5][C:2]=[O:3]'),
        reaction_id=55,
        synnet_ids={84},
        counting_fn=count_one_reagent
    ),
    Reaction(
        reagents=[
            QueryMol('[OH:2]-[$([CX4]),c:1]')
        ],
        product=QueryMol('[Br:3][C,c:1]'),
        reaction_id=56,
        synnet_ids={85},
        counting_fn=count_one_reagent
    ),
    Reaction(
        reagents=[
            QueryMol('[OH:2]-[$([CX4]),c:1]')
        ],
        product=QueryMol('[Cl:3][C,c:1]'),
        reaction_id=57,
        synnet_ids={86},
        counting_fn=count_one_reagent
    ),
    Reaction(
        reagents=[
            QueryMol('[OH,O-:3][S$(S([CX4])):2](=[O:4])=[O:5]')
        ],
        product=QueryMol('[Cl:6][S:2](=[O:5])=[O:4]'),
        reaction_id=58,
        synnet_ids={87},
        counting_fn=count_one_reagent
    ),
    Reaction(
        reagents=[
            QueryMol('[Cl,I,Br:7][c:1]1[c:2][c:3][c:4][c:5][c:6]1')
        ],
        product=QueryMol('[N:9]#[C:8][c:1]1[c:2][c:3][c:4][c:5][c:6]1'),
        reaction_id=59,
        synnet_ids={90},
        counting_fn=count_one_reagent
    )
]

# Check that reaction IDs are unique and count from 1 to the number of reactions
reaction_ids = {reaction.id for reaction in REAL_REACTIONS}

if reaction_ids != set(range(1, len(REAL_REACTIONS) + 1)):
    raise ValueError('Reaction IDs are not unique and/or do not count from 1 to the number of reactions.')
