"""Contains classes representing chemical reactions, including reagents and products."""
import re
from collections import deque
from functools import cache
from typing import Any, Callable, Iterable, Literal, Optional, Union

from rdkit import Chem
from rdkit.Chem import AllChem
from scipy.special import comb


Molecule = Union[str, Chem.Mol]  # Either a SMILES string or an RDKit Mol object


# TODO: does this need to be a path checker instead of a chain checker?
class CarbonPathChecker:
    """Checks whether a SMARTS match with two fragments contains a path
     from one fragment to the other with only non-aromatic C-C connections."""

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
            if '*' in atom.GetSmarts()
        }
        self.end_atoms = {
            atom.GetIdx() + first_num_atoms
            for atom in second_mol.GetAtoms()
            if '*' in atom.GetSmarts()
        }

    @staticmethod
    def is_ch2(atom: Chem.Atom) -> bool:
        if atom.GetAtomicNum() != 6 or atom.GetIsAromatic():
            return False

        neighbors = atom.GetNeighbors()

        if len(neighbors) != 4:
            return False

        num_h_neighbors = sum(neighbor.GetAtomicNum() == 1 for neighbor in neighbors)

        return num_h_neighbors == 2

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
            if self.is_ch2(mol.GetAtomWithIdx(atom_index := matched_atoms[start_atom]))
        }
        end_atom_indices = {
            atom_index for end_atom in self.end_atoms
            if self.is_ch2(mol.GetAtomWithIdx(atom_index := matched_atoms[end_atom]))
        }

        # If none of the side chains of a fragment begin with carbon, return False since there is no path of carbons
        if len(start_atom_indices) == 0 or len(end_atom_indices) == 0:
            return False

        # Loop over the atoms that begin side chains
        for start_atom_index in start_atom_indices:
            # Get the starting atom from its index
            atom = mol.GetAtomWithIdx(start_atom_index)

            # Do a breadth-first search from the start atom to try to find a path with only CH2 to an end atom
            queue = deque([atom])
            while queue:
                # Get the next atom in the queue
                atom = queue.pop()

                # Add the atom to the visited set to avoid visiting it again
                visited_atoms.add(atom.GetIdx())

                # Loop through neighboring atoms
                for neighbor in atom.GetNeighbors():
                    # Check if we've reached an end atom and return True if so since we've found a path of carbon atoms
                    if neighbor.GetIdx() in end_atom_indices:
                        return True

                    # Add neighbor atom to the queue if it is not visited and is a CH2
                    if neighbor.GetIdx() not in visited_atoms and self.is_ch2(neighbor):
                        queue.append(neighbor)

        # If we get here, then there is no path that satisfies the carbon chain criterion so return False
        return False

    def __hash__(self) -> int:
        return hash(self.smarts)

    def __eq__(self, other: Any) -> bool:
        return isinstance(other, CarbonPathChecker) and self.smarts == other.smarts

    def __str__(self) -> str:
        return self.__class__.__name__


# TODO: should this restrict just to aromatic rings and not allow other atoms?
def aryl_checker(mol: Chem.Mol, matched_atoms: list[int], r_group_start_atom: Chem.Atom) -> bool:
    """Checks whether an R group of a molecule is an aryl group."""
    # Do a breadth-first search from the R group start atom and check that the group is aryl
    queue = deque([r_group_start_atom])
    visited_atoms = set(matched_atoms)
    has_aromatic_ring = False

    while queue:
        # Get the next atom in the queue
        atom = queue.pop()

        # Check whether the atom is aromatic and if so, whether it is carbon
        if atom.GetIsAromatic():
            has_aromatic_ring = True

            if atom.GetAtomicNum() != 6:
                return False

        # Add atom to visited atoms
        visited_atoms.add(atom.GetIdx())

        # Add unvisited neighbors to the queue
        queue += [neighbor for neighbor in atom.GetNeighbors() if neighbor.GetIdx() not in visited_atoms]

    return has_aromatic_ring


def alkyl_checker(mol: Chem.Mol, matched_atoms: list[int], r_group_start_atom: Chem.Atom) -> bool:
    """Checks whether an R group of a molecule is an alkyl group."""
    # Do a search from the R group start atom and check that the group is alkyl
    queue = deque([r_group_start_atom])
    visited_atoms = set(matched_atoms)

    while queue:
        # Get the next atom in the queue
        atom = queue.pop()

        # Check that atom satisfies alkyl property (must be carbon with all single bonds so 4 neighbors)
        neighbors = atom.GetNeighbors()

        if atom.GetAtomicNum() != 6 or atom.GetIsAromatic() or len(neighbors) != 4:
            return False

        # Add the atom to the visited set
        visited_atoms.add(atom.GetIdx())

        # Add neighboring atoms (skipping hydrogens) to the queue if not already visited
        queue += [
            neighbor
            for neighbor in neighbors
            if neighbor.GetIdx() not in visited_atoms and neighbor.GetAtomicNum() != 1
        ]

    # If we get here, then all atoms in the R group satisfy the alkyl property
    return True


def h_checker(mol: Chem.Mol, matched_atoms: list[int], r_group_start_atom: Chem.Atom) -> bool:
    """Checks whether an R group of a molecule is a hydrogen atom."""
    return r_group_start_atom.GetAtomicNum() == 1


def cycle_checker(mol: Chem.Mol, matched_atoms: list[int], r_group_start_atom: Chem.Atom) -> bool:
    """Checks whether an R group of a molecule is a cycle."""
    # Do a search from the R group start atom to get all the atoms in the R group
    visited_atoms = set(matched_atoms)
    queue = deque([r_group_start_atom])

    while queue:
        # Get the next atom in the queue
        atom = queue.pop()

        # Check if atom is in ring
        if not atom.IsInRing():
            return False

        # Add the atom to the visited set
        visited_atoms.add(atom.GetIdx())

        # Add the unvisited neighbors that are not hydrogen
        queue += [
            neighbor
            for neighbor in atom.GetNeighbors()
            if neighbor.GetIdx() not in visited_atoms and neighbor.GetAtomicNum() != 1
        ]

    # Get R group atoms
    r_group_atoms = (visited_atoms - set(matched_atoms)) | {r_group_start_atom.GetIdx()}

    # Get molecule ring info
    rings = mol.GetRingInfo()

    # Ensure that at least one ring contains all the atoms in the R group
    for atom_ring in rings.AtomRings():
        if r_group_atoms <= set(atom_ring):
            return True

    return False


class RGroupChecker:
    """Checks whether each R group in a molecule satisfies at least one checker."""

    def __init__(self, smarts: str, checkers: set[Callable[[Chem.Mol, list[int], Chem.Atom], bool]]) -> None:
        self.smarts = smarts
        self.checkers = checkers

        mol = Chem.MolFromSmarts(self.smarts)

        # Get the indices of all the "*" atoms, i.e., the beginning of side chains
        self.r_group_start_atom_smarts_indices = {
            atom.GetIdx()
            for atom in mol.GetAtoms()
            if '*' in atom.GetSmarts()
        }

    def __call__(self, mol: Chem.Mol, matched_atoms: list[int]) -> bool:
        """Checks whether a molecule that has matched to the SMARTS has R groups that satisfy the checkers.

        Each R group must match at least one checker.

        :param mol: The Mol object of the molecule that matches the SMARTS query.
        :param matched_atoms: A list of indices of atoms in the molecule that match the SMARTS query.
                              The ith element of this list is the index of the atom in the molecule
                              that matches the atom at index i in the SMARTS query.
                              (Note: This is technically an RDKit vect object, but it can be treated like a list.)
        :return: Whether each R group in the matched molecule satisfies at least one checker.
        """
        r_group_start_atom_mol_indices = {
            matched_atoms[r_group_start_atom_smarts_index]
            for r_group_start_atom_smarts_index in self.r_group_start_atom_smarts_indices
        }

        return all(
            any(
                checker(mol, matched_atoms, mol.GetAtomWithIdx(r_group_start_atom_mol_index))
                for checker in self.checkers
            ) for r_group_start_atom_mol_index in r_group_start_atom_mol_indices
        )

    def __hash__(self) -> int:
        return hash(self.smarts) + sum(hash(checker) for checker in self.checkers)

    def __eq__(self, other: Any) -> bool:
        return isinstance(other, RGroupChecker) and self.smarts == other.smarts and self.checkers == other.checkers

    def __str__(self) -> str:
        return f'{self.__class__.__name__}({",".join(sorted(checker.__name__ for checker in self.checkers))})'


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
    return re.sub(r'\[([^:]+)(:\d+)]', r'[\1]', smarts)


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
                 checker_class: Optional[type] = None,
                 checker_kwargs: Optional[dict[str, Any]] = None) -> None:
        """Initializes the QueryMol.

        :param smarts: A SMARTS string representing the molecular query.
        :param checker_class: An extra checker class for any molecules that match the SMARTS query.
        :param checker_kwargs: Additional keyword arguments (outside of smarts) for the checker.
        """
        self.smarts_with_atom_mapping = smarts
        self.smarts = strip_atom_mapping(smarts)
        self.query_mol = Chem.MolFromSmarts(self.smarts)

        self.params = Chem.SubstructMatchParameters()

        if checker_class is not None:
            if checker_kwargs is None:
                checker_kwargs = dict()

            self.checker = checker_class(smarts=self.smarts, **checker_kwargs)
            self.params.setExtraFinalCheck(self.checker)
        else:
            self.checker = None

        # Set/list of SMILES that are allowed to match this QueryMol
        self._allow_set = None
        self._allow_list = None

    @property
    def allowed_smiles(self) -> list[str]:
        """Gets a sorted list of allowed SMILES for this QueryMol."""
        return self._allow_list

    @allowed_smiles.setter
    def allowed_smiles(self, allowed_smiles: Iterable[str]) -> None:
        """Sets the set/list of allowed SMILES for this QueryMol.

        :param allowed_smiles: An iterable of SMILES that are allowed for this QueryMol.
        """
        if self._allow_set is not None:
            raise ValueError('Allowed SMILES has already been set.')

        self._allow_set = set(allowed_smiles)
        self._allow_list = sorted(self._allow_set)

    # TODO: cache this (might need to switch to only using strings)
    @cache
    def has_substruct_match(self, mol: Molecule) -> bool:
        """Determines whether the provided molecule includes this QueryMol as a substructure.

        :param mol: A molecule, which can either be a SMILES string or an RDKit Mol object.
        :return: True if the molecule includes this QueryMol as a substructure, False otherwise.
        """
        mol = convert_to_mol(mol, add_hs=True)

        return mol.HasSubstructMatch(self.query_mol, self.params)

    @cache
    def has_match(self, smiles: str, mode: Literal['substructure', 'allow_set']) -> bool:
        """Determines whether the provided molecule matches this QueryMol.

        Note: Caching assumes that self.allow_set does not change.

        :param smiles: A SMILES string.
        :param mode: Method of determining whether the molecule matches this QueryMol.
                     substructure: Checks whether this QueryMol is a substructure of the molecule using SMARTS matching.
                     allow_set: Checks whether this QueryMol's allow_list includes the SMILES of the molecule.
                                Note: Uses exact string matching with no SMILES canonicalization.
        :return: True if the molecule matches this QueryMol, False otherwise.
        """
        match mode:
            case 'substructure':
                return self.has_substruct_match(smiles)
            case 'allow_set':
                if self._allow_set is None:
                    raise ValueError('Cannot use "allow_set" mode when self.allow_set is None.')

                return smiles in self._allow_set
            case _:
                raise ValueError(f'Mode "{mode}" is not supported.')

    # TODO: this doesn't work when using allow_list for has_match since QueryMols with the same SMARTS
    # TODO: in different reactions will have different allow_lists
    def __hash__(self) -> int:
        """Gets the hash of the QueryMol. Note: The hash depends on the SMARTS *without* atom mapping."""
        return hash(self.smarts)

    def __eq__(self, other: Any) -> bool:
        """Determines equality with another object. Note: The equality depends on the SMARTS *without* atom mapping."""
        return isinstance(other, QueryMol) and self.smarts == other.smarts and self.checker == other.checker

    def __str__(self) -> str:
        """Gets a string representation of the QueryMol."""
        if self.checker is None:
            return self.smarts

        return f'{self.smarts}_{self.checker}'


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


def set_allowed_reaction_smiles(reactions: list[Reaction],
                                reaction_to_reagents: dict[int, dict[int, set[str]]]) -> None:
    """Sets the allowed SMILES for each reaction/reagent in a list of Reactions.

    Note: Modifies Reactions in place.

    :param reactions: A list of Reactions whose allowed SMILES will be set.
    :param reaction_to_reagents: A dictionary mapping from reaction ID to reagent index to a set of allowed SMILES.
    """
    for reaction in reactions:
        for reagent_index, reagent in enumerate(reaction.reagents):
            reagent.allowed_smiles = reaction_to_reagents[reaction.id][reagent_index]
