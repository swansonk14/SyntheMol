"""Contains classes representing chemical reactions, including reagents and products."""
import re
from functools import cache
from typing import Any, Iterable, Optional, Union

from rdkit import Chem
from rdkit.Chem import AllChem


Molecule = Union[str, Chem.Mol]  # Either a SMILES string or an RDKit Mol object


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

    def __init__(self, smarts: str) -> None:
        """Initializes the QueryMol.

        :param smarts: A SMARTS string representing the molecular query.
        """
        self.smarts_with_atom_mapping = smarts
        self.smarts = strip_atom_mapping(smarts)
        self.query_mol = Chem.MolFromSmarts(self.smarts)

        # SMILES that are allowed to match this QueryMol
        self._allow_set = self._allow_list = self._allow_tuple = None

    @property
    def allowed_smiles(self) -> Optional[list[str]]:
        """Gets a sorted list of allowed SMILES for this QueryMol."""
        return self._allow_list

    @allowed_smiles.setter
    def allowed_smiles(self, allowed_smiles: Iterable[str]) -> None:
        """Sets the set/list of allowed SMILES for this QueryMol.

        Note: Filters allowed SMILES based on match to SMARTS.

        :param allowed_smiles: An iterable of SMILES that are allowed for this QueryMol.
        """
        if self._allow_set is not None:
            raise ValueError('Allowed SMILES has already been set.')

        self._allow_set = {smiles for smiles in allowed_smiles if self.has_substruct_match(smiles)}
        self._allow_list = sorted(self._allow_set)
        self._allow_tuple = tuple(self._allow_list)

    # TODO: cache this (might need to switch to only using strings)
    # TODO: replace uses of this with uses of has_match
    @cache
    def has_substruct_match(self, mol: Molecule) -> bool:
        """Determines whether the provided molecule includes this QueryMol as a substructure.

        :param mol: A molecule, which can either be a SMILES string or an RDKit Mol object.
        :return: True if the molecule includes this QueryMol as a substructure, False otherwise.
        """
        mol = convert_to_mol(mol, add_hs=True)

        return mol.HasSubstructMatch(self.query_mol)

    @cache
    def has_match(self, smiles: str) -> bool:
        """Determines whether the provided molecule matches this QueryMol.

        Always checks for SMARTS match and additionally checks for match in self.allowed_set
        if self.allowed_set is not None.

        Note: Caching assumes that self.allow_set does not change.

        :param smiles: A SMILES string.
        :return: True if the molecule matches this QueryMol, False otherwise.
        """
        if self._allow_set is not None and smiles not in self._allow_set:
            return False

        return self.has_substruct_match(smiles)

    # TODO: this doesn't work when using allow_list for has_match since QueryMols with the same SMARTS
    # TODO: in different reactions will have different allow_lists
    def __hash__(self) -> int:
        """Gets the hash of the QueryMol. Note: The hash depends on the SMARTS *without* atom mapping."""
        allow_tuple = self._allow_tuple if self._allow_tuple is not None else tuple()

        return hash((self.smarts, *allow_tuple))

    def __eq__(self, other: Any) -> bool:
        """Determines equality with another object. Note: The equality depends on the SMARTS *without* atom mapping."""
        return isinstance(other, QueryMol) and \
               self.smarts == other.smarts and \
               self._allow_tuple == other._allow_tuple

    def __str__(self) -> str:
        """Gets a string representation of the QueryMol."""
        return self.smarts


class Reaction:
    """Contains a chemical reaction including SMARTS for the reagents, product, and reaction and helper functions."""

    def __init__(self,
                 reagents: list[QueryMol],
                 product: QueryMol,
                 reaction_id: Optional[int] = None) -> None:
        """Initializes the Reaction.

        :param reagents: A list of QueryMols containing the reagents of the reaction.
        :param product: A QueryMol containing the product of the reaction.
        :param reaction_id: The ID of the reaction.
        """
        self.reagents = reagents
        self.product = product
        self.id = reaction_id

        self.reaction = AllChem.ReactionFromSmarts(
            f'{".".join(f"({reagent.smarts_with_atom_mapping})" for reagent in self.reagents)}'
            f'>>({self.product.smarts_with_atom_mapping})'
        )

    @property
    def num_reagents(self) -> int:
        return len(self.reagents)

    # TODO: document
    def run_reactants(self, reactants: list[Molecule]) -> tuple[tuple[Chem.Mol, ...], ...]:
        return self.reaction.RunReactants([convert_to_mol(reactant, add_hs=True) for reactant in reactants])


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
