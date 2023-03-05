"""QueryMol class for storing a SMARTS query and associated helper functions."""
from functools import cache
from typing import Iterable

from rdkit import Chem

from SyntheMol.utils import convert_to_mol, strip_atom_mapping


class QueryMol:
    """Contains a molecule query in the form of a SMARTS string along with helper functions."""

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
    def allowed_smiles(self) -> list[str] | None:
        """Gets a sorted list of allowed SMILES for this QueryMol."""
        return self._allow_list

    @allowed_smiles.setter
    def allowed_smiles(self, allowed_smiles: Iterable[str]) -> None:
        """Sets the allowed SMILES for this QueryMol. Can only be set once.

        Note: Filters the allowed SMILES to only include those that match the SMARTS of this QueryMol.

        :param allowed_smiles: An iterable of SMILES that are allowed for this QueryMol.
        """
        if self._allow_set is not None:
            raise ValueError('Allowed SMILES has already been set.')

        allowed_smiles = set(allowed_smiles)
        self._allow_set = {smiles for smiles in allowed_smiles if self.has_substruct_match(smiles)}
        self._allow_list = sorted(self._allow_set)
        self._allow_tuple = tuple(self._allow_list)

    def has_substruct_match(self, smiles: str) -> bool:
        """Determines whether the provided molecule includes this QueryMol as a substructure.

        Note: Only checks SMARTS match, not self.allowed_set.

        :param smiles: A SMILES string
        :return: True if the molecule includes this QueryMol as a substructure, False otherwise.
        """
        mol = convert_to_mol(smiles, add_hs=True)

        return mol.HasSubstructMatch(self.query_mol)

    @cache
    def has_match(self, smiles: str) -> bool:
        """Determines whether the provided molecule matches this QueryMol.

        Always checks for SMARTS match and additionally checks for match in self.allowed_set
        if self.allowed_set is not None.

        Note: Caching assumes that self._allow_set does not change.

        :param smiles: A SMILES string.
        :return: True if the molecule matches this QueryMol, False otherwise.
        """
        if self._allow_set is not None and smiles not in self._allow_set:
            return False

        return self.has_substruct_match(smiles)

    def __str__(self) -> str:
        """Gets a string representation of the QueryMol."""
        return f'QueryMol({self.smarts})'
