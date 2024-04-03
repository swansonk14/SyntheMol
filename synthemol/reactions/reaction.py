"""Reaction class and helper functions."""
from typing import Any

from rdkit import Chem
from rdkit.Chem import AllChem

from synthemol.reactions.query_mol import QueryMol
from synthemol.utils import convert_to_mol, MOLECULE_TYPE


class Reaction:
    """A chemical reaction including SMARTS for the reactants, product, and reaction along with helper functions."""

    def __init__(
        self,
        reactants: list[QueryMol],
        product: QueryMol,
        chemical_space: str | None = None,
        reaction_id: int | None = None,
        sub_reaction_id: int | None = None,
        post_reaction: "Reaction" | None = None
    ) -> None:
        """Initializes the Reaction.

        :param reactants: A list of QueryMols containing the reactants of the reaction.
        :param product: A QueryMol containing the product of the reaction.
        :param chemical_space: The chemical space of the reaction (e.g., Enamine or WuXi).
        :param reaction_id: The ID of the reaction.
        :param sub_reaction_id: The ID of the sub-reaction.
        :param post_reaction: An optional Reaction to run after the main reaction on the product (e.g., BOC cleavage).
        """
        self.reactants = reactants
        self.product = product
        self.chemical_space = chemical_space
        self.reaction_id = reaction_id
        self.sub_reaction_id = sub_reaction_id
        self.id = (
            str(reaction_id)
            if sub_reaction_id is None
            else f"{reaction_id}_{sub_reaction_id}"
        )
        self.post_reaction = post_reaction

        self.reaction_smarts = (
            f'{".".join(f"({reactant.smarts_with_atom_mapping})" for reactant in self.reactants)}'
            f">>({self.product.smarts_with_atom_mapping})"
        )
        self.reaction = AllChem.ReactionFromSmarts(self.reaction_smarts)

    def get_reactant_matches(self, smiles: str) -> list[int]:
        """Gets the indices of the reactants that match the provided SMILES.

        :param smiles: The SMILES to match.
        :return: A list of indices of the reactants that match the provided SMILES.
        """
        return [
            reactant_index
            for reactant_index, reactant in enumerate(self.reactants)
            if reactant.has_match(smiles)
        ]

    @property
    def num_reactants(self) -> int:
        """Gets the number of reactants in the reaction."""
        return len(self.reactants)

    def run_reactants(self, reactants: list[MOLECULE_TYPE]) -> list[str]:
        """Runs the reaction on the provided reactants.

        :param reactants: A list of reactants.
        :return: A list of product SMILES.
        """
        # Convert reactants to mols (and add Hs)
        reactant_mols = [
            convert_to_mol(reactant, add_hs=True) for reactant in reactants
        ]

        # Run reaction on reactants
        product_mols = self.reaction.RunReactants(reactant_mols)

        # Ensure each product has one molecule
        assert all(len(product_mol) == 1 for product_mol in product_mols)

        # Convert product mols to SMILES (and remove Hs)
        products = [
            Chem.MolToSmiles(Chem.RemoveHs(product_mol[0]))
            for product_mol in product_mols
        ]

        # Optionally run post-reaction
        if self.post_reaction is not None:
            products = [self.post_reaction.run_reactants([product]) for product in products]

        return products

    def __str__(self) -> str:
        """Gets the string representation of the Reaction."""
        return f"{self.__class__.__name__}(space={self.chemical_space}, id={self.id}, reaction={self.reaction_smarts})"

    def __repr__(self) -> str:
        """Gets the representation of the Reaction."""
        return str(self)

    def __hash__(self) -> int:
        """Gets the hash of the Reaction."""
        return hash((self.chemical_space, self.id, self.reaction_smarts))

    def __eq__(self, other: Any) -> bool:
        """Determines whether the Reaction is equal to another object."""
        return (
            isinstance(other, self.__class__)
            and self.chemical_space == other.chemical_space
            and self.id == other.id
            and self.reaction_smarts == other.reaction_smarts
        )
