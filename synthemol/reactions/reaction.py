"""Reaction class and helper functions."""
from typing import Any, Union

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
        post_reactions: tuple["Reaction", ...] = (),
    ) -> None:
        """Initializes the Reaction.

        :param reactants: A list of QueryMols containing the reactants of the reaction.
        :param product: A QueryMol containing the product of the reaction.
        :param chemical_space: The chemical space of the reaction (e.g., Enamine or WuXi).
        :param reaction_id: The ID of the reaction.
        :param sub_reaction_id: The ID of the sub-reaction.
        :param post_reactions: A list of Reactions to run after the main reaction on the product
            (e.g., BOC cleavage or ester hydrolysis).
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
        self.post_reactions = post_reactions

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

    def has_match(self, reactants: list[MOLECULE_TYPE]) -> bool:
        """Determines whether the provided reactants match the reaction.

        :param reactants: A list of reactants in order.
        :return: True if the reactants match the reaction, otherwise False.
        """
        # Ensure the number of reactants matches
        if len(reactants) != len(self.reactants):
            return False

        # Ensure each reactant matches the corresponding QueryMol
        return all(
            self.reactants[reactant_index].has_match(reactant)
            for reactant_index, reactant in enumerate(reactants)
        )

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

        # Convert product mols to unique SMILES (and remove Hs)
        products = list(
            dict.fromkeys(
                Chem.MolToSmiles(Chem.RemoveHs(product_mol[0]))
                for product_mol in product_mols
            )
        )

        # Run any post-reactions on the products
        for post_reaction in self.post_reactions:
            # Run each product through the post-reaction and keep post-product if any, otherwise keep original product
            final_products = []
            for product in products:
                # Continue applying post-reaction until it no longer matches (e.g., in case of multiple Boc groups)
                post_reaction_loops = 0
                while post_reaction.has_match([product]):
                    post_products = post_reaction.run_reactants([product])

                    # If there is only one post-product, use it, otherwise keep the original product
                    if len(post_products) == 1:
                        product = post_products[0]
                        post_reaction_loops += 1
                    elif post_reaction_loops > 50:
                        raise ValueError("Post-reaction looped too many times")
                    else:
                        break

                # Add the final product to the list
                final_products.append(product)

            # Remove duplicate SMILES
            products = list(dict.fromkeys(final_products))

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
