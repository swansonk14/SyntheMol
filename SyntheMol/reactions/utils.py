"""Utility functions for SyntheMol reactions."""
from SyntheMol.reactions import Reaction



def set_allowed_smiles(
        reactions: list[Reaction],
        reaction_to_reactant_to_building_blocks: dict[int, dict[int, set[str]]]
) -> None:
    """Sets the allowed SMILES for each reactant in each Reaction in a list of Reactions.

    Note: Modifies Reactions in place.

    :param reactions: A list of Reactions whose allowed SMILES will be set.
    :param reaction_to_reactant_to_building_blocks: A dictionary mapping from reaction ID
                                                    to reactant index to a set of allowed SMILES.
    """
    for reaction in reactions:
        for reactant_index, reactant in enumerate(reaction.reactants):
            reactant.allowed_smiles = reaction_to_reactant_to_building_blocks[reaction.id][reactant_index]
