"""SyntheMol.reactions package."""
from SyntheMol.reactions.custom import CUSTOM_REACTIONS
from SyntheMol.reactions.query_mol import QueryMol
from SyntheMol.reactions.reaction import Reaction
from SyntheMol.reactions.real import REAL_REACTIONS
from SyntheMol.reactions.utils import (
    load_and_set_allowed_reaction_building_blocks,
    set_all_building_blocks,
    set_allowed_reaction_building_blocks
)

if CUSTOM_REACTIONS is None:
    REACTIONS = REAL_REACTIONS
else:
    REACTIONS = CUSTOM_REACTIONS
