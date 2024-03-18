"""synthemol.reactions package."""
from synthemol.reactions.custom import CUSTOM_REACTIONS
from synthemol.reactions.query_mol import QueryMol
from synthemol.reactions.reaction import Reaction
from synthemol.reactions.real import MOST_COMMON_REAL_REACTIONS, REAL_REACTIONS
from synthemol.reactions.utils import (
    load_and_set_allowed_reaction_building_blocks,
    set_all_building_blocks,
)
from synthemol.reactions.wuxi import WUXI_REACTIONS

CHEMICAL_SPACE_TO_REACTIONS: dict[str, tuple[Reaction, ...]] = {
    "real": REAL_REACTIONS,
    "real_most_common": MOST_COMMON_REAL_REACTIONS,
    "wuxi": WUXI_REACTIONS,
    "custom": CUSTOM_REACTIONS,
}
