"""SMARTS representations custom reactions."""
from synthemol.reactions.query_mol import QueryMol
from synthemol.reactions.reaction import Reaction

# To use custom chemical reactions instead of Enamine REAL reactions, replace None with a list of Reaction objects.
# If CUSTOM_REACTIONS is None, synthemol will default to the reactions in real.py.
CUSTOM_REACTIONS: tuple[Reaction] | None = None
