"""SMARTS representations custom reactions."""
from SyntheMol.reactions.query_mol import QueryMol
from SyntheMol.reactions.reaction import Reaction

# To use custom chemical reactions instead of Enamine REAL reactions, replace None with a list of Reaction objects.
# If CUSTOM_REACTIONS is None, SyntheMol will default to the reactions in real.py.
CUSTOM_REACTIONS = None
