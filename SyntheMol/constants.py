"""Contains constants shared throughout SyntheMol."""
from typing import Literal

CHEMBL_SMILES_COL = 'Smiles'
REAL_SPACE_SIZE = 31507987117  # As of August 30, 2022, in the 2022 q1-2 REAL space
REAL_REACTION_COL = 'reaction'
REAL_BUILDING_BLOCK_COLS = ['reagent1', 'reagent2', 'reagent3', 'reagent4']
REAL_BUILDING_BLOCK_ID_COL = 'Reagent_ID'
REAL_TYPE_COL = 'Type'
REAL_SMILES_COL = 'smiles'
SMILES_COL = 'smiles'
MODEL_TYPE = Literal['rf', 'mlp', 'chemprop']
FINGERPRINT_TYPE = Literal['morgan', 'rdkit'] | None
