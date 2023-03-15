"""Contains constants shared throughout SyntheMol."""
from pathlib import Path
from typing import Literal

from rdkit.Chem import Mol
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.neural_network import MLPClassifier, MLPRegressor


CHEMBL_SMILES_COL = 'Smiles'
REAL_SPACE_SIZE = 31507987117  # As of August 30, 2022, in the 2022 q1-2 REAL space
REAL_REACTION_COL = 'reaction'
REAL_BUILDING_BLOCK_COLS = ['reagent1', 'reagent2', 'reagent3', 'reagent4']
REAL_BUILDING_BLOCK_ID_COL = 'ID'
REAL_TYPE_COL = 'Type'
REAL_SMILES_COL = 'smiles'
SMILES_COL = 'smiles'
SCORE_COL = 'score'
MODEL_TYPE = Literal['random_forest', 'mlp', 'chemprop']
FINGERPRINT_TYPES = Literal['morgan', 'rdkit']
MOLECULE_TYPE = str | Mol  # Either a SMILES string or an RDKit Mol object
SKLEARN_MODEL_TYPES = RandomForestClassifier | RandomForestRegressor | MLPClassifier | MLPRegressor
SKLEARN_MODEL_NAME_TYPES = Literal['random_forest', 'mlp']
MODEL_TYPES = Literal['random_forest', 'mlp', 'chemprop']
DATASET_TYPES = Literal['classification', 'regression']

# Path where data files are stored
DATA_DIR = Path(__file__).parent / 'files'

# If using custom building blocks, replace BUILDING_BLOCKS_PATH and set REACTION_TO_BUILDING_BLOCKS_PATH to None
BUILDING_BLOCKS_PATH = DATA_DIR / 'building_blocks.csv'
REACTION_TO_BUILDING_BLOCKS_PATH = DATA_DIR / 'reaction_to_building_blocks.pkl'
