"""Contains constants shared throughout synthemol."""
from importlib import resources
from typing import Literal

from rdkit.Chem import Mol
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.neural_network import MLPClassifier, MLPRegressor


CHEMBL_SMILES_COL = 'Smiles'
REAL_SPACE_SIZE = 31507987117  # As of August 30, 2022, in the 2022 q1-2 REAL space
REAL_REACTION_COL = 'reaction'
REAL_BUILDING_BLOCK_COLS = ['reagent1', 'reagent2', 'reagent3', 'reagent4']
REAL_BUILDING_BLOCK_ID_COL = 'reagent_id'
REAL_TYPE_COL = 'Type'
REAL_SMILES_COL = 'smiles'
SMILES_COL = 'smiles'
ACTIVITY_COL = 'activity'
SCORE_COL = 'score'
MODEL_TYPE = Literal['random_forest', 'mlp', 'chemprop']
FINGERPRINT_TYPES = Literal['morgan', 'rdkit']
MOLECULE_TYPE = str | Mol  # Either a SMILES string or an RDKit Mol object
SKLEARN_MODEL_TYPES = RandomForestClassifier | RandomForestRegressor | MLPClassifier | MLPRegressor
SKLEARN_MODEL_NAME_TYPES = Literal['random_forest', 'mlp']
MODEL_TYPES = Literal['random_forest', 'mlp', 'chemprop']
DATASET_TYPES = Literal['classification', 'regression']
OPTIMIZATION_TYPES = Literal['maximize', 'minimize']

# Path where data files are stored
with resources.path('synthemol', 'resources') as resources_dir:
    DATA_DIR = resources_dir

# If using custom building blocks, replace BUILDING_BLOCKS_PATH and set REACTION_TO_BUILDING_BLOCKS_PATH to None
BUILDING_BLOCKS_PATH = DATA_DIR / 'real' / 'building_blocks.csv'
REACTION_TO_BUILDING_BLOCKS_PATH = DATA_DIR / 'real' / 'reaction_to_building_blocks.pkl'
