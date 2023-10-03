"""Contains constants shared throughout synthemol."""
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
WUXI_BUILDING_BLOCK_SMILES_COL = 'smiles'
WUXI_BUILDING_BLOCK_ID_COL = 'wxid'
WUXI_BUILDING_BLOCK_SUBSET_COL = 'subset'
WUXI_CORES_CORE_COL = 'CORE'
WUXI_CORES_SMILES_COL = 'SMILES'
WUXI_CORES_ID_COL = 'WXID'
ID_COL = 'ID'
SMILES_COL = 'smiles'
ACTIVITY_COL = 'activity'
SCORE_COL = 'score'
ROLLOUT_COL = 'rollout_num'
MODEL_TYPE = Literal['random_forest', 'mlp', 'chemprop']
FINGERPRINT_TYPES = Literal['none', 'morgan', 'rdkit']
MOLECULE_TYPE = str | Mol  # Either a SMILES string or an RDKit Mol object
SKLEARN_MODEL_TYPES = RandomForestClassifier | RandomForestRegressor | MLPClassifier | MLPRegressor
SKLEARN_MODEL_NAME_TYPES = Literal['random_forest', 'mlp']
MODEL_TYPES = Literal['random_forest', 'mlp', 'chemprop', 'qed']
DATASET_TYPES = Literal['classification', 'regression']
RL_MODEL_TYPES = Literal['rdkit', 'chemprop_pretrained', 'chemprop_scratch']
RL_PREDICTION_TYPES = Literal['classification', 'regression']
OPTIMIZATION_TYPES = Literal['maximize', 'minimize']

# Path where data files are stored
DATA_DIR = Path(__file__).parent / 'files'

# If using custom building blocks, replace BUILDING_BLOCKS_PATH and set REACTION_TO_BUILDING_BLOCKS_PATH to None
BUILDING_BLOCKS_PATH = DATA_DIR / 'building_blocks.csv'
REACTION_TO_BUILDING_BLOCKS_PATH = DATA_DIR / 'reaction_to_building_blocks.pkl'
