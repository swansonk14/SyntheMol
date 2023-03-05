"""Contains constants shared throughout SyntheMol."""
from typing import Literal

from rdkit.Chem import Mol
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.neural_network import MLPClassifier, MLPRegressor


CHEMBL_SMILES_COL = 'Smiles'
REAL_SPACE_SIZE = 31507987117  # As of August 30, 2022, in the 2022 q1-2 REAL space
REAL_REACTION_COL = 'reaction'
REAL_BUILDING_BLOCK_COLS = ['reagent1', 'reagent2', 'reagent3', 'reagent4']
REAL_BUILDING_BLOCK_ID_COL = 'Reagent_ID'
REAL_TYPE_COL = 'Type'
REAL_SMILES_COL = 'smiles'
SMILES_COL = 'smiles'
SCORE_COL = 'score'
MODEL_TYPE = Literal['rf', 'mlp', 'chemprop']
FINGERPRINT_TYPES = Literal['morgan', 'rdkit'] | None
MOLECULE_TYPE = str | Mol  # Either a SMILES string or an RDKit Mol object
SKLEARN_MODEL_TYPES = RandomForestClassifier | RandomForestRegressor | MLPClassifier | MLPRegressor
SKLEARN_MODEL_NAME_TYPES = Literal['random_forest', 'mlp']
MODEL_TYPES = Literal['random_forest', 'mlp', 'chemprop']
DATASET_TYPES = Literal['classification', 'regression']
