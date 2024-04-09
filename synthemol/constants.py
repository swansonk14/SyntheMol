"""Contains constants shared throughout synthemol."""
from importlib import resources
from typing import Literal

from rdkit.Chem import Mol
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor


CHEMICAL_SPACES = Literal["real", "wuxi", "custom"]
CHEMBL_SMILES_COL = "Smiles"
REAL_SPACE_SIZE = 31507987117  # As of August 30, 2022, in the 2022 q1-2 REAL space
REAL_REACTION_COL = "reaction"
REAL_BUILDING_BLOCK_COLS = ["reagent1", "reagent2", "reagent3", "reagent4"]
REAL_BUILDING_BLOCK_ID_COL = "reagent_id"
REAL_TYPE_COL = "Type"
REAL_SMILES_COL = "smiles"
REAL_ID_COL = "id"
WUXI_GALAXI_SIZE = 16146071436  # As of the December 31, 2022, version of WuXi GalaXi
WUXI_SMILES_COL = "smiles"
WUXI_ID_COL = "id"
WUXI_BUILDING_BLOCK_SMILES_COL = "smiles"
WUXI_BUILDING_BLOCK_ID_COL = "wxid"
WUXI_BUILDING_BLOCK_SUBSET_COL = "subset"
WUXI_CORES_CORE_COL = "CORE"
WUXI_CORES_SMILES_COL = "SMILES"
WUXI_CORES_ID_COL = "WXID"
ID_COL = "reagent_id"
SMILES_COL = "smiles"
ACTIVITY_COL = "activity"
SCORE_COL = "score"
ROLLOUT_COL = "rollout_num"
FINGERPRINT_TYPES = Literal["rdkit", "morgan"]
MOLECULE_TYPE = str | Mol  # Either a SMILES string or an RDKit Mol object
SKLEARN_MODEL_TYPES = RandomForestClassifier | RandomForestRegressor
SKLEARN_MODEL_NAME_TYPES = Literal["random_forest"]
MODEL_TYPES = Literal["random_forest", "chemprop"]
SCORE_TYPES = Literal["random_forest", "chemprop", "qed", "clogp"]
DATASET_TYPES = Literal["classification", "regression"]
RL_MODEL_TYPES = Literal["mlp_rdkit", "chemprop", "chemprop_rdkit"]
RL_PREDICTION_TYPES = Literal["classification", "regression"]
OPTIMIZATION_TYPES = Literal["maximize", "minimize"]
OLD_REACTION_ORDER = [
    275592,
    22,
    11,
    527,
    2430,
    2708,
    240690,
    2230,
    2718,
    40,
    1458,
    271948,
    27,
]
OLD_REACTIONS = set(OLD_REACTION_ORDER)

# Path where data files are stored
with resources.path("synthemol", "resources") as resources_dir:
    DATA_DIR = resources_dir

# If using custom building blocks, replace BUILDING_BLOCKS_PATH and set REACTION_TO_BUILDING_BLOCKS_PATH to None
BUILDING_BLOCKS_PATH = DATA_DIR / "real" / "building_blocks.csv"
REACTION_TO_BUILDING_BLOCKS_PATH = DATA_DIR / "real" / "reaction_to_building_blocks.pkl"
