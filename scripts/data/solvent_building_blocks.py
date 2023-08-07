from pathlib import Path
from tap import tapify 

import numpy as np
import pandas as pd
import torch
from chemfunc import compute_fingerprints



def combine_features(building_blocks_csv: Path, solvent_csv: Path, smiles_columns: list[str]):
    building_blocks_col, solvent_col = smiles_columns
    building_blocks = pd.read_csv(building_blocks_csv)
    building_block_smiles = building_blocks[building_blocks_col].tolist()

    solvents = pd.read_csv(solvent_csv)
    solvents_smiles = solvents[solvent_col].unique().tolist()

    block_features = compute_fingerprints(building_block_smiles, "morgan")
    solvent_features= compute_fingerprints(solvents_smiles, "morgan")

    block_features = np.repeat(block_features, len(solvent_features), axis=0)  
    solvent_features = np.tile(solvent_features, (len(block_features), 1)) 
    block_solvent_features = np.concatenate(block_features, solvent_features, axis=1)
    np.savez_compressed('block_solvent_features.npz', block_solvent_features)


if __name__ == '__main__':
    tapify(combine_features)



