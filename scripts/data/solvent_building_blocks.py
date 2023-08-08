from pathlib import Path
from tap import tapify 

import numpy as np
import pandas as pd
import torch
from chemfunc import compute_fingerprints



def combine_features(building_blocks_csv: Path, solvent_csv: Path, smiles_columns: list[str], save_dir: Path):
    building_blocks_col, solvent_col = smiles_columns
    building_blocks = pd.read_csv(building_blocks_csv)
    building_block_smiles = building_blocks[building_blocks_col].tolist()

    solvents = pd.read_csv(solvent_csv)

    solvents_smiles = solvents[solvent_col].unique().tolist()
    print(f'num_solvents: {len(solvents_smiles)}')
    solvent_fingerprint= compute_fingerprints(solvents_smiles, "morgan")
    print(f"solvent_features: {solvent_fingerprint.shape}")
    block_fingerprint = compute_fingerprints(building_block_smiles, "morgan")
    
    
    block_features = np.repeat(block_fingerprint, len(solvent_fingerprint), axis=0)  
    print(f'block_shape: {block_features.shape}')
    solvent_features = np.tile(solvent_fingerprint, (block_fingerprint.shape[0], 1)) 
    print(f'tile_shape: {solvent_features.shape}')
    block_solvent_features = np.concatenate((block_features, solvent_features), axis=1)
    np.savez_compressed(f'{save_dir}block_solvent_features.npz', block_solvent_features)


if __name__ == '__main__':
    tapify(combine_features)



