#TODO: Multi-thread this using pool for multiple features?
from cmath import nan
import numpy as np
from rdkit import Chem
from rdkit.DataStructs import ConvertToNumpyArray
from rdkit.Chem import AllChem
from typing import Union
from tqdm import tqdm
from functools import partial
from multiprocessing import Pool
import sys


np.set_printoptions(threshold=sys.maxsize)


def compute_rdkit_2d_normalized_feature(mol: Union[str, Chem.Mol]) -> np.ndarray:

    from descriptastorus.descriptors import rdDescriptors, rdNormalizedDescriptors
    smiles = Chem.MolToSmiles(mol, isomericSmiles=True) if type(mol) != str else mol
    generator = rdNormalizedDescriptors.RDKit2DNormalized()
    feature = generator.process(smiles)[1:]
    feature = np.where(np.isnan(feature) == True, 0, feature)

    
    return feature

def compute_rdkit_2d_normalized_features(mols: list[Union[str, Chem.Mol]]) ->np.ndarray:
    feature = partial(compute_rdkit_2d_normalized_feature)

    with Pool() as pool:
        features = np.array(list(tqdm(pool.imap(feature,mols), total = len(mols), desc = 'RDKit Features')))
    features = features.astype(np.float32)
    return features