"""Functions to compute Morgan fingerprints for molecules."""
from functools import partial
from multiprocessing import Pool
from typing import Union

import numpy as np
from descriptastorus.descriptors import rdDescriptors, rdNormalizedDescriptors
from rdkit import Chem
from rdkit.DataStructs import ConvertToNumpyArray
from rdkit.Chem import AllChem
from tqdm import tqdm


# TODO: load all of this from chem_utils instead
def compute_morgan_fingerprint(mol: Union[str, Chem.Mol],
                               radius: int = 2,
                               num_bits: int = 2048) -> np.ndarray:
    """Generates a binary Morgan fingerprint for a molecule.

    :param mol: A molecule (i.e., either a SMILES string or an RDKit molecule).
    :param radius: Morgan fingerprint radius.
    :param num_bits: Number of bits in Morgan fingerprint.
    :return: A 1D boolean numpy array (num_bits,) containing the binary Morgan fingerprint.
    """
    mol = Chem.MolFromSmiles(mol) if type(mol) == str else mol
    morgan_vec = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=num_bits)
    morgan_fp = np.zeros((1,))
    ConvertToNumpyArray(morgan_vec, morgan_fp)
    morgan_fp = morgan_fp.astype(bool)

    return morgan_fp


def compute_rdkit_fingerprint(mol: Union[str, Chem.Mol]) -> np.ndarray:
    """Generates RDKit 2D normalized features for a molecule.

    :param mol: A molecule (i.e., either a SMILES or an RDKit molecule).
    :return: A 1D numpy array containing the RDKit 2D normalized features.
    """
    smiles = Chem.MolToSmiles(mol, isomericSmiles=True) if type(mol) != str else mol
    generator = rdNormalizedDescriptors.RDKit2DNormalized()
    features = generator.process(smiles)[1:]

    return features


def compute_fingerprints(mols: list[Union[str, Chem.Mol]], fingerprint_type: str, **kwargs) -> np.ndarray:
    """Generates molecular fingerprints for each molecule in a list of molecules (in parallel).

    :param mols: A list of molecules (i.e., either a SMILES string or an RDKit molecule).
    :param fingerprint_type: The type of fingerprint to compute.
    :return: A 2D numpy array (num_molecules, num_features) containing the fingerprints for each molecule.
    """
    if fingerprint_type == 'morgan':
        fingerprint_fn = compute_morgan_fingerprint
    elif fingerprint_type == 'rdkit':
        fingerprint_fn = compute_rdkit_fingerprint
    else:
        raise ValueError(f'Fingerprint type "{fingerprint_type}" is not supported.')

    fingerprint_fn = partial(fingerprint_fn, **kwargs)

    with Pool() as pool:
        fingerprints = np.array(list(tqdm(pool.imap(fingerprint_fn, mols),
                                          total=len(mols), desc=f'{fingerprint_type} fingerprints')))

    return fingerprints
