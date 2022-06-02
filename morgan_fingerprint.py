"""Functions to compute Morgan fingerprints for molecules."""
from functools import partial
from multiprocessing import Pool
from typing import Union

import numpy as np
from rdkit import Chem
from rdkit.DataStructs import ConvertToNumpyArray
from rdkit.Chem import AllChem
from tqdm import tqdm


def compute_morgan_fingerprint(mol: Union[str, Chem.Mol],
                               radius: int = 2,
                               num_bits: int = 2048) -> np.ndarray:
    """
    Generates a binary Morgan fingerprint for a molecule.

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


def compute_morgan_fingerprints(mols: list[Union[str, Chem.Mol]],
                                radius: int = 2,
                                num_bits: int = 2048) -> np.ndarray:
    """
    Generates binary Morgan fingerprints for each molecule in a list of molecules (in parallel).

    :param mols: A list of molecules (i.e., either a SMILES string or an RDKit molecule).
    :param radius: Morgan fingerprint radius.
    :param num_bits: Number of bits in Morgan fingerprint.
    :return: A 2D boolean numpy array (num_molecules, num_bits) containing
             the binary Morgan fingerprints for each molecule.
    """
    morgan_fn = partial(compute_morgan_fingerprint, radius=radius, num_bits=num_bits)

    with Pool() as pool:
        morgan_fps = np.array(list(tqdm(pool.imap(morgan_fn, mols), total=len(mols), desc='Morgan fingerprints')))

    return morgan_fps
