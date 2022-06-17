"""Functions to compute Morgan fingerprints for molecules."""
from functools import partial
from multiprocessing import Pool
from typing import Callable, Union

import numpy as np
from descriptastorus.descriptors import rdNormalizedDescriptors
from rdkit import Chem
from rdkit.DataStructs import ConvertToNumpyArray
from rdkit.Chem import AllChem
from tqdm import tqdm


Molecule = Union[str, Chem.Mol]
FingerprintGenerator = Callable[[Molecule], np.ndarray]
FINGERPRINT_GENERATOR_REGISTRY = {}
MORGAN_RADIUS = 2
MORGAN_NUM_BITS = 2048


def register_fingerprint_generator(fingerprint_type: str) -> Callable[[FingerprintGenerator], FingerprintGenerator]:
    """Creates a decorator which registers a fingerprint generator in a global dictionary to enable access by name.

    :param fingerprint_type: The name to use to access the fingerprint generator.
    :return: A decorator which will add a fingerprint generator to the registry using the specified name.
    """
    def decorator(fingerprint_generator: FingerprintGenerator) -> FingerprintGenerator:
        FINGERPRINT_GENERATOR_REGISTRY[fingerprint_type] = fingerprint_generator
        return fingerprint_generator

    return decorator


def get_fingerprint_generator(fingerprint_type: str) -> FingerprintGenerator:
    """Gets a registered fingerprint generator by name.

    :param fingerprint_type: The name of the fingerprint generator.
    :return: The desired fingerprint generator.
    """
    if fingerprint_type not in FINGERPRINT_GENERATOR_REGISTRY:
        raise ValueError(f'Features generator "{fingerprint_type}" could not be found. '
                         f'If this generator relies on RDKit features, you may need to install descriptastorus.')

    return FINGERPRINT_GENERATOR_REGISTRY[fingerprint_type]


def get_available_fingerprint_generators() -> list[str]:
    """Returns a list of names of available fingerprint generators."""
    return list(FINGERPRINT_GENERATOR_REGISTRY.keys())


# TODO: load all of this from chem_utils instead
@register_fingerprint_generator('morgan')
def compute_morgan_fingerprint(mol: Union[str, Chem.Mol],
                               radius: int = MORGAN_RADIUS,
                               num_bits: int = MORGAN_NUM_BITS) -> np.ndarray:
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


@register_fingerprint_generator('rdkit')
def compute_rdkit_fingerprint(mol: Molecule) -> np.ndarray:
    """Generates RDKit 2D normalized features for a molecule.

    :param mol: A molecule (i.e., either a SMILES or an RDKit molecule).
    :return: A 1D numpy array containing the RDKit 2D normalized features.
    """
    smiles = Chem.MolToSmiles(mol, isomericSmiles=True) if type(mol) != str else mol
    generator = rdNormalizedDescriptors.RDKit2DNormalized()
    rdkit_fp = generator.process(smiles)[1:]
    rdkit_fp = np.where(np.isnan(rdkit_fp), 0, rdkit_fp)
    rdkit_fp = rdkit_fp.astype(np.float32)

    return rdkit_fp


def compute_fingerprint(mol: Molecule, fingerprint_type: str) -> np.ndarray:
    """Generates a molecular fingerprint for a molecule.

    :param mol: A molecule (i.e., either a SMILES string or an RDKit molecule).
    :param fingerprint_type: THe type of fingerprint to compute.
    :return: A 1D numpy array (num_features) containing the fingerprint for the molecule.
    """
    fingerprint_generator = get_fingerprint_generator(fingerprint_type)
    fingerprint = fingerprint_generator(mol)

    return fingerprint


def compute_fingerprints(mols: list[Molecule], fingerprint_type: str) -> np.ndarray:
    """Generates molecular fingerprints for each molecule in a list of molecules (in parallel).

    :param mols: A list of molecules (i.e., either a SMILES string or an RDKit molecule).
    :param fingerprint_type: The type of fingerprint to compute.
    :return: A 2D numpy array (num_molecules, num_features) containing the fingerprints for each molecule.
    """
    fingerprint_generator = get_fingerprint_generator(fingerprint_type)

    with Pool() as pool:
        fingerprints = np.array(list(tqdm(pool.imap(fingerprint_generator, mols),
                                          total=len(mols), desc=f'{fingerprint_type} fingerprints')))

    return fingerprints
