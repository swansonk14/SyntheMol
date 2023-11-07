"""Contains training and predictions functions for Chemprop models."""
from pathlib import Path

import numpy as np
import torch
from chemprop.args import TrainArgs
from chemprop.models import MoleculeModel
from chemprop.utils import load_checkpoint, load_scalers
from sklearn.preprocessing import StandardScaler


def chemprop_build_model(
        dataset_type: str,
        use_rdkit_fingerprints: bool = False,
        rdkit_fingerprint_size: int = 200,
        property_name: str = "task"
) -> MoleculeModel:
    """Builds a Chemprop model.

    :param dataset_type: The type of dataset (classification or regression).
    :param use_rdkit_fingerprints: Whether to use RDKit fingerprints as features.
    :param rdkit_fingerprint_size: The size of the RDKit fingerprint vector.
    :param property_name: The name of the property being predicted.
    :return: A Chemprop model.
    """
    arg_list = [
       '--data_path', 'foo.csv',
       '--dataset_type', dataset_type,
       '--save_dir', 'foo',
       '--quiet'
   ]

    if use_rdkit_fingerprints:
        arg_list += ['--features_generator', 'rdkit_2d_normalized', '--no_features_scaling']

    args = TrainArgs().parse_args(arg_list)
    args.task_names = [property_name]

    if use_rdkit_fingerprints:
        args.features_size = rdkit_fingerprint_size

    # Ensure reproducibility
    torch.manual_seed(0)

    # Build model
    model = MoleculeModel(args)

    return model


def chemprop_load(
        model_path: Path,
        device: torch.device = torch.device('cpu')
) -> MoleculeModel:
    """Loads a Chemprop model.

    :param model_path: A path to a Chemprop model.
    :param device: The device on which to load the model.
    :return: A Chemprop model.
    """
    return load_checkpoint(
        path=str(model_path),
        device=device
    ).eval()


def chemprop_load_scaler(
        model_path: Path
) -> StandardScaler:
    """Loads a Chemprop model's data scaler.

    :param model_path: A path to a Chemprop model.
    :return: A data scaler.
    """
    return load_scalers(path=str(model_path))[0]


def chemprop_predict_on_molecule(
        model: MoleculeModel,
        smiles: str,
        fingerprint: np.ndarray | None = None,
        scaler: StandardScaler | None = None
) -> float:
    """Predicts the property of a molecule using a Chemprop model.

    :param model: A Chemprop model.
    :param smiles: A SMILES string.
    :param fingerprint: A 1D array of molecular fingerprints (if applicable).
    :param scaler: A data scaler (if applicable).
    :return: The prediction on the molecule.
    """
    # Make prediction
    pred = model(
        batch=[[smiles]],
        features_batch=[fingerprint] if fingerprint is not None else None
    ).item()

    # Scale prediction if applicable
    if scaler is not None:
        pred = scaler.inverse_transform([[pred]])[0][0]

    return float(pred)


def chemprop_predict_on_molecule_ensemble(
        models: list[MoleculeModel],
        smiles: str,
        fingerprint: np.ndarray | None = None,
        scalers: list[StandardScaler] | None = None
) -> float:
    """Predicts the property of a molecule using an ensemble of Chemprop models.

    :param models: An ensemble of Chemprop models.
    :param smiles: A SMILES string.
    :param fingerprint: A 1D array of molecular fingerprints (if applicable).
    :param scalers: An ensemble of data scalers (if applicable).
    :return: The ensemble prediction on the molecule.
    """
    return float(np.mean([
        chemprop_predict_on_molecule(
            model=model,
            smiles=smiles,
            fingerprint=fingerprint,
            scaler=scaler
        ) for model, scaler in zip(models, scalers)
    ]))
