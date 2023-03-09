"""Contains training and predictions functions for Chemprop models."""
from pathlib import Path

import numpy as np
import torch
from sklearn.metrics import average_precision_score, mean_absolute_error
from tqdm import trange

from chemprop.args import TrainArgs
from chemprop.data import MoleculeDataLoader, MoleculeDatapoint, MoleculeDataset
from chemprop.train import get_loss_func, predict as _chemprop_predict, train as _chemprop_train
from chemprop.models import MoleculeModel
from chemprop.utils import build_optimizer, build_lr_scheduler, load_checkpoint, save_checkpoint


def chemprop_load(
        model_path: Path
) -> MoleculeModel:
    """Loads a Chemprop model.

    :param model_path: A path to a Chemprop model.
    :return: A Chemprop model.
    """
    return load_checkpoint(path=str(model_path)).eval()


def chemprop_predict_ensemble_on_molecule(
        models: list[MoleculeModel],
        smiles: str,
        fingerprint: np.ndarray | None = None
) -> float:
    """Predicts the property of a molecule using an ensemble of Chemprop models.

    :param models: An ensemble of Chemprop models.
    :param smiles: A SMILES string.
    :param fingerprint: A 1D array of molecular fingerprints (if applicable).
    :return: The ensemble prediction on the molecule.
    """
    return float(np.mean([
        model(
            batch=[[smiles]],
            features_batch=[fingerprint] if fingerprint is not None else None
        ).item()
        for model in models
    ]))


def chemprop_predict(
        model: MoleculeModel,
        smiles: list[str],
        fingerprints: np.ndarray | None = None
) -> np.ndarray:
    """Predicts molecular properties using a Chemprop model.

    :param model: A Chemprop model.
    :param smiles: A list of SMILES strings.
    :param fingerprints: A 2D array of molecular fingerprints (num_molecules, num_features).
    :return: A 1D array of predicted properties (num_molecules,)."""
    # Set up data loader
    data_loader = chemprop_build_data_loader(
        smiles=smiles,
        fingerprints=fingerprints
    )

    # Make predictions
    preds = np.array(_chemprop_predict(model=model, data_loader=data_loader))[:, 0]

    return preds


def chemprop_build_data_loader(
        smiles: list[str],
        fingerprints: np.ndarray | None = None,
        properties: list[int] | None = None,
        shuffle: bool = False
) -> MoleculeDataLoader:
    """Builds a chemprop MoleculeDataLoader.

    :param smiles: A list of SMILES strings.
    :param fingerprints: A 2D array of molecular fingerprints (num_molecules, num_features).
    :param properties: A list of molecular properties (num_molecules,).
    :param shuffle: Whether to shuffle the data loader.
    :return: A Chemprop data loader.
    """
    if fingerprints is None:
        fingerprints = [None] * len(smiles)

    if properties is None:
        properties = [None] * len(smiles)
    else:
        properties = [[float(prop)] for prop in properties]

    return MoleculeDataLoader(
        dataset=MoleculeDataset([
            MoleculeDatapoint(
                smiles=[smiles],
                targets=prop,
                features=fingerprint,
            ) for smiles, fingerprint, prop in zip(smiles, fingerprints, properties)
        ]),
        num_workers=0,  # Needed for deterministic behavior and faster when training/testing on CPU only
        shuffle=shuffle
    )


def chemprop_train(
        dataset_type: str,
        train_smiles: list[str],
        val_smiles: list[str],
        fingerprint_type: str | None,
        train_fingerprints: np.ndarray | None,
        val_fingerprints: np.ndarray | None,
        property_name: str,
        train_properties: list[int],
        val_properties: list[int],
        epochs: int,
        save_path: Path
) -> MoleculeModel:
    """Trains and saves a Chemprop model.

    :param dataset_type: The type of dataset (classification or regression).
    :param train_smiles: A list of SMILES strings for training.
    :param val_smiles: A list of SMILES strings for validation.
    :param fingerprint_type: The type of fingerprint to use (morgan or rdkit).
    :param train_fingerprints: A 2D array of molecular fingerprints for training (num_molecules, num_features).
    :param val_fingerprints: A 2D array of molecular fingerprints for validation (num_molecules, num_features).
    :param property_name: The name of the property to predict.
    :param train_properties: A list of molecular properties for training (num_molecules,).
    :param val_properties: A list of molecular properties for validation (num_molecules,).
    :param epochs: The number of epochs to train for.
    :param save_path: The path to save the model to.
    """
    # Create args
    arg_list = [
        '--data_path', 'foo.csv',
        '--dataset_type', dataset_type,
        '--save_dir', 'foo',
        '--epochs', str(epochs),
        '--quiet'
    ]

    match fingerprint_type:
        case 'morgan':
            arg_list += ['--features_generator', 'morgan']
        case 'rdkit':
            arg_list += ['--features_generator', 'rdkit_2d_normalized', '--no_features_scaling']
        case None:
            pass
        case _:
            raise ValueError(f'Fingerprint type "{fingerprint_type}" is not supported.')

    args = TrainArgs().parse_args(arg_list)
    args.task_names = [property_name]
    args.train_data_size = len(train_smiles)

    if fingerprint_type is not None:
        args.features_size = train_fingerprints.shape[1]

    # Ensure reproducibility
    torch.manual_seed(0)
    torch.use_deterministic_algorithms(True)

    # Build data loaders
    train_data_loader = chemprop_build_data_loader(
        smiles=train_smiles,
        fingerprints=train_fingerprints,
        properties=train_properties,
        shuffle=True
    )

    # Build model
    model = MoleculeModel(args)
    print(model)

    # Get loss function, optimizer, and learning rate scheduler
    loss_func = get_loss_func(args)
    optimizer = build_optimizer(model, args)
    scheduler = build_lr_scheduler(optimizer, args)

    # Run training
    save_path = str(save_path)
    best_score = float('inf') if args.minimize_score else -float('inf')
    val_metric = 'PRC-AUC' if dataset_type == 'classification' else 'MAE'
    best_epoch = n_iter = 0
    for epoch in trange(args.epochs):
        print(f'Epoch {epoch}')
        n_iter = _chemprop_train(
            model=model,
            data_loader=train_data_loader,
            loss_func=loss_func,
            optimizer=optimizer,
            scheduler=scheduler,
            args=args,
            n_iter=n_iter
        )

        val_probs = chemprop_predict(
            model=model,
            smiles=val_smiles,
            fingerprints=val_fingerprints
        )

        if dataset_type == 'classification':
            val_score = average_precision_score(val_properties, val_probs)
            new_best_val_score = val_score > best_score
        elif dataset_type == 'regression':
            val_score = mean_absolute_error(val_properties, val_probs)
            new_best_val_score = val_score < best_score
        else:
            raise ValueError(f'Dataset type "{dataset_type}" is not supported.')

        if new_best_val_score:
            best_score, best_epoch = val_score, epoch
            save_checkpoint(path=save_path, model=model, args=args)

    # Evaluate on test set using model with the best validation score
    print(f'Best validation {val_metric} = {best_score:.6f} on epoch {best_epoch}')
    model = load_checkpoint(save_path, device=args.device)

    return model
