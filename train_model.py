"""Trains a machine learning classifier model."""
import pickle
from pathlib import Path
from typing import Literal

import numpy as np
import pandas as pd
import torch
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.model_selection import train_test_split
from sklearn.neural_network import MLPClassifier
from tap import Tap
from torch.optim.lr_scheduler import ExponentialLR
from tqdm import trange

from chem_utils.molecular_fingerprints import compute_fingerprints
from chemprop.args import TrainArgs
from chemprop.data import MoleculeDataLoader, MoleculeDatapoint, MoleculeDataset
from chemprop.train import get_loss_func, predict, train
from chemprop.models import MoleculeModel
from chemprop.utils import build_optimizer, build_lr_scheduler, load_checkpoint, save_checkpoint


class Args(Tap):
    data_path: Path  # Path to CSV file containing data.
    save_path: Path  # Path to a PKL or PT file where the trained model will be saved.
    model_type: Literal['rf', 'mlp', 'chemprop']  # Type of model to train. 'rf' = random forest. 'mlp' = multilayer perceptron.
    fingerprint_type: Literal['morgan', 'rdkit']  # Type of fingerprints to use as input features.
    smiles_column: str = 'smiles'  # The name of the column containing SMILES.
    activity_column: str = 'activity'  # The name of the column containing binary activity values.

    def process_args(self) -> None:
        self.save_path.parent.mkdir(parents=True, exist_ok=True)


def build_sklearn_model(train_fingerprints: np.ndarray,
                        test_fingerprints: np.ndarray,
                        train_activities: list[int],
                        test_activities: list[int],
                        save_path: Path,
                        model_type: Literal['rf', 'mlp']) -> None:
    """Trains, evaluates, and saves a scikit-learn model."""
    # Build model
    if model_type == 'rf':
        model = RandomForestClassifier(n_jobs=-1, random_state=0)
    elif model_type == 'mlp':
        model = MLPClassifier(hidden_layer_sizes=(100, 100, 100), random_state=0)
    else:
        raise ValueError(f'Model type "{model_type}" is not supported.')

    print(model)

    # Train model
    model.fit(train_fingerprints, train_activities)

    # Evaluate model
    test_probs = model.predict_proba(test_fingerprints)[:, 1]
    print(f'Test ROC-AUC = {roc_auc_score(test_activities, test_probs):.3f}')
    print(f'Test PRC-AUC = {average_precision_score(test_activities, test_probs):.3f}')

    # Save model
    with open(save_path, 'wb') as f:
        pickle.dump(model, f)


def build_chemprop_model(train_smiles: list[str],
                         val_smiles: list[str],
                         test_smiles: list[str],
                         train_fingerprints: np.ndarray,
                         val_fingerprints: np.ndarray,
                         test_fingerprints: np.ndarray,
                         train_activities: list[int],
                         val_activities: list[int],
                         test_activities: list[int],
                         save_path: Path) -> None:
    """Trains, evaluates, and saves a chemprop model."""
    # Create args
    args = TrainArgs().parse_args(['--data_path', 'foo.csv', '--dataset_type', 'classification', '--quiet'])
    args.task_names = ['activity']
    args.train_data_size = len(train_smiles)

    # Build data loaders
    train_data_loader = MoleculeDataLoader(
        dataset=MoleculeDataset([
            MoleculeDatapoint(
                smiles=[smiles],
                targets=[float(activity)],
                features=fingerprint,
            ) for smiles, fingerprint, activity in zip(train_smiles, train_fingerprints, train_activities)
        ]),
        batch_size=args.batch_size,
        num_workers=0,
        shuffle=True,
        seed=args.seed
    )
    val_data_loader = MoleculeDataLoader(
        dataset=MoleculeDataset([
            MoleculeDatapoint(
                smiles=[smiles],
                targets=[float(activity)],
                features=fingerprint,
            ) for smiles, fingerprint, activity in zip(val_smiles, val_fingerprints, val_activities)
        ]),
        batch_size=args.batch_size,
        num_workers=0,
        shuffle=False,
        seed=args.seed
    )
    test_data_loader = MoleculeDataLoader(
        dataset=MoleculeDataset([
            MoleculeDatapoint(
                smiles=[smiles],
                targets=[float(activity)],
                features=fingerprint,
            ) for smiles, fingerprint, activity in zip(test_smiles, test_fingerprints, test_activities)
        ]),
        batch_size=args.batch_size,
        num_workers=0,
        shuffle=False,
        seed=args.seed
    )

    # Build model
    torch.manual_seed(args.seed)
    model = MoleculeModel(args)
    model = model.to(args.device)
    print(model)

    # Get loss function, optimizer, and learning rate scheduler
    loss_func = get_loss_func(args)
    optimizer = build_optimizer(model, args)
    scheduler = build_lr_scheduler(optimizer, args)

    # Run training
    best_score = float('inf') if args.minimize_score else -float('inf')
    best_epoch = n_iter = 0
    for epoch in trange(args.epochs):
        print(f'Epoch {epoch}')
        n_iter = train(
            model=model,
            data_loader=train_data_loader,
            loss_func=loss_func,
            optimizer=optimizer,
            scheduler=scheduler,
            args=args,
            n_iter=n_iter
        )

        if isinstance(scheduler, ExponentialLR):
            scheduler.step()

        val_probs = predict(
            model=model,
            data_loader=val_data_loader
        )
        val_score = average_precision_score(val_activities, val_probs)

        if val_score > best_score:
            best_score, best_epoch = val_score, epoch
            save_checkpoint(
                path=save_path,
                model=model,
                args=args
            )

    # Evaluate on test set using model with best validation score
    print(f'Best validation {args.metric} = {best_score:.6f} on epoch {best_epoch}')
    model = load_checkpoint(save_path, device=args.device)

    test_probs = predict(
        model=model,
        data_loader=test_data_loader
    )
    print(f'Test ROC-AUC = {roc_auc_score(test_activities, test_probs):.3f}')
    print(f'Test PRC-AUC = {average_precision_score(test_activities, test_probs):.3f}')


def train_model(args: Args) -> None:
    """Trains a machine learning classifier model."""
    # Load data
    data = pd.read_csv(args.data_path)
    print(f'Data size = {len(data):,}')

    # Get SMILES
    smiles = list(data[args.smiles_column])

    # Get fingerprints
    fingerprints = compute_fingerprints(smiles, fingerprint_type=args.fingerprint_type)

    # Get activity
    activities = list(data[args.activity_column])

    # Split into train and test
    train_smiles, val_test_smiles, train_fingerprints, val_test_fingerprints, train_activities, val_test_activities = \
        train_test_split(
            smiles,
            fingerprints,
            activities,
            test_size=0.2,
            random_state=0
        )
    val_smiles, test_smiles, val_fingerprints, test_fingerprints, val_activities, test_activities = \
        train_test_split(
            val_test_smiles,
            val_test_fingerprints,
            val_test_activities,
            test_size=0.5,
            random_state=0
        )

    print(f'Train size = {len(train_smiles):,}')
    print(f'Validation size = {len(val_smiles):,}')
    print(f'Test size = {len(test_smiles):,}')

    # Build and train model
    if args.model_type == 'chemprop':
        build_chemprop_model(
            train_smiles=train_smiles,
            val_smiles=val_smiles,
            test_smiles=test_smiles,
            train_fingerprints=train_fingerprints,
            val_fingerprints=val_fingerprints,
            test_fingerprints=test_fingerprints,
            train_activities=train_activities,
            val_activities=val_activities,
            test_activities=test_activities,
            save_path=args.save_path
        )
    else:
        build_sklearn_model(
            train_fingerprints=train_fingerprints,
            test_fingerprints=test_fingerprints,
            train_activities=train_activities,
            test_activities=test_activities,
            save_path=args.save_path,
            model_type=args.model_type
        )


if __name__ == '__main__':
    train_model(Args().parse_args())
