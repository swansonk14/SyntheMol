"""Trains a machine learning classifier model."""
import pickle
from pathlib import Path
from random import Random
from typing import Literal, Optional

import numpy as np
import pandas as pd
import torch
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.neural_network import MLPClassifier
from tap import Tap
from tqdm import trange

from chem_utils.molecular_fingerprints import compute_fingerprints
from chemprop.args import TrainArgs
from chemprop.data import MoleculeDataLoader, MoleculeDatapoint, MoleculeDataset
from chemprop.train import get_loss_func, predict as _chemprop_predict, train as chemprop_train
from chemprop.models import MoleculeModel
from chemprop.utils import build_optimizer, build_lr_scheduler, load_checkpoint, save_checkpoint


class Args(Tap):
    data_path: Path  # Path to CSV file containing data.
    save_dir: Path  # Path to a directory where the trained model(s) and results will be saved.
    model_type: Literal['rf', 'mlp', 'chemprop']  # Type of model to train. 'rf' = random forest. 'mlp' = multilayer perceptron.
    fingerprint_type: Optional[Literal['morgan', 'rdkit']] = None  # Type of fingerprints to use as input features.
    smiles_column: str = 'smiles'  # The name of the column containing SMILES.
    activity_column: str = 'activity'  # The name of the column containing binary activity values.
    num_models: int = 1  # The number of models to train using 10-fold cross-validation.

    def process_args(self) -> None:
        self.save_dir.mkdir(parents=True, exist_ok=True)


def chemprop_predict(*args, **kwargs) -> np.ndarray:
    """Converts chemprop predictions from list of lists to a 1D numpy array (assumes single prediction task)."""
    return np.array(_chemprop_predict(*args, **kwargs))[:, 0]


def build_sklearn_model(train_fingerprints: np.ndarray,
                        test_fingerprints: np.ndarray,
                        train_activities: list[int],
                        test_activities: list[int],
                        save_path: Path,
                        model_type: Literal['rf', 'mlp']) -> dict[str, float]:
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

    # Save model
    with open(save_path, 'wb') as f:
        pickle.dump(model, f)

    # Make predictions
    test_probs = model.predict_proba(test_fingerprints)[:, 1]

    # Evaluate predictions
    scores = {
        'roc_auc': roc_auc_score(test_activities, test_probs),
        'prc_auc': average_precision_score(test_activities, test_probs)
    }

    return scores


def build_chemprop_data_loader(smiles: list[str],
                               fingerprints: Optional[np.ndarray],
                               activities: Optional[list[int]] = None,
                               shuffle: bool = False) -> MoleculeDataLoader:
    """Builds a chemprop MoleculeDataLoader."""
    if fingerprints is None:
        fingerprints = [None] * len(smiles)

    if activities is None:
        activities = [None] * len(smiles)
    else:
        activities = [[float(activity)] for activity in activities]

    return MoleculeDataLoader(
        dataset=MoleculeDataset([
            MoleculeDatapoint(
                smiles=[smiles],
                targets=activity,
                features=fingerprint,
            ) for smiles, fingerprint, activity in zip(smiles, fingerprints, activities)
        ]),
        num_workers=0,  # Needed for deterministic behavior and faster when training/testing on CPU only
        shuffle=shuffle
    )


def build_chemprop_model(train_smiles: list[str],
                         val_smiles: list[str],
                         test_smiles: list[str],
                         fingerprint_type: Optional[str],
                         train_fingerprints: Optional[np.ndarray],
                         val_fingerprints: Optional[np.ndarray],
                         test_fingerprints: Optional[np.ndarray],
                         train_activities: list[int],
                         val_activities: list[int],
                         test_activities: list[int],
                         save_path: Path) -> dict[str, float]:
    """Trains, evaluates, and saves a chemprop model."""
    # Create args
    arg_list = [
        '--data_path', 'foo.csv',
        '--dataset_type', 'classification',
        '--save_dir', 'foo',
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
    args.task_names = ['activity']
    args.train_data_size = len(train_smiles)

    if fingerprint_type is not None:
        args.features_size = train_fingerprints.shape[1]

    # Ensure reproducibility
    torch.manual_seed(0)
    torch.use_deterministic_algorithms(True)

    # Build data loaders
    train_data_loader = build_chemprop_data_loader(
        smiles=train_smiles,
        fingerprints=train_fingerprints,
        activities=train_activities,
        shuffle=True
    )
    val_data_loader = build_chemprop_data_loader(
        smiles=val_smiles,
        fingerprints=val_fingerprints,
        activities=val_activities
    )
    test_data_loader = build_chemprop_data_loader(
        smiles=test_smiles,
        fingerprints=test_fingerprints,
        activities=test_activities
    )

    # Build model
    model = MoleculeModel(args)
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
        n_iter = chemprop_train(
            model=model,
            data_loader=train_data_loader,
            loss_func=loss_func,
            optimizer=optimizer,
            scheduler=scheduler,
            args=args,
            n_iter=n_iter
        )

        val_probs = chemprop_predict(model=model, data_loader=val_data_loader)
        val_score = average_precision_score(val_activities, val_probs)

        if val_score > best_score:
            best_score, best_epoch = val_score, epoch
            save_checkpoint(path=save_path, model=model, args=args)

    # Evaluate on test set using model with the best validation score
    print(f'Best validation {args.metric} = {best_score:.6f} on epoch {best_epoch}')
    model = load_checkpoint(save_path, device=args.device)

    # Make predictions
    test_probs = chemprop_predict(model=model, data_loader=test_data_loader)

    # Evaluate predictions
    scores = {
        'roc_auc': roc_auc_score(test_activities, test_probs),
        'prc_auc': average_precision_score(test_activities, test_probs)
    }

    return scores


def train_model(args: Args) -> None:
    """Trains a machine learning classifier model."""
    # Check compatibility of model and fingerprint type
    if args.model_type != 'chemprop' and args.fingerprint_type is None:
        raise ValueError('Must define fingerprint_type if using sklearn model.')

    # Load data
    data = pd.read_csv(args.data_path)
    print(f'Data size = {len(data):,}')

    # Compute fingerprints
    if args.fingerprint_type is not None:
        fingerprints = compute_fingerprints(data[args.smiles_column], fingerprint_type=args.fingerprint_type)
    else:
        fingerprints = None

    # Set up cross-validation
    num_folds = 10
    indices = np.tile(np.arange(num_folds), 1 + len(data) // num_folds)[:len(data)]
    random = Random(0)
    random.shuffle(indices)

    assert 1 <= args.num_models <= num_folds

    # Run cross-validation
    all_scores = []
    for model_num in trange(args.num_models, desc='cross-val'):
        print(f'Model {model_num}')

        test_index = model_num
        val_index = (model_num + 1) % num_folds

        test_mask = indices == test_index
        val_mask = indices == val_index
        train_mask = ~(test_mask | val_mask)

        test_data = data[test_mask]
        val_data = data[val_mask]
        train_data = data[train_mask]

        if fingerprints is not None:
            test_fingerprints = fingerprints[test_mask]
            val_fingerprints = fingerprints[val_mask]
            train_fingerprints = fingerprints[train_mask]
        else:
            test_fingerprints = val_fingerprints = train_fingerprints = None

        # Build and train model
        if args.model_type == 'chemprop':
            scores = build_chemprop_model(
                train_smiles=train_data[args.smiles_column],
                val_smiles=val_data[args.smiles_column],
                test_smiles=test_data[args.smiles_column],
                fingerprint_type=args.fingerprint_type,
                train_fingerprints=train_fingerprints,
                val_fingerprints=val_fingerprints,
                test_fingerprints=test_fingerprints,
                train_activities=train_data[args.activity_column],
                val_activities=val_data[args.activity_column],
                test_activities=test_data[args.activity_column],
                save_path=args.save_dir / f'model_{model_num}.pt'
            )
        else:
            scores = build_sklearn_model(
                train_fingerprints=train_fingerprints,
                test_fingerprints=test_fingerprints,
                train_activities=train_data[args.activity_column],
                test_activities=test_data[args.activity_column],
                save_path=args.save_dir / f'model_{model_num}.pkl',
                model_type=args.model_type
            )

        # Print scores
        for score_name, score_value in scores.items():
            print(f'Test {score_name} = {score_value:.3f}')
        print()

        all_scores.append(scores)

    # Process and save scores
    score_names = list(all_scores[0])

    all_scores = pd.DataFrame(all_scores)
    all_scores['Model'] = [f'Model {model_num}' for model_num in range(args.num_models)]
    all_scores = all_scores[['Model'] + score_names]
    all_scores.to_csv(args.save_dir / 'scores.csv', index=False)

    # Process and save summary scores
    summary_scores = {}
    for score_name in score_names:
        summary_scores[f'{score_name}_mean'] = np.mean(all_scores[score_name])
        summary_scores[f'{score_name}_std'] = np.std(all_scores[score_name])

    summary_scores = pd.DataFrame([summary_scores])
    summary_scores.to_csv(args.save_dir / 'summary_scores.csv', index=False)


if __name__ == '__main__':
    train_model(Args().parse_args())
