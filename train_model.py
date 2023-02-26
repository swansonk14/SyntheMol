"""Trains a machine learning classifier model."""
import pickle
import time
from pathlib import Path
from random import Random
from typing import Literal, Optional

import numpy as np
import pandas as pd
import torch
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.metrics import average_precision_score, mean_absolute_error, r2_score, roc_auc_score
from sklearn.neural_network import MLPClassifier, MLPRegressor
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
    dataset_type: Literal['classification', 'regression']  # Type of dataset.
    model_type: Literal['rf', 'mlp', 'chemprop']  # Type of model to train. 'rf' = random forest. 'mlp' = multilayer perceptron.
    fingerprint_type: Optional[Literal['morgan', 'rdkit']] = None  # Type of fingerprints to use as input features.
    property_column: str  # The name of the column containing property values.
    smiles_column: str = 'smiles'  # The name of the column containing SMILES.
    num_models: int = 1  # The number of models to train using 10-fold cross-validation.
    epochs: int = 30  # The number of epochs to train the chemprop model.

    def process_args(self) -> None:
        self.save_dir.mkdir(parents=True, exist_ok=True)


def sklearn_predict(model: RandomForestClassifier | RandomForestRegressor | MLPClassifier | MLPRegressor,
                    fingerprints: np.ndarray) -> np.ndarray:
    """Predicts properties using a scikit-learn model."""
    if isinstance(model, RandomForestClassifier) or isinstance(model, MLPClassifier):
        return model.predict_proba(fingerprints)[:, 1]
    elif isinstance(model, RandomForestRegressor) or isinstance(model, MLPRegressor):
        return model.predict(fingerprints)

    raise ValueError(f'Model type {type(model)} is not supported.')


def chemprop_predict(*args, **kwargs) -> np.ndarray:
    """Converts chemprop predictions from list of lists to a 1D numpy array (assumes single prediction task)."""
    return np.array(_chemprop_predict(*args, **kwargs))[:, 0]


def build_sklearn_model(dataset_type: str,
                        train_fingerprints: np.ndarray,
                        test_fingerprints: np.ndarray,
                        train_properties: list[int],
                        test_properties: list[int],
                        save_path: Path,
                        model_type: Literal['rf', 'mlp']) -> tuple[dict[str, float], np.ndarray]:
    """Trains, evaluates, and saves a scikit-learn model."""
    # Build model
    if model_type == 'rf':
        if dataset_type == 'classification':
            rf_class = RandomForestClassifier
        elif dataset_type == 'regression':
            rf_class = RandomForestRegressor
        else:
            raise ValueError(f'Dataset type "{dataset_type}" is not supported.')

        model = rf_class(n_jobs=-1, random_state=0)
    elif model_type == 'mlp':
        if dataset_type == 'classification':
            mlp_class = MLPClassifier
        elif dataset_type == 'regression':
            mlp_class = MLPRegressor
        else:
            raise ValueError(f'Dataset type "{dataset_type}" is not supported.')

        model = mlp_class(hidden_layer_sizes=(100, 100, 100), random_state=0)
    else:
        raise ValueError(f'Model type "{model_type}" is not supported.')

    print(model)

    # Train model
    model.fit(train_fingerprints, train_properties)

    # Save model
    with open(save_path, 'wb') as f:
        pickle.dump(model, f)

    # Make test predictions
    test_preds = sklearn_predict(model=model, fingerprints=test_fingerprints)

    # Make and evaluate test predictions
    if dataset_type == 'classification':
        scores = {
            'roc_auc': roc_auc_score(test_properties, test_preds),
            'prc_auc': average_precision_score(test_properties, test_preds)
        }
    elif dataset_type == 'regression':
        scores = {
            'mae': mean_absolute_error(test_properties, test_preds),
            'r2': r2_score(test_properties, test_preds),
        }
    else:
        raise ValueError(f'Dataset type "{dataset_type}" is not supported.')

    return scores, test_preds


def build_chemprop_data_loader(smiles: list[str],
                               fingerprints: Optional[np.ndarray],
                               properties: Optional[list[int]] = None,
                               shuffle: bool = False) -> MoleculeDataLoader:
    """Builds a chemprop MoleculeDataLoader."""
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


def build_chemprop_model(dataset_type: str,
                         train_smiles: list[str],
                         val_smiles: list[str],
                         test_smiles: list[str],
                         fingerprint_type: Optional[str],
                         train_fingerprints: Optional[np.ndarray],
                         val_fingerprints: Optional[np.ndarray],
                         test_fingerprints: Optional[np.ndarray],
                         property_name: str,
                         train_properties: list[int],
                         val_properties: list[int],
                         test_properties: list[int],
                         save_path: Path,
                         epochs: int) -> tuple[dict[str, float], np.ndarray]:
    """Trains, evaluates, and saves a chemprop model."""
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
    train_data_loader = build_chemprop_data_loader(
        smiles=train_smiles,
        fingerprints=train_fingerprints,
        properties=train_properties,
        shuffle=True
    )
    val_data_loader = build_chemprop_data_loader(
        smiles=val_smiles,
        fingerprints=val_fingerprints,
        properties=val_properties
    )
    test_data_loader = build_chemprop_data_loader(
        smiles=test_smiles,
        fingerprints=test_fingerprints,
        properties=test_properties
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
    val_metric = 'PRC-AUC' if dataset_type == 'classification' else 'MAE'
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

    # Make test predictions
    test_preds = chemprop_predict(model=model, data_loader=test_data_loader)

    # Evaluate predictions
    if dataset_type == 'classification':
        scores = {
            'roc_auc': roc_auc_score(test_properties, test_preds),
            'prc_auc': average_precision_score(test_properties, test_preds)
        }
    elif dataset_type == 'regression':
        scores = {
            'mae': mean_absolute_error(test_properties, test_preds),
            'r2': r2_score(test_properties, test_preds)
        }
    else:
        raise ValueError(f'Dataset type "{dataset_type}" is not supported.')

    return scores, test_preds


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
        start_time = time.time()

        if args.model_type == 'chemprop':
            scores, test_preds = build_chemprop_model(
                dataset_type=args.dataset_type,
                train_smiles=train_data[args.smiles_column],
                val_smiles=val_data[args.smiles_column],
                test_smiles=test_data[args.smiles_column],
                fingerprint_type=args.fingerprint_type,
                train_fingerprints=train_fingerprints,
                val_fingerprints=val_fingerprints,
                test_fingerprints=test_fingerprints,
                property_name=args.property_column,
                train_properties=train_data[args.property_column],
                val_properties=val_data[args.property_column],
                test_properties=test_data[args.property_column],
                save_path=args.save_dir / f'model_{model_num}.pt',
                epochs=args.epochs
            )
        else:
            scores, test_preds = build_sklearn_model(
                dataset_type=args.dataset_type,
                train_fingerprints=train_fingerprints,
                test_fingerprints=test_fingerprints,
                train_properties=train_data[args.property_column],
                test_properties=test_data[args.property_column],
                save_path=args.save_dir / f'model_{model_num}.pkl',
                model_type=args.model_type
            )

        scores['time_seconds'] = time.time() - start_time

        # Save test predictions
        test_df = pd.DataFrame({
            'smiles': test_data[args.smiles_column],
            args.property_column: test_data[args.property_column],
            'prediction': test_preds
        })
        test_df.to_csv(args.save_dir / f'model_{model_num}_test_preds.csv', index=False)

        # Print scores
        for score_name, score_value in scores.items():
            print(f'Test {score_name} = {score_value:.3f}')
        print()

        all_scores.append(scores)

    # Process and save scores
    score_names = list(all_scores[0])

    all_scores = pd.DataFrame(all_scores)
    all_scores['model'] = [f'model_{model_num}' for model_num in range(args.num_models)]
    all_scores = all_scores[['model'] + score_names]
    # TODO: just put mean and std in here
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
