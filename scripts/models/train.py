"""Trains a machine learning property prediction model."""
import time
from pathlib import Path
from random import Random

import numpy as np
import pandas as pd
from chemfunc import compute_fingerprints
from tqdm import trange

from chemprop_models import (
    chemprop_predict,
    chemprop_train
)
from evaluate import evaluate
from sklearn_models import sklearn_train
from synthemol.constants import DATASET_TYPES, FINGERPRINT_TYPES, MODEL_TYPES, SMILES_COL
from synthemol.models import sklearn_predict


def train(
        data_path: Path,
        save_dir: Path,
        dataset_type: DATASET_TYPES,
        model_type: MODEL_TYPES,
        property_column: str,
        smiles_column: str = SMILES_COL,
        fingerprint_type: FINGERPRINT_TYPES | None = None,
        num_models: int = 1,
        epochs: int = 30,
        num_workers: int = 0,
        use_gpu: bool = False
) -> None:
    """Trains a machine learning property prediction model.

    :param data_path: Path to CSV file containing data.
    :param save_dir: Path to a directory where the trained model(s) and results will be saved.
    :param dataset_type: Type of dataset.
    :param model_type: Type of model to train.
    :param property_column: The name of the column containing property values.
    :param smiles_column: The name of the column containing SMILES.
    :param fingerprint_type: Type of fingerprints to use as input features.
    :param num_models: The number of models to train using 10-fold cross-validation.
    :param epochs: The number of epochs to train the chemprop model.
    :param num_workers: Number of workers for the data loader (only applicable to chemprop model type).
    :param use_gpu: Whether to use GPU (only applicable to chemprop model type).
    """
    # Check compatibility of model and fingerprint type
    if model_type != 'chemprop' and fingerprint_type is None:
        raise ValueError('Must define fingerprint_type if using sklearn model.')

    # Load data
    data = pd.read_csv(data_path)
    print(f'Data size = {len(data):,}')

    # Compute fingerprints
    if fingerprint_type is not None:
        fingerprints = compute_fingerprints(data[smiles_column], fingerprint_type=fingerprint_type)
    else:
        fingerprints = None

    # Set up cross-validation
    num_folds = 10
    indices = np.tile(np.arange(num_folds), 1 + len(data) // num_folds)[:len(data)]
    random = Random(0)
    random.shuffle(indices)

    assert 1 <= num_models <= num_folds

    # Create save directory
    save_dir.mkdir(parents=True, exist_ok=True)

    # Run cross-validation
    all_scores = []
    for model_num in trange(num_models, desc='cross-val'):
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

        if model_type == 'chemprop':
            # Train and save chemprop model
            model = chemprop_train(
                dataset_type=dataset_type,
                train_smiles=train_data[smiles_column],
                val_smiles=val_data[smiles_column],
                fingerprint_type=fingerprint_type,
                train_fingerprints=train_fingerprints,
                val_fingerprints=val_fingerprints,
                property_name=property_column,
                train_properties=train_data[property_column],
                val_properties=val_data[property_column],
                epochs=epochs,
                save_path=save_dir / f'model_{model_num}.pt',
                num_workers=num_workers,
                use_gpu=use_gpu
            )

            print(model)

            # Make test predictions with chemprop model
            test_preds = chemprop_predict(
                model=model,
                smiles=test_data[smiles_column],
                fingerprints=test_fingerprints,
                num_workers=num_workers,
                use_gpu=use_gpu
            )
        else:
            # Train and save sklearn model
            model = sklearn_train(
                model_type=model_type,
                dataset_type=dataset_type,
                fingerprints=train_fingerprints,
                properties=train_data[property_column],
                save_path=save_dir / f'model_{model_num}.pkl'
            )

            print(model)

            # Make test predictions with sklearn model
            test_preds = sklearn_predict(
                model=model,
                fingerprints=test_fingerprints
            )

        # Evaluate test predictions
        scores = evaluate(
            true=test_data[property_column],
            preds=test_preds,
            dataset_type=dataset_type
        )

        # Record train/eval time
        scores['time_seconds'] = time.time() - start_time

        # Save test predictions
        test_df = pd.DataFrame({
            SMILES_COL: test_data[smiles_column],
            property_column: test_data[property_column],
            'prediction': test_preds
        })
        test_df.to_csv(save_dir / f'model_{model_num}_test_preds.csv', index=False)

        # Print scores
        for score_name, score_value in scores.items():
            print(f'Test {score_name} = {score_value:.3f}')
        print()

        all_scores.append(scores)

    # Process and save scores
    score_names = list(all_scores[0])

    all_scores = pd.DataFrame(all_scores)
    all_scores['model'] = [f'model_{model_num}' for model_num in range(num_models)]
    all_scores = all_scores[['model'] + score_names]
    all_scores.to_csv(save_dir / 'scores.csv', index=False)

    # Process and save summary scores
    summary_scores = {}
    for score_name in score_names:
        summary_scores[f'{score_name}_mean'] = np.mean(all_scores[score_name])
        summary_scores[f'{score_name}_std'] = np.std(all_scores[score_name])

    summary_scores = pd.DataFrame([summary_scores])
    summary_scores.to_csv(save_dir / 'summary_scores.csv', index=False)


if __name__ == '__main__':
    from tap import tapify

    tapify(train)
