"""Map building blocks to model prediction scores."""
import json
from pathlib import Path

import pandas as pd
from tap import tapify

from SyntheMol.constants import FINGERPRINT_TYPES, MODEL_TYPE, SMILES_COL
from SyntheMol.models.predict_model import predict_ensemble  # TODO: fix this import


def map_building_block_to_scores(
        building_blocks_path: Path,
        model_path: Path,
        save_path: Path,
        model_type: MODEL_TYPE,
        fingerprint_type: FINGERPRINT_TYPES = None,
        smiles_column: str = SMILES_COL
) -> None:
    """Map building blocks to prediction scores.

    :param building_blocks_path: Path to a CSV file containing building blocks.
    :param model_path: Path to a directory of model checkpoints or to a specific PKL or PT file containing a trained model.
    :param save_path: Path to a JSON file where a dictionary mapping building blocks to model scores will be saved.
    :param model_type: Type of model to train. 'rf' = random forest. 'mlp' = multilayer perceptron.
    :param fingerprint_type: Type of fingerprints to use as input features.
    :param smiles_column: Name of the column containing SMILES.
    """
    # Load building blocks
    building_blocks = sorted(set(pd.read_csv(building_blocks_path)[smiles_column]))

    # Make predictions
    scores = predict_ensemble(
        smiles=building_blocks,
        fingerprint_type=fingerprint_type,
        model_type=model_type,
        model_path=model_path
    )

    # Map building blocks to predictions
    building_block_to_model_score = dict(zip(building_blocks, scores))

    # Save building block to model score mapping
    save_path.parent.mkdir(parents=True, exist_ok=True)

    with open(save_path, 'w') as f:
        json.dump(building_block_to_model_score, f, indent=4, sort_keys=True)


if __name__ == '__main__':
    tapify(map_building_block_to_scores)
