"""Utility functions for generating molecules."""
from functools import cache
from pathlib import Path
from typing import Callable

import numpy as np
import pandas as pd
import torch
from chemfunc import compute_fingerprint
from rdkit import Chem
from rdkit.Chem.QED import qed

from synthemol.constants import (
    FINGERPRINT_TYPES,
    MODEL_TYPES,
    ROLLOUT_COL,
    SCORE_COL,
    SMILES_COL
)
from synthemol.generate.node import Node
from synthemol.models import (
    chemprop_load,
    chemprop_load_scaler,
    chemprop_predict_on_molecule_ensemble,
    sklearn_load,
    sklearn_predict_on_molecule_ensemble
)


def create_model_scorer(
        model_type: MODEL_TYPES,
        model_path: Path,
        fingerprint_type: FINGERPRINT_TYPES,
        device: torch.device = torch.device('cpu')
) -> Callable[[str], float]:
    """Creates a function that scores a molecule using a model or ensemble of models.

    :param model_type: The type of model to use.
    :param model_path: Path to a directory of model checkpoints (ensemble) or to a specific PKL or PT file.
    :param fingerprint_type: The type of fingerprint to use as input features for the model.
    :param device: The device on which to run the model.
    """
    # Check compatibility of model and fingerprint type
    if model_type in {'random_forest', 'mlp'} and fingerprint_type == 'none':
        raise ValueError('Must define fingerprint_type if using a scikit-learn model.')

    # Get model paths
    if model_type != 'qed' and model_path.is_dir():
        model_paths = list(model_path.glob('**/*.pt' if model_type == 'chemprop' else '**/*.pkl'))

        if len(model_paths) == 0:
            raise ValueError(f'Could not find any models in directory {model_path}.')
    else:
        model_paths = [model_path]

    # Load models and set up scoring function
    if model_type == 'qed':
        def model_scorer(smiles: str) -> float:
            return qed(Chem.MolFromSmiles(smiles))
    elif model_type == 'chemprop':
        # Ensure reproducibility
        torch.manual_seed(0)

        if device.type == 'cpu':
            torch.use_deterministic_algorithms(True)

        # Load models
        models = [chemprop_load(model_path=model_path, device=device) for model_path in model_paths]
        scalers = [chemprop_load_scaler(model_path=model_path) for model_path in model_paths]

        # Set up model scoring function for ensemble of chemprop models
        def model_scorer(smiles: str) -> float:
            # Compute fingerprint
            if fingerprint_type != 'none':
                fingerprint = compute_fingerprint(smiles, fingerprint_type=fingerprint_type)
            else:
                fingerprint = None

            # Make prediction
            return chemprop_predict_on_molecule_ensemble(
                models=models,
                smiles=smiles,
                fingerprint=fingerprint,
                scalers=scalers
            )
    else:
        # Load scikit-learn models
        models = [sklearn_load(model_path=model_path) for model_path in model_paths]

        # Set up model scoring function for ensemble of scikit-learn models
        def model_scorer(smiles: str) -> float:
            # Compute fingerprint
            fingerprint = compute_fingerprint(smiles, fingerprint_type=fingerprint_type)

            # Make prediction
            return sklearn_predict_on_molecule_ensemble(
                models=models,
                fingerprint=fingerprint
            )

    return model_scorer


def create_model_scoring_fn(
        model_types: list[MODEL_TYPES],
        model_paths: list[Path],
        fingerprint_types: list[FINGERPRINT_TYPES],
        model_weights: tuple[float] = (1.0,),
        device: torch.device = torch.device('cpu'),
        smiles_to_score: dict[str, float] | None = None
) -> Callable[[str], float]:
    """Creates a function that scores a molecule using a weight combination of models or ensembles of models.

    :param model_types: List of types of models provided by model_paths.
    :param model_paths: List of paths with each path pointing to a directory of model checkpoints (ensemble)
                        or to a specific PKL or PT file containing a trained model.
                        Note: All models must have a single output.
    :param fingerprint_types: List of types of fingerprints to use as input features for the model_paths.
    :param model_weights: List of weights for each model/ensemble in model_paths for defining the reward function.
    :param device: The device on which to run the model.
    :param smiles_to_score: An optional dictionary mapping SMILES to precomputed scores.
    :return: A function that scores a molecule using a weight combination of models or ensembles of models.
    """
    # Set up model scorer functions
    model_scorers = [
        create_model_scorer(
            model_type=model_type,
            model_path=model_path,
            fingerprint_type=fingerprint_type,
            device=device
        )
        for model_type, model_path, fingerprint_type in zip(model_types, model_paths, fingerprint_types)
    ]

    # Build model scoring function including precomputed building block scores
    @cache
    def model_scoring_fn(smiles: str) -> float:
        # If SMILES is in precomputed scores, return the score
        if smiles_to_score is not None and smiles in smiles_to_score:
            return smiles_to_score[smiles]

        # Otherwise, compute the score using a weighted combination of the model scorers
        return sum(
            model_weight * model_scorer(smiles=smiles)
            for model_weight, model_scorer in zip(model_weights, model_scorers)
        )

    return model_scoring_fn


def save_generated_molecules(
        nodes: list[Node],
        building_block_id_to_smiles: dict[int, str],
        save_path: Path
) -> None:
    """Save generated molecules to a CSV file.

    :param nodes: A list of Nodes containing molecules. Only nodes with a single molecule are saved.
    :param building_block_id_to_smiles: A dictionary mapping building block IDs to SMILES.
    :param save_path: A path to a CSV file where the molecules will be saved.
    """
    # Convert construction logs from lists to dictionaries
    construction_dicts = []
    max_reaction_num = 0
    reaction_num_to_max_reactant_num = {}

    for node in nodes:
        construction_dict = {'num_reactions': len(node.construction_log)}
        max_reaction_num = max(max_reaction_num, len(node.construction_log))

        for reaction_index, reaction_log in enumerate(node.construction_log):
            reaction_num = reaction_index + 1
            construction_dict[f'reaction_{reaction_num}_id'] = reaction_log.reaction_id

            reaction_num_to_max_reactant_num[reaction_num] = max(
                reaction_num_to_max_reactant_num.get(reaction_num, 0),
                len(reaction_log.reactant_ids)
            )

            for reactant_index, reactant_id in enumerate(reaction_log.reactant_ids):
                reactant_num = reactant_index + 1
                construction_dict[f'building_block_{reaction_num}_{reactant_num}_id'] = reactant_id
                construction_dict[f'building_block_{reaction_num}_{reactant_num}_smiles'] = building_block_id_to_smiles.get(reactant_id, '')

        construction_dicts.append(construction_dict)

    # Specify column order for CSV file
    columns = [SMILES_COL, 'node_id', 'num_expansions', ROLLOUT_COL, SCORE_COL, 'Q_value', 'num_reactions']

    for reaction_num in range(1, max_reaction_num + 1):
        columns.append(f'reaction_{reaction_num}_id')

        for reactant_num in range(1, reaction_num_to_max_reactant_num[reaction_num] + 1):
            columns.append(f'building_block_{reaction_num}_{reactant_num}_id')
            columns.append(f'building_block_{reaction_num}_{reactant_num}_smiles')

    # Save data
    save_path.parent.mkdir(parents=True, exist_ok=True)
    data = pd.DataFrame(
        data=[
            {
                SMILES_COL: '.'.join(node.molecules),
                'node_id': node.node_id,
                'num_expansions': node.N,
                ROLLOUT_COL: node.rollout_num,
                SCORE_COL: node.P,
                'Q_value': node.Q(),
                **construction_dict
            }
            for node, construction_dict in zip(nodes, construction_dicts)
        ],
        columns=columns
    )
    data.to_csv(save_path, index=False)
