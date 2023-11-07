"""Contains a class for scoring molecules during generation."""
from abc import ABC, abstractmethod
from pathlib import Path

import torch
from chemfunc import compute_fingerprint
from rdkit import Chem
from rdkit.Chem.QED import qed

from synthemol.constants import FINGERPRINT_TYPES, MODEL_TYPES
from synthemol.models import (
    chemprop_load,
    chemprop_load_scaler,
    chemprop_predict_on_molecule_ensemble,
    sklearn_load,
    sklearn_predict_on_molecule_ensemble
)


class Scorer(ABC):
    """Base class for scoring molecules."""

    @abstractmethod
    def __call__(self, smiles: str) -> float:
        """Scores a molecule.

        :param smiles: A SMILES string.
        :return: The score of the molecule.
        """
        pass


class QEDScorer(Scorer):
    """Scores molecules using QED."""

    def __call__(self, smiles: str) -> float:
        """Scores a molecule using QED.

        :param smiles: A SMILES string.
        :return: The score of the molecule.
        """
        return qed(Chem.MolFromSmiles(smiles))


class SKLearnScorer(Scorer):
    """Scores molecules using a scikit-learn model or ensemble of models."""

    def __init__(
            self,
            model_path: Path,
            fingerprint_type: FINGERPRINT_TYPES
    ) -> None:
        """Initialize the scorer.

        :param model_path: Path to a directory of model checkpoints (ensemble) or to a specific PT file.
        :param fingerprint_type: The type of fingerprint to use as input features for the model.
        """
        # Get model paths
        if model_path.is_dir():
            model_paths = list(model_path.glob('**/*.pkl'))

            if len(model_paths) == 0:
                raise ValueError(f'Could not find any models in directory {model_path}.')
        else:
            model_paths = [model_path]

        # Save fingerprint type
        self.fingerprint_type = fingerprint_type

        # Load scikit-learn models
        self.models = [sklearn_load(model_path=model_path) for model_path in model_paths]

    def __call__(self, smiles: str) -> float:
        """Scores a molecule using a scikit-learn model or ensemble of models.

        :param smiles: A SMILES string.
        :return: The score of the molecule.
        """
        # Compute fingerprint
        fingerprint = compute_fingerprint(smiles, fingerprint_type=self.fingerprint_type)

        # Make prediction
        return sklearn_predict_on_molecule_ensemble(
            models=self.models,
            fingerprint=fingerprint
        )


class ChempropScorer(Scorer):
    """Scores molecules using a Chemprop model or ensemble of models."""

    def __init__(
            self,
            model_path: Path,
            fingerprint_type: FINGERPRINT_TYPES,
            device: torch.device = torch.device('cpu')
    ) -> None:
        """Initialize the scorer.

        :param model_path: Path to a directory of model checkpoints (ensemble) or to a specific PT file.
        :param fingerprint_type: The type of fingerprint to use as input features for the model.
        :param device: The device on which to run the model.
        """
        # Get model paths
        if model_path.is_dir():
            model_paths = list(model_path.glob('**/*.pt'))

            if len(model_paths) == 0:
                raise ValueError(f'Could not find any models in directory {model_path}.')
        else:
            model_paths = [model_path]

        # Save fingerprint type
        self.fingerprint_type = fingerprint_type

        # Ensure reproducibility
        torch.manual_seed(0)

        if device.type == 'cpu':
            torch.use_deterministic_algorithms(True)

        # Load models
        self.models = [chemprop_load(model_path=model_path, device=device) for model_path in model_paths]
        self.scalers = [chemprop_load_scaler(model_path=model_path) for model_path in model_paths]

    def __call__(self, smiles: str) -> float:
        """Scores a molecule using a Chemprop model or ensemble of models.

        :param smiles: A SMILES string.
        :return: The score of the molecule.
        """
        # Compute fingerprint
        if self.fingerprint_type != 'none':
            fingerprint = compute_fingerprint(smiles, fingerprint_type=self.fingerprint_type)
        else:
            fingerprint = None

        # Make prediction
        return chemprop_predict_on_molecule_ensemble(
            models=self.models,
            smiles=smiles,
            fingerprint=fingerprint,
            scalers=self.scalers
        )


# TODO: no model_path or fingerprint_type needed for qed
def create_model_scorer(
        model_type: MODEL_TYPES,
        model_path: Path,
        fingerprint_type: FINGERPRINT_TYPES,
        device: torch.device = torch.device('cpu')
) -> Scorer:
    """Creates a scorer object that scores a molecule.

    :param model_type: The type of model to use.
    :param model_path: Path to a directory of model checkpoints (ensemble) or to a specific PKL or PT file.
    :param fingerprint_type: The type of fingerprint to use as input features for the model.
    :param device: The device on which to run the model.
    """
    # Check compatibility of model and fingerprint type
    if model_type in {'random_forest', 'mlp'} and fingerprint_type == 'none':
        raise ValueError('Must define fingerprint_type if using a scikit-learn model.')

    # Load models and set up scoring function
    if model_type == 'qed':
        scorer = QEDScorer()
    elif model_type == 'chemprop':
        scorer = ChempropScorer(
            model_path=model_path,
            fingerprint_type=fingerprint_type,
            device=device
        )
    else:
        scorer = SKLearnScorer(
            model_path=model_path,
            fingerprint_type=fingerprint_type
        )

    return scorer


class MoleculeScorer:
    def __init__(
            self,
            model_types: list[MODEL_TYPES],
            model_paths: list[Path],
            fingerprint_types: list[FINGERPRINT_TYPES],
            model_weights: tuple[float, ...] = (1.0,),
            device: torch.device = torch.device('cpu'),
            smiles_to_scores: dict[str, list[float]] | None = None
    ) -> None:
        """Initialize the MoleculeScorer, which contains a collection of one or more individual scorers.

        :param model_types: List of types of models provided by model_paths.
        :param model_paths: List of paths with each path pointing to a directory of model checkpoints (ensemble)
                            or to a specific PKL or PT file containing a trained model.
                            Note: All models must have a single output.
        :param fingerprint_types: List of types of fingerprints to use as input features for the model_paths.
        :param model_weights: Weights for each model/ensemble in model_paths for defining the reward function.
        :param device: The device on which to run the model.
        :param smiles_to_scores: An optional dictionary mapping SMILES to precomputed scores.
        """
        # Create individual scorers
        self.scorers = [
            create_model_scorer(
                model_type=model_type,
                model_path=model_path,
                fingerprint_type=fingerprint_type,
                device=device
            )
            for model_type, model_path, fingerprint_type in zip(model_types, model_paths, fingerprint_types)
        ]

        # Save model weights and smiles to scores
        self._model_weights = model_weights
        self.smiles_to_individual_scores = smiles_to_scores

        # Initialize a cache for molecule scores
        self.smiles_to_score: dict[str, float] = {}

    @property
    def model_weights(self) -> tuple[float, ...]:
        """Returns the model weights."""
        return self._model_weights

    @model_weights.setter
    def model_weights(self, model_weights: tuple[float, ...]) -> None:
        """Sets the model weights and clears the score cache."""
        self._model_weights = model_weights

        # Clear cache
        self.smiles_to_score = {}

    def __call__(self, smiles: str) -> float:
        """Scores a molecule using a collection of scorers.

        :param smiles: A SMILES string.
        :return: The score of the molecule.
        """
        # Check if score is cached
        if smiles in self.smiles_to_score:
            return self.smiles_to_score[smiles]

        # Look up individual scores
        if smiles in self.smiles_to_individual_scores:
            individual_scores = self.smiles_to_individual_scores[smiles]
        else:
            individual_scores = [scorer(smiles) for scorer in self.scorers]

        # Compute weighted average score
        score = sum(
            individual_score * weight
            for individual_score, weight in zip(individual_scores, self.model_weights)
        )

        # Cache score
        self.smiles_to_score[smiles] = score

        return score
