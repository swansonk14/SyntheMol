"""Contains a class for scoring molecules during generation."""
from abc import ABC, abstractmethod
from pathlib import Path

import torch
from chemfunc import compute_fingerprint
from rdkit import Chem
from rdkit.Chem.Crippen import MolLogP
from rdkit.Chem.QED import qed


from synthemol.constants import FINGERPRINT_TYPES, SCORE_TYPES
from synthemol.generate.score_weights import ScoreWeights
from synthemol.models import (
    chemprop_load,
    chemprop_load_scaler,
    chemprop_predict_on_molecule_ensemble,
    sklearn_load,
    sklearn_predict_on_molecule_ensemble,
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


class CLogPScorer(Scorer):
    """Scores molecules using CLogP."""

    def __call__(self, smiles: str) -> float:
        """Scores a molecule using CLogP.

        :param smiles: A SMILES string.
        :return: The score of the molecule.
        """
        return MolLogP(Chem.MolFromSmiles(smiles))


class SKLearnScorer(Scorer):
    """Scores molecules using a scikit-learn model or ensemble of models."""

    def __init__(self, model_path: Path, fingerprint_type: FINGERPRINT_TYPES) -> None:
        """Initialize the scorer.

        :param model_path: Path to a directory of model checkpoints (ensemble) or to a specific PT file.
        :param fingerprint_type: The type of fingerprint to use as input features for the model.
        """
        # Get model paths
        if model_path.is_dir():
            model_paths = list(model_path.glob("**/*.pkl"))

            if len(model_paths) == 0:
                raise ValueError(
                    f"Could not find any models in directory {model_path}."
                )
        else:
            model_paths = [model_path]

        # Save fingerprint type
        self.fingerprint_type = fingerprint_type

        # Load scikit-learn models
        self.models = [
            sklearn_load(model_path=model_path) for model_path in model_paths
        ]

    def __call__(self, smiles: str) -> float:
        """Scores a molecule using a scikit-learn model or ensemble of models.

        :param smiles: A SMILES string.
        :return: The score of the molecule.
        """
        # Compute fingerprint
        fingerprint = compute_fingerprint(
            smiles, fingerprint_type=self.fingerprint_type
        )

        # Make prediction
        return sklearn_predict_on_molecule_ensemble(
            models=self.models, fingerprint=fingerprint
        )


class ChempropScorer(Scorer):
    """Scores molecules using a Chemprop model or ensemble of models."""

    def __init__(
        self,
        model_path: Path,
        fingerprint_type: FINGERPRINT_TYPES | None = None,
        device: torch.device = torch.device("cpu"),
    ) -> None:
        """Initialize the scorer.

        :param model_path: Path to a directory of model checkpoints (ensemble) or to a specific PT file.
        :param fingerprint_type: The type of fingerprint to use as input features for the model.
        :param device: The device on which to run the model.
        """
        # Get model paths
        if model_path.is_dir():
            model_paths = list(model_path.glob("**/*.pt"))

            if len(model_paths) == 0:
                raise ValueError(
                    f"Could not find any models in directory {model_path}."
                )
        else:
            model_paths = [model_path]

        # Save fingerprint type
        self.fingerprint_type = fingerprint_type

        # Ensure reproducibility
        torch.manual_seed(0)

        if device.type == "cpu":
            torch.use_deterministic_algorithms(True)

        # Load models
        self.models = [
            chemprop_load(model_path=model_path, device=device)
            for model_path in model_paths
        ]
        self.scalers = [
            chemprop_load_scaler(model_path=model_path) for model_path in model_paths
        ]

    def __call__(self, smiles: str) -> float:
        """Scores a molecule using a Chemprop model or ensemble of models.

        :param smiles: A SMILES string.
        :return: The score of the molecule.
        """
        # Compute fingerprint
        if self.fingerprint_type is not None:
            fingerprint = compute_fingerprint(
                smiles, fingerprint_type=self.fingerprint_type
            )
        else:
            fingerprint = None

        # Make prediction
        return chemprop_predict_on_molecule_ensemble(
            models=self.models,
            smiles=smiles,
            fingerprint=fingerprint,
            scalers=self.scalers,
        )


def create_scorer(
    score_type: SCORE_TYPES,
    model_path: Path | None = None,
    fingerprint_type: FINGERPRINT_TYPES | None = None,
    device: torch.device = torch.device("cpu"),
) -> Scorer:
    """Creates a scorer object that scores a molecule.

    :param score_type: The type of score to use.
    :param model_path: For score types that are model-based ("random_forest" and "chemprop"), the corresponding
            model path should be a path to a directory of model checkpoints (ensemble)
            or to a specific PKL or PT file containing a trained model with a single output.
            For score types that are not model-based, the corresponding model path must be None.
    :param fingerprint_type: For score types that are model-based and require fingerprints as input, the corresponding
            fingerprint type should be the type of fingerprint (e.g., "rdkit").
            For model-based scores that don't require fingerprints or non-model-based scores,
            the corresponding fingerprint type must be None.
    :param device: The device on which to run the scorer.
    """
    if score_type == "qed":
        if model_path is not None:
            raise ValueError("QED does not use a model path.")

        if fingerprint_type is not None:
            raise ValueError("QED does not use fingerprints.")

        scorer = QEDScorer()
    if score_type == "clogp":
        if model_path is not None:
            raise ValueError("CLogP does not use a model path.")

        if fingerprint_type is not None:
            raise ValueError("CLogP does not use fingerprints.")

        scorer = CLogPScorer()
    elif score_type == "chemprop":
        if model_path is None:
            raise ValueError("Chemprop requires a model path.")

        scorer = ChempropScorer(
            model_path=model_path, fingerprint_type=fingerprint_type, device=device
        )
    elif score_type == "random_forest":
        if model_path is None:
            raise ValueError("Random forest requires a model path.")

        if fingerprint_type is None:
            raise ValueError("Random forest requires a fingerprint type.")

        scorer = SKLearnScorer(model_path=model_path, fingerprint_type=fingerprint_type)
    else:
        raise ValueError(f"Score type {score_type} is not supported.")

    return scorer


class MoleculeScorer:
    def __init__(
        self,
        score_types: list[SCORE_TYPES],
        score_weights: ScoreWeights,
        model_paths: list[Path | None] | None = None,
        fingerprint_types: list[FINGERPRINT_TYPES | None] | None = None,
        device: torch.device = torch.device("cpu"),
        smiles_to_scores: dict[str, list[float]] | None = None,
    ) -> None:
        """Initialize the MoleculeScorer, which contains a collection of one or more individual scorers.

        :param score_types: List of types of scores to score molecules.
        :param score_weights: Weights for each scorer for defining the reward function.
        :param model_paths: For score types that are model-based ("random_forest" and "chemprop"), the corresponding
            model path should be a path to a directory of model checkpoints (ensemble)
            or to a specific PKL or PT file containing a trained model with a single output.
            For score types that are not model-based, the corresponding model path must be None.
            If all score types are not model-based, this argument can be None.
        :param fingerprint_types: For score types that are model-based and require fingerprints as input, the corresponding
            fingerprint type should be the type of fingerprint (e.g., "rdkit").
            For model-based scores that don't require fingerprints or non-model-based scores,
            the corresponding fingerprint type must be None.
            If all score types do not require fingerprints, this argument can be None.
        :param device: The device on which to run the scorer.
        :param smiles_to_scores: An optional dictionary mapping SMILES to precomputed scores.
        """
        # Save parameters
        self.score_weights = score_weights
        self.smiles_to_individual_scores = smiles_to_scores

        # Handle None model_paths and fingerprint_types
        if model_paths is None:
            model_paths = [None] * len(score_types)

        if fingerprint_types is None:
            fingerprint_types = [None] * len(score_types)

        # Create individual scorers
        self.scorers = [
            create_scorer(
                score_type=score_type,
                model_path=model_path,
                fingerprint_type=fingerprint_type,
                device=device,
            )
            for score_type, model_path, fingerprint_type in zip(
                score_types, model_paths, fingerprint_types
            )
        ]

        # Initialize a cache for molecule scores
        self.smiles_to_score: dict[str, float] = {}

    @property
    def num_scores(self) -> int:
        """Returns the number of scores."""
        return self.score_weights.num_weights

    def compute_individual_scores(self, smiles: str) -> list[float]:
        """Computes the individual scores of a molecule (with caching).

        :param smiles: A SMILES string.
        :return: A list of individual scores.
        """
        if smiles in self.smiles_to_individual_scores:
            individual_scores = self.smiles_to_individual_scores[smiles]
        else:
            individual_scores = [scorer(smiles) for scorer in self.scorers]
            self.smiles_to_individual_scores[smiles] = individual_scores

        return individual_scores

    def __call__(self, smiles: str) -> float:
        """Scores a molecule using a collection of scorers.

        :param smiles: A SMILES string.
        :return: The score of the molecule.
        """
        # Check if score is cached
        if smiles in self.smiles_to_score:
            return self.smiles_to_score[smiles]

        # Look up individual scores or compute the scores and cache them
        individual_scores = self.compute_individual_scores(smiles=smiles)

        # Compute weighted average score
        score = sum(
            individual_score * weight
            for individual_score, weight in zip(
                individual_scores, self.score_weights.weights
            )
        )

        # Cache score
        self.smiles_to_score[smiles] = score

        return score
