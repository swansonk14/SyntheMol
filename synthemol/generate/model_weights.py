"""Contains a class for holding and updating model weights for multiparameter models."""
from typing import Iterable


class ModelWeights:
    def __init__(self, base_model_weights: Iterable[float], immutable: bool) -> None:
        """Initialize the model weights object.

        :param base_model_weights: The initial model weights.
        :param immutable: Whether the model weights are immutable.
        """
        self._base_model_weights = tuple(base_model_weights)
        self._immutable = immutable
        self._model_weights = self._base_model_weights

    @property
    def immutable(self) -> bool:
        """Returns whether the model weights are immutable."""
        return self._immutable

    @property
    def weights(self) -> tuple[float, ...]:
        """Returns the model weights."""
        return self._model_weights

    @weights.setter
    def weights(self, model_weights: Iterable[float]) -> None:
        """Sets the model weights.

        :param model_weights: The new model weights.
        """
        if self.immutable:
            raise ValueError("Model weights are immutable.")

        self._model_weights = tuple(model_weights)

    @property
    def num_weights(self) -> int:
        """Returns the number of model weights."""
        return len(self.weights)
