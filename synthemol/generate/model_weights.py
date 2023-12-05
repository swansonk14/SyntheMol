"""Contains a class for holding and updating model weights for multiparameter models."""
from typing import Iterable


class ModelWeights:
    def __init__(
        self,
        base_model_weights: Iterable[float],
        immutable: bool,
        model_names: Iterable[str] | None = None,
    ) -> None:
        """Initialize the model weights object.

        :param base_model_weights: The initial model weights.
        :param immutable: Whether the model weights are immutable.
        :param model_names: The names of the models. If None, defaults to "Model 1", "Model 2", etc.
        """
        self._base_model_weights = tuple(base_model_weights)
        self._immutable = immutable
        self._model_weights = self._base_model_weights

        # Extract model names or set default model names
        if model_names is not None:
            self._model_names = tuple(model_names)

            if len(self._model_names) != len(self._base_model_weights):
                raise ValueError(
                    f"Number of model names ({len(self._model_names):,}) must match "
                    f"number of model weights ({len(self._base_model_weights):,})."
                )
        else:
            self._model_names = tuple(
                f"Model {i + 1}" for i in range(len(self._base_model_weights))
            )

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

    @property
    def model_names(self) -> tuple[str, ...]:
        """Returns the model names."""
        return self._model_names
