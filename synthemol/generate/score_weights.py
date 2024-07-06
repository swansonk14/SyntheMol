"""Contains a class for holding and updating score weights for multiparameter models."""
from typing import Iterable, Literal


class ScoreWeights:
    def __init__(
        self,
        base_score_weights: Iterable[float],
        immutable: bool,
        score_names: Iterable[str] | None = None,
        score_signs: Iterable[Literal[1, -1]] | None = None,
    ) -> None:
        """Initialize the score weights object.

        :param base_score_weights: The initial score weights.
        :param immutable: Whether the score weights are immutable.
        :param score_names: The names of the scores. If None, defaults to "Score 1", "Score 2", etc.
        :param score_signs: The signs (+1 or -1) of the scores. If None, defaults to all +1.
        """
        self._base_score_weights = tuple(base_score_weights)
        self._immutable = immutable
        self._score_weights = self._base_score_weights
        self._score_signs = score_signs

        # Extract score names or set default score names
        if score_names is not None:
            self._score_names = tuple(score_names)

            if len(self._score_names) != len(self._base_score_weights):
                raise ValueError(
                    f"Number of score names ({len(self._score_names):,}) must match "
                    f"number of score weights ({len(self._base_score_weights):,})."
                )
        else:
            self._score_names = tuple(
                f"Score {i + 1}" for i in range(len(self._base_score_weights))
            )

        # Set default score signs
        if self._score_signs is None:
            self._score_signs = (1,) * self.num_weights

    @property
    def immutable(self) -> bool:
        """Returns whether the score weights are immutable."""
        return self._immutable

    @property
    def score_signs(self) -> tuple[int, ...]:
        """Returns the score signs."""
        return self._score_signs

    @property
    def signed_weights(self) -> tuple[float, ...]:
        """Returns the signed score weights."""
        return tuple(sign * weight for sign, weight in zip(self._score_signs, self._score_weights))

    @property
    def weights(self) -> tuple[float, ...]:
        """Returns the score weights."""
        return self._score_weights

    @weights.setter
    def weights(self, score_weights: Iterable[float]) -> None:
        """Sets the score weights.

        :param score_weights: The new score weights.
        """
        if self.immutable:
            raise ValueError("Score weights are immutable.")

        self._score_weights = tuple(score_weights)

    @property
    def num_weights(self) -> int:
        """Returns the number of score weights."""
        return len(self.weights)

    @property
    def score_names(self) -> tuple[str, ...]:
        """Returns the score names."""
        return self._score_names
