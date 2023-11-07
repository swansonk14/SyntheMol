"""Contains a class for a multilayer perceptron (MLP)."""
import torch
import torch.nn as nn


class MLP(nn.Module):
    """A multilayer perceptron model."""

    def __init__(
        self,
        input_dim: int,
        hidden_dim: int,
        output_dim: int,
        num_layers: int,
        sigmoid: bool,
        device: torch.device = torch.device('cpu'),
    ) -> None:
        """Initialize the model.

        :param input_dim: The dimensionality of the input to the model.
        :param hidden_dim: The dimensionality of the hidden layers.
        :param output_dim: The dimensionality of the output of the model.
        :param num_layers: The number of layers.
        :param sigmoid: Whether to apply a sigmoid function to the output (during inference only).
        :param device: The device to use for the model.
        """
        super(MLP, self).__init__()

        assert num_layers > 1

        # Create layer dimensions
        layer_dims = [input_dim] + [hidden_dim] * (num_layers - 1) + [output_dim]

        # Create layers
        self.layers = nn.ModuleList(
            [
                nn.Linear(layer_dims[i], layer_dims[i + 1])
                for i in range(len(layer_dims) - 1)
            ]
        )

        self.sigmoid = sigmoid

        # Create activation function
        self.activation = nn.ReLU()

        # Set device
        self.device = device

    def forward(self, X: torch.Tensor) -> torch.Tensor:
        """Runs the model on the data.

        :param X: A tensor containing the input (batch_size, input_dim).
        :return: A tensor containing the model's prediction (batch_size, output_dim).
        """
        # Move data to device
        X = X.to(self.device)

        # Apply layers
        for i, layer in enumerate(self.layers):
            X = layer(X)

            if i < len(self.layers) - 1:
                X = self.activation(X)

        # Apply sigmoid during inference if desired
        if self.sigmoid and not self.training:
            X = torch.sigmoid(X)

        return X
