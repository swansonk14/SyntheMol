"""Contains reinforcement learning models for use in generating molecules."""
import torch
import torch.nn as nn
from tqdm import trange

from chemfunc.molecular_fingerprints import compute_fingerprints


class MLP(nn.Module):
    """A multilayer perceptron model."""

    def __init__(
        self,
        input_dim: int,
        hidden_dim: int,
        output_dim: int,
        num_layers: int
    ) -> None:
        """Initialize the model.

        :param input_dim: The dimensionality of the input to the model.
        :param hidden_dim: The dimensionality of the hidden layers.
        :param output_dim: The dimensionality of the output of the model.
        :param num_layers: The number of layers.
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

        # Create activation function
        self.activation = nn.ReLU()

    def forward(self, X: torch.Tensor) -> torch.Tensor:
        """Runs the model on the data.

        :param X: A tensor containing the input (batch_size, input_dim).
        :return: A tensor containing the model's prediction (batch_size, output_dim).
        """
        # Apply layers
        for i, layer in enumerate(self.layers):
            X = layer(X)

            if i < len(self.layers) - 1:
                X = self.activation(X)

        return X


class RLModel:
    """A reinforcement learning model for predicting the reward of a molecule or set of molecules."""

    def __init__(
            self,
            max_num_molecules: int = 3,
            features_size: int = 200,
            batch_size: int = 32,
            num_epochs: int = 5,
            hidden_dim: int = 100,
            num_layers: int = 2,
            learning_rate: float = 1e-3,
            num_workers: int = 0
    ) -> None:
        """Initializes the model.

        :param max_num_molecules: The maximum number of molecules to process at a time.
        :param features_size: The size of the features for each molecule.
        :param batch_size: The batch size.
        :param num_epochs: The number of epochs to train for.
        :param hidden_dim: The dimensionality of the hidden layers.
        :param num_layers: The number of layers.
        :param learning_rate: The learning rate.
        :param num_workers: The number of workers to use for data loading.
        """
        self.max_num_molecules = max_num_molecules
        self.features_size = features_size
        self.total_features_size = max_num_molecules * features_size
        self.batch_size = batch_size
        self.num_epochs = num_epochs
        self.num_workers = num_workers

        self.molecule_tuples: list[tuple[str]] = []
        self.rewards: list[float] = []
        self.features: list[torch.Tensor] = []

        self.smiles_to_features = {}

        self.model = MLP(
            input_dim=self.total_features_size,
            hidden_dim=hidden_dim,
            output_dim=1,
            num_layers=num_layers
        )
        self.optimizer = torch.optim.Adam(self.model.parameters(), lr=learning_rate)
        self.loss_fn = nn.HuberLoss()

    def compute_rdkit_features(self, molecule_tuples: list[tuple[str]]) -> torch.Tensor:
        """Computes the RDKit features for each molecule in each tuple of molecules.

        :param molecule_tuples: A list of tuples of SMILES strings representing one or more molecules.
        :return: A tensor containing the features for each molecule in each tuple of molecules (num_tuples, total_features_size).
        """
        # Get all molecules
        molecule_nums = [len(molecules) for molecules in molecule_tuples]
        all_molecules = [molecule for molecules in molecule_tuples for molecule in molecules]

        # Determine unknown molecules
        unknown_molecules = [molecule for molecule in all_molecules if molecule not in self.smiles_to_features]

        # Compute features for unknown molecules and add to dictionary
        if unknown_molecules:
            unknown_features = compute_fingerprints(mols=unknown_molecules, fingerprint_type='rdkit')
            unknown_features = torch.from_numpy(unknown_features)

            for i, molecule in enumerate(unknown_molecules):
                self.smiles_to_features[molecule] = unknown_features[i]

        # Get all molecules and their features
        all_features = torch.stack([self.smiles_to_features[molecule] for molecule in all_molecules])

        # Set up feature vectors for each combination of molecules
        features = torch.zeros((len(molecule_tuples), self.total_features_size))
        index = 0
        for i, molecule_num in enumerate(molecule_nums):
            molecules_features = all_features[index:index + molecule_num].flatten()
            features[i, :len(molecules_features)] = molecules_features
            index += molecule_num

        return features

    def buffer(self, molecules: tuple[str], reward: float) -> None:
        """Adds a training example to the buffer.

        :param molecules: A tuple of SMILES strings representing one or more molecules.
        :param reward: The reward of those molecules (i.e., the score of the final molecule constructed
                       from these molecules and potentially others).
        """
        # Add molecules and reward
        self.molecule_tuples.append(molecules)
        self.rewards.append(reward)

        # Compute features
        self.features.append(self.compute_rdkit_features([molecules])[0])

    def train(self) -> None:
        """Trains the model on the examples in the buffer."""
        # Set model to train mode
        self.model.train()

        # TODO: upsample new examples and downsample old examples

        # Create dataset
        dataset = torch.utils.data.TensorDataset(torch.stack(self.features), torch.tensor(self.rewards))

        # Create dataloader
        dataloader = torch.utils.data.DataLoader(
            dataset=dataset, batch_size=self.batch_size, shuffle=True, num_workers=self.num_workers
        )

        # Loop over epochs
        for _ in trange(self.num_epochs, desc='Training RL model'):
            # Loop over batches of molecule features and rewards
            for batch_features, batch_rewards in dataloader:
                # Make predictions
                predictions = self.model(batch_features).squeeze(dim=-1)

                # Compute loss
                loss = self.loss_fn(predictions, batch_rewards)

                # Backpropagate
                self.model.zero_grad()
                loss.backward()
                self.optimizer.step()

    def evaluate(self) -> dict[str, float]:
        """Evaluates the model on the examples in the buffer.

        :return: A dictionary of metrics.
        """
        # Make predictions on all examples in the buffer
        predictions = self.predict(features=torch.stack(self.features))

        # Convert rewards to tensor
        rewards = torch.tensor(self.rewards)

        # Evaluate predictions
        results = {
            'RL Mean Squared Error': torch.nn.functional.mse_loss(predictions, rewards),
            'RL Loss': self.loss_fn(predictions, rewards).item()
        }

        return results

    def predict(
            self,
            molecule_tuples: list[tuple[str]] | None = None,
            features: torch.Tensor | None = None
    ) -> torch.Tensor:
        """Predicts the reward for each tuple of molecules.

        Must provide one of `molecule_tuples` or `features` but not both.

        :param molecule_tuples: A list of tuples of SMILES strings representing one or more molecules.
        :param features: A tensor containing the features for each molecule in each tuple of molecules
                         (num_tuples, total_features_size).
        :return: A 1D tensor of rewards for each tuple of molecules.
        """
        # Check inputs
        if (molecule_tuples is None) == (features is None):
            raise ValueError('Must provide either `molecule_tuples` or `features` but not both.')

        # Set model to eval mode
        self.model.eval()

        # Compute features in bulk across all molecules
        if features is None:
            features = self.compute_rdkit_features(molecule_tuples)

        # Create dataset
        dataset = torch.utils.data.TensorDataset(features)

        # Create dataloader
        dataloader = torch.utils.data.DataLoader(
            dataset=dataset, batch_size=self.batch_size, shuffle=False, num_workers=self.num_workers
        )

        # Loop over batches of molecules and make reward predictions
        rewards = []

        with torch.no_grad():
            for (batch_features,) in dataloader:
                # Predict rewards
                batch_rewards = self.model(batch_features)

                # Add rewards to list
                rewards.extend(batch_rewards.flatten().tolist())

        # Convert to PyTorch
        rewards = torch.tensor(rewards)

        return rewards
