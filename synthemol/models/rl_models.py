"""Contains reinforcement learning models for use in generating molecules."""
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Iterator

import torch
import torch.nn as nn
from chemfunc.molecular_fingerprints import compute_fingerprints
from chemprop.data import MoleculeDataLoader, MoleculeDatapoint, MoleculeDataset
from sklearn.metrics import mean_squared_error, r2_score
from tqdm import trange

from synthemol.models import chemprop_load


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

        # Apply sigmoid during prediction
        if not self.training:
            X = torch.sigmoid(X)

        return X


class RLModel(ABC):
    """A reinforcement learning model for predicting the reward of a molecule or set of molecules."""

    def __init__(
            self,
            num_workers: int = 0,
            num_epochs: int = 5,
            batch_size: int = 50,
            learning_rate: float = 1e-3
    ) -> None:
        """Initializes the model.

        :param num_workers: The number of workers to use for data loading.
        :param num_epochs: The number of epochs to train for.
        :param batch_size: The batch size.
        :param learning_rate: The learning rate.
        """
        self.num_workers = num_workers
        self.num_epochs = num_epochs
        self.batch_size = batch_size
        self.learning_rate = learning_rate

        self.molecule_tuples: list[tuple[str]] = []
        self.rewards: list[float] = []

        self.optimizer = torch.optim.Adam(self.model.parameters(), lr=self.learning_rate)
        self.loss_fn = nn.BCEWithLogitsLoss()  # TODO: have option for regression

    def buffer(self, molecules: tuple[str], reward: float) -> None:
        """Adds a training example to the buffer.

        :param molecules: A tuple of SMILES strings representing one or more molecules.
        :param reward: The reward of those molecules (i.e., the score of the final molecule constructed
                       from these molecules and potentially others).
        """
        self.molecule_tuples.append(molecules)
        self.rewards.append(reward)

    def train(self) -> None:
        """Trains the model on the examples in the buffer."""
        # Set model to train mode
        self.model.train()

        # TODO: upsample new examples and downsample old examples

        # Get dataloader
        dataloader = self.get_dataloader(
            molecule_tuples=self.molecule_tuples,
            rewards=self.rewards,
            shuffle=True
        )

        # Loop over epochs
        for _ in trange(self.num_epochs, desc='Training RL model'):
            # Loop over batches of molecule features and rewards
            for batch_data, batch_rewards in dataloader:
                # Make predictions
                predictions = self.model(batch_data).squeeze(dim=-1)

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
        predictions = self.predict(self.molecule_tuples)

        # Convert rewards to tensor
        rewards = torch.tensor(self.rewards)

        # Convert to numpy
        predictions_numpy = predictions.detach().numpy()
        rewards_numpy = rewards.detach().numpy()

        # Evaluate predictions
        results = {
            'RL Loss': self.loss_fn(predictions, rewards).item(),
            'RL Mean Squared Error': mean_squared_error(rewards_numpy, predictions_numpy),
            'RL R^2': r2_score(rewards_numpy, predictions_numpy),
        }

        return results

    def predict(self, molecule_tuples: list[tuple[str]]) -> torch.Tensor:
        """Predicts the reward for each tuple of molecules.

        :param molecule_tuples: A list of tuples of SMILES strings representing one or more molecules.
        :return: A 1D tensor of rewards for each tuple of molecules.
        """
        # Set model to eval mode
        self.model.eval()

        # Get dataloader
        dataloader = self.get_dataloader(
            molecule_tuples=molecule_tuples,
            shuffle=False
        )

        # Loop over batches of molecules and make reward predictions
        rewards = []

        with torch.no_grad():
            from tqdm import tqdm
            for (batch_data,) in tqdm(dataloader):
                # Predict rewards
                batch_rewards = self.model(batch_data)

                # Add rewards to list
                rewards.extend(batch_rewards.flatten().tolist())

        # Convert to PyTorch
        rewards = torch.tensor(rewards)

        return rewards

    @property
    def buffer_size(self) -> int:
        """Returns the number of examples in the buffer."""
        return len(self.molecule_tuples)

    @property
    @abstractmethod
    def model(self) -> nn.Module:
        """Returns the model."""
        pass

    @abstractmethod
    def get_dataloader(
            self,
            molecule_tuples: list[tuple[str]],
            rewards: list[float] | None = None,
            shuffle: bool = False
    ) -> torch.utils.data.DataLoader:
        """Returns a dataloader for the given molecules.

        :param molecule_tuples: A list of tuples of SMILES strings representing one or more molecules.
        :param rewards: A list of rewards for each tuple of molecules.
        :param shuffle: Whether to shuffle the data.
        :return: A dataloader.
        """
        pass


class RLModelRDKit(RLModel):
    def __init__(
            self,
            max_num_molecules: int = 3,
            features_size: int = 200,
            hidden_dim: int = 100,
            num_layers: int = 2,
            num_workers: int = 0,
            num_epochs: int = 5,
            batch_size: int = 50,
            learning_rate: float = 1e-3,
    ) -> None:
        """Initializes the model.

        :param max_num_molecules: The maximum number of molecules to process at a time.
        :param features_size: The size of the features for each molecule.
        :param hidden_dim: The dimensionality of the hidden layers.
        :param num_layers: The number of layers.
        :param num_workers: The number of workers to use for data loading.
        :param num_epochs: The number of epochs to train for.
        :param batch_size: The batch size.
        :param learning_rate: The learning rate.
        """
        self.max_num_molecules = max_num_molecules
        self.features_size = features_size
        self.total_features_size = max_num_molecules * features_size

        self.smiles_to_features: dict[str, torch.Tensor] = {}

        self._model = MLP(
            input_dim=self.total_features_size,
            hidden_dim=hidden_dim,
            output_dim=1,
            num_layers=num_layers
        )

        super().__init__(
            num_workers=num_workers,
            num_epochs=num_epochs,
            batch_size=batch_size,
            learning_rate=learning_rate
        )

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

    @property
    def model(self) -> nn.Module:
        """Returns the model."""
        return self._model

    def get_dataloader(
            self,
            molecule_tuples: list[tuple[str]],
            rewards: list[float] | None = None,
            shuffle: bool = False
    ) -> torch.utils.data.DataLoader:
        """Returns a dataloader for the given molecules.

        :param molecule_tuples: A list of tuples of SMILES strings representing one or more molecules.
        :param rewards: A list of rewards for each tuple of molecules.
        :param shuffle: Whether to shuffle the data.
        :return: A dataloader.
        """
        # Compute features in bulk across all molecules
        features = self.compute_rdkit_features(molecule_tuples)

        # Create dataset input
        dataset_input = (features,) if rewards is None else (features, torch.tensor(rewards))

        # Create dataset
        dataset = torch.utils.data.TensorDataset(*dataset_input)

        # Create dataloader
        dataloader = torch.utils.data.DataLoader(
            dataset=dataset,
            batch_size=self.batch_size,
            shuffle=shuffle,
            num_workers=self.num_workers
        )

        return dataloader


class RLMoleculeDataLoader(MoleculeDataLoader):
    def __iter__(self) -> Iterator[tuple[torch.Tensor, ...]]:
        """Creates an iterator which returns (x, y) or (x,) pairs
        where x is a batch of molecules and y is a batch of rewards."""
        for batch in super(RLMoleculeDataLoader, self).__iter__():
            data = batch.batch_graph()
            targets = batch.targets()

            if targets[0] is None:
                yield (data,)
            else:
                yield data, torch.tensor([target[0] for target in targets])


class RLModelChemprop(RLModel):
    def __init__(
            self,
            model_path: Path,
            num_workers: int = 0,
            num_epochs: int = 5,
            batch_size: int = 50,
            learning_rate: float = 1e-3,
    ) -> None:
        """Initializes the model.

        :param model_path: The path to pretrained Chemprop checkpoint file.
        :param num_workers: The number of workers to use for data loading.
        :param num_epochs: The number of epochs to train for.
        :param batch_size: The batch size.
        :param learning_rate: The learning rate.
        """
        # If model_path is a directory, take the first model
        if model_path.is_dir():
            self.model_path = sorted(model_path.glob('**/*.pt'))[0]
        else:
            self.model_path = model_path

        self._model = chemprop_load(
            model_path=self.model_path
        )

        super().__init__(
            num_workers=num_workers,
            num_epochs=num_epochs,
            batch_size=batch_size,
            learning_rate=learning_rate
        )

    @property
    def model(self) -> nn.Module:
        """Returns the model."""
        return self._model

    def get_dataloader(
            self,
            molecule_tuples: list[tuple[str]],
            rewards: list[float] | None = None,
            shuffle: bool = False
    ) -> torch.utils.data.DataLoader:
        """Returns a dataloader for the given molecules.

        :param molecule_tuples: A list of tuples of SMILES strings representing one or more molecules.
        :param rewards: A list of rewards for each tuple of molecules.
        :param shuffle: Whether to shuffle the data.
        :return: A dataloader.
        """
        if rewards is None:
            rewards = [None] * len(molecule_tuples)
        else:
            rewards = [[float(reward)] for reward in rewards]

        return RLMoleculeDataLoader(
            dataset=MoleculeDataset([
                MoleculeDatapoint(
                    smiles=['.'.join(molecule_tuple)],  # Merge multiple molecules into one string
                    targets=reward,
                ) for molecule_tuple, reward in zip(molecule_tuples, rewards)
            ]),
            num_workers=self.num_workers,
            shuffle=shuffle
        )
