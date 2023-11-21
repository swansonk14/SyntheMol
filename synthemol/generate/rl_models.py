"""Contains reinforcement learning models for use in generating molecules."""
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Any, Literal

import numpy as np
import torch
import torch.nn as nn
from chemfunc.molecular_fingerprints import compute_fingerprints
from chemprop.models import MoleculeModel
from chemprop.features import BatchMolGraph, MolGraph
from scipy.stats import pearsonr, spearmanr
from sklearn.metrics import mean_squared_error, r2_score
from tqdm import tqdm, trange

from synthemol.constants import RL_PREDICTION_TYPES
from synthemol.generate.model_weights import ModelWeights
from synthemol.generate.node import Node
from synthemol.models import chemprop_build_model, chemprop_load, MLP


# Caching for Chemprop MolGraphs
SMILES_TO_MOL_GRAPH = {}


class RLModel(ABC):
    """A reinforcement learning model for predicting the reward of a molecule or set of molecules."""

    def __init__(
        self,
        prediction_type: RL_PREDICTION_TYPES,
        model_weights: ModelWeights,
        model_paths: list[Path] | None = None,
        max_num_molecules: int = 3,
        features_size: int = 200,
        num_workers: int = 0,
        num_epochs: int = 5,
        batch_size: int = 50,
        learning_rate: float = 1e-3,
        device: torch.device = torch.device("cpu"),
        extended_evaluation: bool = False,
    ) -> None:
        """Initializes the model.

        :param prediction_type: The type of prediction made by the RL model, which determines the loss function.
                                'classification' = binary classification. 'regression' = regression.
        :param model_weights: The weights of the models for each property.
        :param model_paths: A list of paths to PT files or directories of PT files containing models to load if using
                            pretrained models. For directories, the first PT file in the directory will be used.
                            If None, models are initialized randomly.
        :param max_num_molecules: The maximum number of molecules to process at a time.
        :param features_size: The size of the features for each molecule.
        :param num_workers: The number of workers to use for data loading.
        :param num_epochs: The number of epochs to train for.
        :param batch_size: The batch size.
        :param learning_rate: The learning rate.
        :param device: The device to use for the model.
        :param extended_evaluation: Whether to perform extended evaluation of the RL models using all combinations
                                    of numbers of source/target building blocks.
        """
        if model_paths is not None and len(model_paths) != model_weights.num_weights:
            raise ValueError(
                f"Number of model paths ({len(model_paths):,}) must match number "
                f"of model weights ({model_weights.num_weights:,})."
            )

        self.prediction_type = prediction_type
        self.model_weights = model_weights
        self.model_paths = model_paths
        self.max_num_molecules = max_num_molecules
        self.features_size = features_size
        self.total_features_size = max_num_molecules * features_size
        self.num_workers = num_workers
        self.num_epochs = num_epochs
        self.batch_size = batch_size
        self.learning_rate = learning_rate
        self.device = device
        self.extended_evaluation = extended_evaluation

        self.train_source_nodes: list[Node] = []
        self.train_target_nodes: list[Node] = []
        self.test_source_nodes: list[Node] = []
        self.test_target_nodes: list[Node] = []

        self.smiles_to_features: dict[str, torch.Tensor] = {}

        # Build or load models
        if model_paths is None:
            # Build models
            self.models = [self.build_model() for _ in range(model_weights.num_weights)]
        else:
            # Get model path or first model path from each ensemble
            model_paths = [
                model_path
                if model_path.is_file()
                else sorted(model_path.glob("**/*.pt"))[0]
                for model_path in model_paths
            ]

            # Load models
            self.models = [self.load_model(model_path) for model_path in model_paths]

        # Set optimizer
        self.optimizers = [
            torch.optim.Adam(model.parameters(), lr=self.learning_rate)
            for model in self.models
        ]

        # Set loss function
        if self.prediction_type == "classification":
            self.loss_fn = nn.BCEWithLogitsLoss()
            self.eval_loss_fn = nn.BCELoss()
        elif self.prediction_type == "regression":
            self.loss_fn = self.eval_loss_fn = nn.MSELoss()
        else:
            raise ValueError(
                f"Prediction type {self.prediction_type} is not supported."
            )

    def compute_rdkit_features(
        self, molecule_tuples: list[tuple[str]], average_across_tuple: bool = False
    ) -> torch.Tensor:
        """Computes the RDKit features for molecules in each tuple of molecules.

        :param molecule_tuples: A list of tuples of SMILES strings representing one or more molecules.
        :param average_across_tuple: Whether to average the features across the molecules in each tuple.
        :return: A tensor containing the features for molecules in each tuple of molecules.
                 If average_across_tuple is True, then the tensor will have shape (num_tuples, features_size).
                 If average_across_tuple is False, then the tensor will have shape (num_tuples, total_features_size)
                 where total_features_size = max_num_molecules * features_size.
        """
        # Get all molecules
        molecule_nums = [len(molecules) for molecules in molecule_tuples]
        all_molecules = [
            molecule for molecules in molecule_tuples for molecule in molecules
        ]

        # Determine unknown molecules
        unknown_molecules = [
            molecule
            for molecule in all_molecules
            if molecule not in self.smiles_to_features
        ]

        # Compute features for unknown molecules and add to dictionary
        if unknown_molecules:
            unknown_features = compute_fingerprints(
                mols=unknown_molecules, fingerprint_type="rdkit"
            )
            unknown_features = torch.from_numpy(unknown_features)

            for i, molecule in enumerate(unknown_molecules):
                self.smiles_to_features[molecule] = unknown_features[i]

        # Get all molecules and their features
        all_features = torch.stack(
            [self.smiles_to_features[molecule] for molecule in all_molecules]
        )

        # Set up feature vectors for each combination of molecules
        features_dim = (
            self.features_size if average_across_tuple else self.total_features_size
        )
        features = torch.zeros((len(molecule_tuples), features_dim))
        index = 0
        for i, molecule_num in enumerate(molecule_nums):
            # Get features for all molecules in the tuple
            molecules_features = all_features[
                index : index + molecule_num
            ]  # (molecule_num, features_size)

            # Average features across molecules in the tuple if desired
            if average_across_tuple:
                molecules_features = torch.mean(
                    molecules_features, dim=0
                )  # (features_size,)
            else:
                molecules_features = (
                    molecules_features.flatten()
                )  # (total_features_size,)

            # Set features for tuple
            features[i, : len(molecules_features)] = molecules_features
            index += molecule_num

        return features

    def buffer(self, source_node: Node, target_node: Node) -> None:
        """Adds a new source/target node pair to the buffer (test set).

        :param source_node: A node on the pathway to target_node.
        :param target_node: The target node, which is the best node found during the search
                            along the pathway from source_node.
        """
        self.test_source_nodes.append(source_node)
        self.test_target_nodes.append(target_node)

    def test_to_train(self) -> None:
        """Moves the test set of new source/target node pairs to the training set."""
        self.train_source_nodes += self.test_source_nodes
        self.train_target_nodes += self.test_target_nodes
        self.test_source_nodes = []
        self.test_target_nodes = []

    def train(self) -> None:
        """Trains the model on the examples in the training set."""
        # Set models to train mode
        self.train_mode()

        # Get dataloader
        dataloader = self.get_dataloader(
            molecule_tuples=[node.molecules for node in self.train_source_nodes],
            rewards=[node.individual_scores for node in self.train_target_nodes],
            shuffle=True,
        )

        # Loop over epochs
        for _ in trange(self.num_epochs, desc="Training RL model", leave=False):
            # Loop over batches of molecule features and rewards
            for batch_data, batch_rewards in dataloader:
                # Move batch rewards to device
                batch_rewards = batch_rewards.to(self.device)

                # Loop over models and train each one on the batch
                for model_index, (model, optimizer) in enumerate(
                    zip(self.models, self.optimizers)
                ):
                    # Make predictions
                    predictions = self.run_model(
                        model=model, batch_data=batch_data
                    ).squeeze(dim=-1)

                    # Compute loss
                    loss = self.loss_fn(predictions, batch_rewards[:, model_index])

                    # Backpropagate
                    model.zero_grad()
                    loss.backward()
                    optimizer.step()

    # TODO: also evaluate on weighted average property prediction?
    def evaluate(self, split: Literal["train", "test"]) -> dict[str, float]:
        """Evaluates the model on the train or test set.

        :param split: Whether to evaluate on train molecules or test molecules.
        :return: A dictionary of metrics.
        """
        # Select split
        if split == "train":
            source_nodes = self.train_source_nodes
            target_nodes = self.train_target_nodes
        elif split == "test":
            source_nodes = self.test_source_nodes
            target_nodes = self.test_target_nodes
        else:
            raise ValueError(f"Split type {split} is not supported.")

        # Get molecule tuples and rewards (here, rewards are weighted average scores)
        molecule_tuples = [node.molecules for node in source_nodes]
        rewards = [node.individual_scores for node in target_nodes]

        # Make predictions
        predictions = self.predict_individual_scores(
            molecule_tuples=molecule_tuples
        )  # (num_molecules, num_properties)

        # Convert to numpy
        predictions = predictions.numpy()  # (num_molecules, num_properties)
        rewards = np.array(rewards)  # (num_molecules, num_properties)

        # Determine source/target number of building blocks
        source_num_building_blocks = np.array(
            [node.num_building_blocks for node in source_nodes]
        )
        target_num_building_blocks = np.array(
            [node.num_building_blocks for node in target_nodes]
        )

        # Determine source/target number of reactions
        source_num_reactions = np.array([node.num_reactions for node in source_nodes])
        target_num_reactions = np.array([node.num_reactions for node in target_nodes])

        # Get unique source/target number of building blocks (None option for all)
        unique_source_num_building_blocks = [None]
        unique_target_num_building_blocks = [None]

        # Get unique source/target number of reactions (None option for all)
        unique_source_num_reactions = [None]
        unique_target_num_reactions = [None]

        # For extended evaluation, include all combinations of source/target number of building blocks and reactions
        if self.extended_evaluation:
            unique_source_num_building_blocks += sorted(set(source_num_building_blocks))
            unique_target_num_building_blocks += sorted(set(target_num_building_blocks))
            unique_source_num_reactions += sorted(set(source_num_reactions))
            unique_target_num_reactions += sorted(set(target_num_reactions))

        # Set up results dictionary
        results = {}
        split_name = split.title()

        # Get statistics for each property
        for model_index in range(self.model_weights.num_weights):
            # Get predictions and rewards for this property
            model_name = self.model_weights.model_names[model_index]
            model_predictions = predictions[:, model_index]
            model_rewards = rewards[:, model_index]

            # Get statistics for each unique source/target number of building blocks and reactions
            for source_num_bb in unique_source_num_building_blocks:
                for target_num_bb in unique_target_num_building_blocks:
                    for source_num_react in unique_source_num_reactions:
                        for target_num_react in unique_target_num_reactions:
                            # Get mask for source/target number of building blocks and reactions
                            mask = np.ones(len(source_num_building_blocks), dtype=bool)
                            bb_string = ""
                            react_string = ""

                            if source_num_bb is not None:
                                bb_string += f" {source_num_bb} Source BBs"
                                mask &= source_num_building_blocks == source_num_bb

                            if target_num_bb is not None:
                                bb_string += f" {target_num_bb} Target BBs"
                                mask &= target_num_building_blocks == target_num_bb

                            if source_num_react is not None:
                                react_string += f" {source_num_react} Source Reactions"
                                mask &= source_num_reactions == source_num_react

                            if target_num_react is not None:
                                react_string += f" {target_num_react} Target Reactions"
                                mask &= target_num_reactions == target_num_react

                            # Skip if no examples
                            if not np.any(mask):
                                continue

                            # Get predictions and rewards for source/target number of building blocks
                            predictions_masked = model_predictions[mask]
                            rewards_masked = model_rewards[mask]

                            # Create description of this subset
                            description = (
                                f"{model_name} {split_name}{bb_string}{react_string}"
                            )

                            # Evaluate predictions
                            results |= {
                                f"RL {description} Loss": self.eval_loss_fn(
                                    torch.from_numpy(predictions_masked),
                                    torch.from_numpy(rewards_masked),
                                ).item(),
                                f"RL {description} Mean Squared Error": mean_squared_error(
                                    rewards_masked, predictions_masked
                                ),
                            }

                            if (len(np.unique(predictions_masked)) > 2) & (
                                len(np.unique(rewards_masked)) > 2
                            ):
                                results |= {
                                    f"RL {description} R^2": r2_score(
                                        rewards_masked, predictions_masked
                                    ),
                                    f"RL {description} PearsonR": pearsonr(
                                        predictions_masked, rewards_masked
                                    )[0],
                                    f"RL {description} SpearmanR": spearmanr(
                                        predictions_masked, rewards_masked
                                    )[0],
                                }

        return results

    def predict_individual_scores(
        self, molecule_tuples: list[tuple[str]]
    ) -> torch.Tensor:
        """Predicts the individual scores for each tuple of molecules.

        :param molecule_tuples: A list of tuples of SMILES strings representing one or more molecules.
        :return: A 2D tensor containing the predicted individual scores for
                 each tuple of molecules (num_molecules, num_properties).
        """
        # Set models to eval mode
        self.eval_mode()

        # Get dataloader
        dataloader = self.get_dataloader(molecule_tuples=molecule_tuples, shuffle=False)

        # Loop over batches of molecules and make reward predictions
        predictions = []

        with torch.no_grad():
            for batch_data, _ in tqdm(
                dataloader, desc="Predicting RL model", leave=False
            ):
                # Predict rewards
                batch_preds = self.run_models(
                    batch_data=batch_data
                )  # (num_molecules, num_properties)

                # Add predictions to list
                predictions.append(batch_preds.cpu())

        # Concatenate predictions
        predictions = torch.cat(predictions, dim=0)  # (num_molecules, num_properties)

        return predictions

    def predict(self, molecule_tuples: list[tuple[str]]) -> torch.Tensor:
        """Predicts the weighted average reward for each tuple of molecules.

        :param molecule_tuples: A list of tuples of SMILES strings representing one or more molecules.
        :return: A 1D tensor of predicted weighted average rewards for each tuple of molecules.
        """
        # Predict individual scores
        individual_preds = self.predict_individual_scores(
            molecule_tuples=molecule_tuples
        )  # (num_molecules, num_properties)

        # Move weights to tensor
        weights = torch.tensor(self.model_weights.weights).view(
            1, -1
        )  # (1, num_properties)

        # Compute weighted average
        predictions = torch.sum(weights * individual_preds, dim=1)  # (num_molecules,)

        return predictions

    def save(self, path: Path) -> None:
        """Saves the model to the given path.

        :param path: The path to a PT file where the model will be saved.
        """
        torch.save(self, path)

    @classmethod
    def load(cls, path: Path) -> "RLModel":
        """Loads the model from the given path.

        :param path: The path to a PT file containing a model to load.
        """
        return torch.load(path)

    @property
    def train_size(self) -> int:
        """Returns the number of examples in the training set."""
        return len(self.train_source_nodes)

    @property
    def test_size(self) -> int:
        """Returns the number of examples in the test set (buffer)."""
        return len(self.test_source_nodes)

    @abstractmethod
    def build_model(self) -> nn.Module:
        """Builds a model (for predicting an individual property) from scratch.

        :return: A model with randomly initialized weights.
        """
        pass

    @abstractmethod
    def load_model(self, model_path: Path) -> nn.Module:
        """Loads a pretrained model (for predicting an individual property).

        :param model_path: The path to a PT file containing a model to load.
        :return: A model with pretrained weights.
        """
        pass

    def train_mode(self) -> None:
        """Set models to train mode."""
        for model in self.models:
            model.train()

    def eval_mode(self) -> None:
        """Set models to eval mode."""
        for model in self.models:
            model.eval()

    @abstractmethod
    def run_model(self, model: nn.Module, batch_data: Any) -> torch.Tensor:
        """Makes predictions using a model.

        :param model: A model.
        :param batch_data: A Tensor containing the data for the batch.
        :return: A 1D tensor containing the model's predictions.
        """
        pass

    def run_models(self, batch_data: Any) -> torch.Tensor:
        """Makes predictions by computing the models' predictions.

        :param batch_data: Data for the batch.
        :return: A 2D tensor containing the models' predictions (batch_size, num_properties).
        """
        # Compute individual model predictions
        predictions = torch.stack(
            [self.run_model(model, batch_data) for model in self.models]
        )  # (num_properties, batch_size)

        # Transpose to (batch_size, num_properties)
        predictions = predictions.transpose(0, 1)

        return predictions

    @property
    def num_models(self) -> int:
        """Returns the number of models."""
        return len(self.models)

    @abstractmethod
    def get_dataloader(
        self,
        molecule_tuples: list[tuple[str]],
        rewards: list[list[float]] | None = None,
        shuffle: bool = False,
    ) -> torch.utils.data.DataLoader:
        """Returns a dataloader for the given molecules.

        :param molecule_tuples: A list of tuples of SMILES strings representing one or more molecules.
        :param rewards: A list of lists of rewards for each tuple of molecules of shape (num_molecules, num_properties).
        :param shuffle: Whether to shuffle the data.
        :return: A dataloader.
        """
        pass


def rl_mlp_collate_fn(
    data: list[tuple[torch.Tensor], tuple[torch.Tensor, torch.Tensor]]
) -> tuple[torch.Tensor, torch.Tensor | None]:
    """Collates a batch of data for the RL MLP model.

    :param data: A list of tuples of features (and possibly rewards).
    :return: A tuple of features and rewards (or None).
    """
    # Get features and rewards
    if len(data[0]) == 1:
        features = torch.stack(
            [element[0] for element in data]
        )  # (num_molecules, features_size)
        rewards = None
    else:
        features, rewards = zip(*data)
        features = torch.stack(features)  # (num_molecules, features_size)
        rewards = torch.stack(rewards)  # (num_molecules, num_properties)

    return features, rewards


class RLModelMLP(RLModel):
    def __init__(
        self,
        prediction_type: RL_PREDICTION_TYPES,
        model_weights: ModelWeights,
        model_paths: list[Path] | None = None,
        max_num_molecules: int = 3,
        features_size: int = 200,
        hidden_dim: int = 100,
        num_layers: int = 2,
        num_workers: int = 0,
        num_epochs: int = 5,
        batch_size: int = 50,
        learning_rate: float = 1e-3,
        device: torch.device = torch.device("cpu"),
        extended_evaluation: bool = False,
    ) -> None:
        """Initializes the RL MLP model.

        :param prediction_type: The type of prediction made by the RL model, which determines the loss function.
                                'classification' = binary classification. 'regression' = regression.
        :param model_weights: The weights of the models for each property.
        :param model_paths: A list of paths to PT files or directories of PT files containing models to load if using
                            pretrained models. For directories, the first PT file in the directory will be used.
                            If None, models are initialized randomly.
        :param max_num_molecules: The maximum number of molecules to process at a time.
        :param features_size: The size of the features for each molecule.
        :param hidden_dim: The dimensionality of the hidden layers.
        :param num_layers: The number of layers.
        :param num_workers: The number of workers to use for data loading.
        :param num_epochs: The number of epochs to train for.
        :param batch_size: The batch size.
        :param learning_rate: The learning rate.
        :param device: The device to use for the model.
        :param extended_evaluation: Whether to perform extended evaluation of the RL models using all combinations
                                    of numbers of source/target building blocks.
        """
        # Store MLP-specific parameters
        self.hidden_dim = hidden_dim
        self.output_dim = 1
        self.num_layers = num_layers
        self._model = MLP(
            input_dim=max_num_molecules * features_size,
            hidden_dim=hidden_dim,
            output_dim=1,
            num_layers=num_layers,
            sigmoid=prediction_type == "classification",
            device=device,
        ).to(device)

        super().__init__(
            prediction_type=prediction_type,
            model_weights=model_weights,
            model_paths=model_paths,
            max_num_molecules=max_num_molecules,
            features_size=features_size,
            num_workers=num_workers,
            num_epochs=num_epochs,
            batch_size=batch_size,
            learning_rate=learning_rate,
            device=device,
            extended_evaluation=extended_evaluation,
        )

    def build_model(self) -> MLP:
        """Builds a model (for predicting an individual property) from scratch.

        :return: A model with randomly initialized weights.
        """
        return MLP(
            input_dim=self.max_num_molecules * self.features_size,
            hidden_dim=self.hidden_dim,
            output_dim=1,
            num_layers=self.num_layers,
            sigmoid=self.prediction_type == "classification",
            device=self.device,
        )

    def load_model(self, model_path: Path) -> MLP:
        """Loads a pretrained model (for predicting an individual property).

        :param model_path: The path to a PT file containing a model to load.
        :return: A model with pretrained weights.
        """
        # TODO: implement this
        # Need to load Chemprop feed-forward model and then match layers in state dicts
        raise NotImplementedError

    def run_model(self, model: MLP, batch_data: torch.Tensor) -> torch.Tensor:
        """Makes predictions using a model.

        :param model: A model.
        :param batch_data: A Tensor containing the data for the batch.
        :return: A 1D tensor containing the model's prediction.
        """
        return model(batch_data).squeeze(dim=-1)

    def get_dataloader(
        self,
        molecule_tuples: list[tuple[str]],
        rewards: list[list[float]] | None = None,
        shuffle: bool = False,
    ) -> torch.utils.data.DataLoader:
        """Returns a dataloader for the given molecules.

        :param molecule_tuples: A list of tuples of SMILES strings representing one or more molecules.
        :param rewards: A list of lists of rewards for each tuple of molecules of shape (num_molecules, num_properties).
        :param shuffle: Whether to shuffle the data.
        :return: A dataloader.
        """
        # Compute features in bulk across all molecules
        features = self.compute_rdkit_features(
            molecule_tuples=molecule_tuples, average_across_tuple=False
        )

        # Set up data based on presence of rewards
        if rewards is None:
            data = (features,)
        else:
            data = (features, torch.tensor(rewards))

        # Create dataset
        dataset = torch.utils.data.TensorDataset(*data)

        # Create dataloader
        dataloader = torch.utils.data.DataLoader(
            dataset=dataset,
            batch_size=self.batch_size,
            shuffle=shuffle,
            num_workers=self.num_workers,
            collate_fn=rl_mlp_collate_fn,
        )

        return dataloader


class RLChempropMoleculeDataset(torch.utils.data.Dataset):
    def __init__(
        self,
        molecule_tuples: list[tuple[str, ...]],
        features: torch.Tensor | None = None,
        rewards: list[list[float]] | None = None,
    ) -> None:
        """Initializes the dataset.

        :param molecule_tuples: A list of tuples of SMILES strings representing one or more molecules.
        :param features: A tensor containing the features for each molecule in each tuple of molecules
                         (num_tuples, total_features_size).
        :param rewards: A list of lists of rewards for each tuple of molecules of shape (num_molecules, num_properties).
        """
        # Store molecule tuples
        self.molecule_tuples = molecule_tuples

        # Handle features
        if features is None:
            self.features = [None] * len(molecule_tuples)
        else:
            self.features = features

        # Handle rewards
        if rewards is None:
            self.rewards = [None] * len(molecule_tuples)
        else:
            self.rewards = rewards

        # Add each molecule to the MolGraph cache
        for molecule_tuple in tqdm(molecule_tuples, desc="Caching MolGraphs"):
            for molecule in molecule_tuple:
                if molecule not in SMILES_TO_MOL_GRAPH:
                    SMILES_TO_MOL_GRAPH[molecule] = MolGraph(molecule)

        # Get MolGraphs for each molecule
        self.mol_graphs_tuples = [
            tuple(SMILES_TO_MOL_GRAPH[molecule] for molecule in molecule_tuple)
            for molecule_tuple in molecule_tuples
        ]

    def __len__(self) -> int:
        """Returns the number of tuples of molecules."""
        return len(self.molecule_tuples)

    def __getitem__(
        self, index: int
    ) -> tuple[tuple[MolGraph, ...], torch.Tensor | None, list[float] | None]:
        """Returns a MolGraph tuple and the corresponding features (or None) and reward (or None)."""
        return self.mol_graphs_tuples[index], self.features[index], self.rewards[index]


def rl_chemprop_collate_fn(
    data: list[tuple[tuple[MolGraph, ...], torch.Tensor | None, list[float] | None]]
) -> tuple[tuple[list[BatchMolGraph], list[np.ndarray] | None], torch.Tensor | None]:
    """Collates data into a batch for the RL Chemprop model.

    :param data: A list of tuples of MolGraph tuples, features, and rewards.
    :return: A tuple with a tuple containing a list of BatchMolGraph and features, and a tensor of rewards (or None).
    """
    # Get molecules, features, and rewards
    mol_graph_tuples, features, rewards = zip(*data)

    # Create BatchMolGraph from MolGraphs
    batch_mol_graph = BatchMolGraph(
        [
            mol_graph
            for mol_graph_tuple in mol_graph_tuples
            for mol_graph in mol_graph_tuple
        ]
    )

    # Create atom and bond scopes for BatchMolGraph based on mol_graph_tuples
    a_scope, b_scope = [], []
    n_atoms = n_bonds = 0
    for mol_graph_tuple in mol_graph_tuples:
        n_atoms_tuple = sum(mol_graph.n_atoms for mol_graph in mol_graph_tuple)
        n_bonds_tuple = sum(mol_graph.n_bonds for mol_graph in mol_graph_tuple)

        a_scope.append((n_atoms, n_atoms_tuple))
        b_scope.append((n_bonds, n_bonds_tuple))

        n_atoms += n_atoms_tuple
        n_bonds += n_bonds_tuple

    # Set BatchMolGraph atom and bond scopes
    batch_mol_graph.a_scope = a_scope
    batch_mol_graph.b_scope = b_scope

    # Set up features
    if features[0] is None:
        features = None
    else:
        features = [feats.numpy() for feats in features]

    # Set up rewards
    if rewards[0] is None:
        rewards = None
    else:
        rewards = torch.tensor(rewards)

    # Set up batch data for Chemprop
    batch_data = ([batch_mol_graph], features)

    return batch_data, rewards


class RLModelChemprop(RLModel):
    def __init__(
        self,
        use_rdkit_features: bool,
        prediction_type: RL_PREDICTION_TYPES,
        model_weights: ModelWeights,
        model_paths: list[Path] | None = None,
        max_num_molecules: int = 3,
        features_size: int = 200,
        num_workers: int = 0,
        num_epochs: int = 5,
        batch_size: int = 50,
        learning_rate: float = 1e-3,
        device: torch.device = torch.device("cpu"),
        extended_evaluation: bool = False,
    ) -> None:
        """Initializes the RL Chemprop model.

        :param use_rdkit_features: Whether to use RDKit fingerprints as features.
        :param prediction_type: The type of prediction made by the RL model, which determines the loss function.
                                'classification' = binary classification. 'regression' = regression.
        :param model_weights: The weights of the models for each property.
        :param model_paths: A list of paths to PT files or directories of PT files containing models to load if using
                            pretrained models. For directories, the first PT file in the directory will be used.
                            If None, models are initialized randomly.
        :param max_num_molecules: The maximum number of molecules to process at a time.
        :param features_size: The size of the features for each molecule.
        :param num_workers: The number of workers to use for data loading.
        :param num_epochs: The number of epochs to train for.
        :param batch_size: The batch size.
        :param learning_rate: The learning rate.
        :param device: The device to use for the model.
        :param extended_evaluation: Whether to perform extended evaluation of the RL models using all combinations
                                    of numbers of source/target building blocks.
        """
        # Store whether to use RDKit features
        self.use_rdkit_features = use_rdkit_features

        super().__init__(
            prediction_type=prediction_type,
            model_weights=model_weights,
            model_paths=model_paths,
            max_num_molecules=max_num_molecules,
            features_size=features_size,
            num_workers=num_workers,
            num_epochs=num_epochs,
            batch_size=batch_size,
            learning_rate=learning_rate,
            device=device,
            extended_evaluation=extended_evaluation,
        )

    def build_model(self) -> MoleculeModel:
        """Builds a model (for predicting an individual property) from scratch.

        :return: A model with randomly initialized weights.
        """
        return chemprop_build_model(
            dataset_type=self.prediction_type,
            use_rdkit_features=self.use_rdkit_features,
            rdkit_features_size=self.features_size,
            property_name="rl_objective",
        ).to(self.device)

    def load_model(self, model_path: Path) -> nn.Module:
        """Loads a pretrained model (for predicting an individual property).

        :param model_path: The path to a PT file containing a model to load.
        :return: A model with pretrained weights.
        """
        return chemprop_load(model_path=model_path, device=self.device)

    def run_model(
        self,
        model: MoleculeModel,
        batch_data: tuple[list[BatchMolGraph], list[np.ndarray] | None],
    ) -> torch.Tensor:
        """Makes predictions using the model.

        :param model: A model.
        :param batch_data: A tuple containing a list of BatchMolGraphs and features for the batch.
        :return: A 1D tensor containing the model's prediction.
        """
        # Unpack batch data
        batch_mol_graphs, features = batch_data

        # Get predictions
        predictions = model(batch=batch_mol_graphs, features_batch=features).squeeze(
            dim=-1
        )

        # Return predictions
        return predictions

    def get_dataloader(
        self,
        molecule_tuples: list[tuple[str]],
        rewards: list[list[float]] | None = None,
        shuffle: bool = False,
    ) -> torch.utils.data.DataLoader:
        """Returns a dataloader for the given molecules.

        :param molecule_tuples: A list of tuples of SMILES strings representing one or more molecules.
        :param rewards: A list of lists of rewards for each tuple of molecules of shape (num_molecules, num_properties).
        :param shuffle: Whether to shuffle the data.
        :return: A dataloader.
        """
        # Compute RDKit features if needed
        if self.use_rdkit_features:
            # Compute features and average across molecules in a tuple for a tensor of shape (num_tuples, features_size)
            features = self.compute_rdkit_features(
                molecule_tuples=molecule_tuples, average_across_tuple=True
            )
        else:
            features = None

        # Create dataset
        dataset = RLChempropMoleculeDataset(
            molecule_tuples=molecule_tuples, features=features, rewards=rewards
        )

        # Create dataloader
        dataloader = torch.utils.data.DataLoader(
            dataset=dataset,
            batch_size=self.batch_size,
            shuffle=shuffle,
            num_workers=self.num_workers,
            collate_fn=rl_chemprop_collate_fn,
        )

        return dataloader
