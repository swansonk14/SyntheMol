"""Generate molecules combinatorially using a Monte Carlo tree search guided by a molecular property predictor."""
from datetime import datetime
from pathlib import Path
from typing import Literal

import numpy as np
import pandas as pd
import torch
import wandb
from tap import tapify

from synthemol.constants import (
    BUILDING_BLOCKS_PATH,
    FINGERPRINT_TYPES,
    ID_COL,
    MODEL_TYPES,
    OPTIMIZATION_TYPES,
    REACTION_TO_BUILDING_BLOCKS_PATH,
    RL_MODEL_TYPES,
    RL_PREDICTION_TYPES,
    SCORE_COL,
    SMILES_COL,
)
from synthemol.generate.model_weights import ModelWeights
from synthemol.generate.rl_models import RLModelChemprop, RLModelMLP
from synthemol.reactions import (
    CHEMICAL_SPACE_TO_REACTIONS,
    load_and_set_allowed_reaction_building_blocks,
    set_all_building_blocks,
)
from synthemol.generate.generator import Generator
from synthemol.generate.scorer import MoleculeScorer
from synthemol.generate.utils import parse_success_threshold, save_generated_molecules


# TODO: once tuple[Literal] fix is done in tap, add ('none',) default to fingerprint_types and add literal to reaction_sets
def generate(
    search_type: Literal["mcts", "rl"],
    save_dir: Path,
    model_types: list[MODEL_TYPES],
    model_paths: list[Path],
    fingerprint_types: list[FINGERPRINT_TYPES],
    model_names: list[str] = None,
    base_model_weights: list[float] | None = None,
    success_thresholds: list[str] | None = None,
    chemical_spaces: tuple[str, ...] = ("real",),
    building_blocks_paths: tuple[Path, ...] = (BUILDING_BLOCKS_PATH,),
    reaction_to_building_blocks_paths: tuple[Path, ...] = (
        REACTION_TO_BUILDING_BLOCKS_PATH,
    ),
    building_blocks_id_column: str = ID_COL,
    building_blocks_score_columns: tuple[str, ...] = (SCORE_COL,),
    building_blocks_smiles_column: str = SMILES_COL,
    max_reactions: int = 1,
    n_rollout: int = 10,
    explore_weight: float = 10.0,
    num_expand_nodes: int | None = None,
    rl_model_type: RL_MODEL_TYPES = "mlp_rdkit",
    rl_pretrained: bool = False,
    rl_prediction_type: RL_PREDICTION_TYPES = "classification",
    rl_base_temperature: float = 0.1,
    rl_temperature_similarity_target: float = 0.5,
    rl_train_frequency: int = 10,
    rl_train_epochs: int = 5,
    rl_extended_evaluation: bool = False,
    num_workers: int = 0,
    use_gpu: bool = False,
    optimization: OPTIMIZATION_TYPES = "maximize",
    rng_seed: int = 0,
    no_building_block_diversity: bool = False,
    store_nodes: bool = False,
    save_frequency: int = 1000,
    verbose: bool = False,
    replicate: bool = False,
    wandb_log: bool = False,
    wandb_project_name: str = "synthemol",
    wandb_run_name: str | None = None,
) -> None:
    """Generate molecules combinatorially using a Monte Carlo tree search guided by a molecular property predictor.

    :param search_type: The type of search to perform. 'mcts' = Monte Carlo tree search. 'rl' = Reinforcement learning.
    :param save_dir: Path to directory where the generated molecules will be saved.
    :param model_types: List of types of models provided by model_paths.
    :param model_paths: List of paths with each path pointing to a directory of model checkpoints (ensemble)
                        or to a specific PKL or PT file containing a trained model.
                        Note: All models must have a single output.
    :param fingerprint_types: List of types of fingerprints to use as input features for the model_paths.
    :param model_names: List of names for each model/ensemble in model_paths.
                        If None, models will be named "Model 1", "Model 2", etc.
    :param base_model_weights: Initial weights for each model/ensemble in model_paths for defining the reward function.
                               If None, defaults to equal weights for each model/ensemble.
    :param success_thresholds: The threshold for each model/ensemble in model_paths for defining success of the form "> 0.5".
                               If provided, the model weights will be dynamically set to maximize joint success
                               across all models/ensembles.
    :param chemical_spaces: A tuple of names of reaction sets to use. 'real' = Enamine REAL Space reactions.
                            'wuxi' = WuXi GalaXi reactions. 'custom' = Custom reactions (in synthemol/reactions/custom.py).
    :param building_blocks_paths: Paths to CSV files containing molecular building blocks.
    :param reaction_to_building_blocks_paths: Paths to PKL files containing mapping from reactions to allowed building blocks.
    :param building_blocks_id_column: Name of the column containing IDs for each building block.
    :param building_blocks_score_columns: Name of columns containing scores for each building block.
    :param building_blocks_smiles_column: Name of the column containing SMILES for each building block.
    :param max_reactions: Maximum number of reactions that can be performed to expand building blocks into molecules.
    :param n_rollout: The number of times to run the generation process.
    :param explore_weight: The hyperparameter that encourages exploration.
    :param num_expand_nodes: The number of child nodes to include when expanding a given node. If None, all child nodes will be included.
    :param rl_model_type: The type of RL model to use. 'mlp_rdkit' = MLP RDKit model.
                          'chemprop' = Chemprop model. 'chemprop_rdkit' = Chemprop RDKit model.
    :param rl_pretrained: Whether to use pretrained model checkpoints from model_paths to initialize the RL models.
                          If True, the first model in each ensemble in model_paths will be used as the initialization.
                          If False, RL models are randomly initialized.
    :param rl_prediction_type: The type of prediction made by the RL model, which determines the loss function.
                               'classification' = binary classification. 'regression' = regression.
    :param rl_base_temperature: The initial temperature parameter for the softmax function used to select building blocks.
                                Higher temperature means more exploration. If rl_temperature_similarity_target is provided,
                                the temperature is adjusted based on generated molecule diversity.
    :param rl_temperature_similarity_target: Adjusts the temperature to obtain the maximally scoring molecules
                                             that are at most this similar to previously generated molecules. Starts with
                                             the temperature provided by rl_base_temperature.
                                             If -1, the temperature is not adjusted.
    :param rl_train_frequency: The number of rollouts between each training step of the RL model.
    :param rl_train_epochs: The number of epochs to train the RL model for each training step.
    :param rl_extended_evaluation: Whether to perform extended evaluation of the RL model after each training step.
    :param num_workers: The number of workers for RL model data loading.
    :param use_gpu: Whether to use GPU for model training/prediction. Only affects PyTorch-based models (not sklearn).
    :param optimization: Whether to maximize or minimize the score.
    :param rng_seed: Seed for random number generators.
    :param no_building_block_diversity: Whether to turn off the score modification that encourages diverse building blocks.
    :param store_nodes: Whether to store in memory all the nodes of the search tree.
                        This doubles the speed of the search but significantly increases
                        the memory usage (e.g., 450 GB for 20,000 rollouts instead of 600 MB).
    :param save_frequency: The number of rollouts between each save of the generated molecules.
    :param verbose: Whether to print out additional information during generation.
    :param replicate: This is necessary to replicate the results from the paper, but otherwise should not be used
                      since it limits the potential choices of building blocks.
    :param wandb_log: Whether to log results to Weights & Biases.
    :param wandb_project_name: The name of the Weights & Biases project to log results to.
    :param wandb_run_name: The name of the Weights & Biases run to log results to.
    """
    # Change type of building blocks score columns for compatibility with Pandas
    building_blocks_score_columns = list(building_blocks_score_columns)

    # Get number of models
    num_models = len(model_types)

    # Set up default model weights as equal weighting
    if base_model_weights is None:
        base_model_weights = [1 / num_models] * num_models

    # Check lengths of model arguments match
    for arg in [
        model_types,
        model_paths,
        fingerprint_types,
        model_names,
        base_model_weights,
        building_blocks_score_columns,
        success_thresholds,
    ]:
        if arg is not None and len(arg) != num_models:
            raise ValueError("Model parameters have the same length.")

    # Check lengths of chemical space arguments match
    if (
        len(
            {
                len(arg)
                for arg in [
                    chemical_spaces,
                    building_blocks_paths,
                    reaction_to_building_blocks_paths,
                ]
            }
        )
        != 1
    ):
        raise ValueError(
            "chemical_spaces, building_blocks_paths, and reaction_to_building_blocks_paths "
            "must have the same length."
        )

    # Set RL temperature similarity target to None if not desired
    if rl_temperature_similarity_target == -1:
        rl_temperature_similarity_target = None

    # Set up model weights, with option for dynamic weights if success_thresholds is provided
    if success_thresholds is not None:
        success_comparators = tuple(
            parse_success_threshold(success_threshold)
            for success_threshold in success_thresholds
        )
        model_weights = ModelWeights(
            base_model_weights=base_model_weights,
            immutable=False,
            model_names=model_names,
        )
    else:
        success_comparators = None
        model_weights = ModelWeights(
            base_model_weights=base_model_weights,
            immutable=True,
            model_names=model_names,
        )

    # Create save directory
    save_dir.mkdir(parents=True, exist_ok=True)

    # Get reactions
    reactions = tuple(
        reaction
        for chemical_space in sorted(set(chemical_spaces))
        for reaction in CHEMICAL_SPACE_TO_REACTIONS[chemical_space]
    )

    print(f"Using {len(reactions):,} reactions")

    # Load building blocks
    print("Loading building blocks...")

    # Optionally alter building blocks loading to precisely replicate previous experiments
    if replicate:
        # Ensure only REAL space
        assert chemical_spaces == ("real",)

        # Ensure only one building blocks score column
        assert len(building_blocks_score_columns) == 1
        building_blocks_score_column = building_blocks_score_columns[0]

        # Change score loading dtype to ensure numerical precision
        real_building_block_data = pd.read_csv(
            building_blocks_paths[0],
            dtype={building_blocks_id_column: str, building_blocks_score_column: str},
        )
        real_building_block_data[
            building_blocks_score_column
        ] = real_building_block_data[building_blocks_score_column].astype(float)

        # Reorder reactions
        old_reactions_order = [
            275592,
            22,
            11,
            527,
            2430,
            2708,
            240690,
            2230,
            2718,
            40,
            1458,
            271948,
            27,
        ]
        reactions = tuple(
            sorted(
                reactions, key=lambda reaction: old_reactions_order.index(reaction.id)
            )
        )

        # Deduplicate building blocks by SMILES
        real_building_block_data.drop_duplicates(
            subset=building_blocks_smiles_column, inplace=True
        )

        # Put building block data in a dictionary
        chemical_space_to_building_block_data = {
            chemical_spaces[0]: real_building_block_data
        }
    # Otherwise, load the building blocks normally
    else:
        chemical_space_to_building_block_data = {
            chemical_space: pd.read_csv(
                building_blocks_path, dtype={building_blocks_id_column: str}
            )
            for chemical_space, building_blocks_path in zip(
                chemical_spaces, building_blocks_paths
            )
        }

    # Print building block stats and check that building block IDs are unique
    for (
        chemical_space,
        building_block_data,
    ) in chemical_space_to_building_block_data.items():
        print(f"Loaded {len(building_block_data):,} {chemical_space} building blocks")

        if building_block_data[building_blocks_id_column].nunique() != len(
            building_block_data
        ):
            raise ValueError(
                f"Building block IDs are not unique in {chemical_space} chemical space."
            )

    # Unique set of building block SMILES
    building_block_smiles: set[str] = set.union(
        *[
            set(building_block_data[building_blocks_smiles_column])
            for building_block_data in chemical_space_to_building_block_data.values()
        ]
    )

    print(f"Found {len(building_block_smiles):,} unique building blocks")

    # Map chemical space to building block SMILES to IDs
    chemical_space_to_building_block_smiles_to_id: dict[str, dict[str, str]] = {
        chemical_space: dict(
            zip(
                building_block_data[building_blocks_smiles_column],
                building_block_data[building_blocks_id_column],
            )
        )
        for chemical_space, building_block_data in chemical_space_to_building_block_data.items()
    }

    # Map chemical space to building block ID to SMILES
    chemical_space_to_building_block_id_to_smiles: dict[str, dict[str, str]] = {
        chemical_space: dict(
            zip(
                building_block_data[building_blocks_id_column],
                building_block_data[building_blocks_smiles_column],
            )
        )
        for chemical_space, building_block_data in chemical_space_to_building_block_data.items()
    }

    # Map building block SMILES to scores
    building_block_smiles_to_scores: dict[str, list[float]] = {
        smiles: list(scores)
        for building_block_data in chemical_space_to_building_block_data.values()
        for smiles, scores in zip(
            building_block_data[building_blocks_smiles_column],
            building_block_data[building_blocks_score_columns].itertuples(index=False),
        )
    }

    # Optionally, set up Weights & Biases logging and log building block stats
    if wandb_log:
        # Set up Weights & Biases run name
        if wandb_run_name is None:
            wandb_run_name = f"{search_type}" + (
                f"_{rl_model_type}" if search_type == "rl" else ""
            )

            for model_type, fingerprint_type in zip(model_types, fingerprint_types):
                wandb_run_name += f"_{model_type}" + (
                    f"_{fingerprint_type}" if fingerprint_type != "none" else ""
                )

            wandb_run_name += f'_{"_".join(chemical_spaces)}'

        # Initialize Weights & Biases logging
        wandb.init(
            project=wandb_project_name,
            name=wandb_run_name,
            config={
                "search_type": search_type,
                "save_dir": save_dir,
                "model_paths": model_paths,
                "model_types": model_types,
                "fingerprint_types": fingerprint_types,
                "model_names": model_names,
                "base_model_weights": base_model_weights,
                "success_thresholds": success_thresholds,
                "chemical_spaces": chemical_spaces,
                "building_blocks_paths": building_blocks_paths,
                "reaction_to_building_blocks_paths": reaction_to_building_blocks_paths,
                "building_blocks_id_column": building_blocks_id_column,
                "building_blocks_score_columns": building_blocks_score_columns,
                "building_blocks_smiles_column": building_blocks_smiles_column,
                "max_reactions": max_reactions,
                "n_rollout": n_rollout,
                "explore_weight": explore_weight,
                "num_expand_nodes": num_expand_nodes,
                "rl_model_type": rl_model_type,
                "rl_pretrained": rl_pretrained,
                "rl_prediction_type": rl_prediction_type,
                "rl_base_temperature": rl_base_temperature,
                "rl_temperature_similarity_target": rl_temperature_similarity_target,
                "rl_train_frequency": rl_train_frequency,
                "rl_train_epochs": rl_train_epochs,
                "rl_extended_evaluation": rl_extended_evaluation,
                "optimization": optimization,
                "rng_seed": rng_seed,
                "no_building_block_diversity": no_building_block_diversity,
                "store_nodes": store_nodes,
                "verbose": verbose,
                "replicate": replicate,
            },
        )

        # For each chemical space, log number of building blocks and score histograms
        for (
            chemical_space,
            building_block_data,
        ) in chemical_space_to_building_block_data.items():
            # Log number of building blocks
            wandb.log(
                {f"{chemical_space} Building Block Count": len(building_block_data)}
            )

            # Build table of building block scores
            wandb_bb_table = wandb.Table(
                data=building_block_data[building_blocks_score_columns].values,
                columns=list(model_weights.model_names),
            )

            # Log building block score histograms
            for model_name in model_weights.model_names:
                histogram_name = f"{chemical_space} Building Block {model_name} Scores"
                wandb.log(
                    {
                        histogram_name: wandb.plot.histogram(
                            wandb_bb_table, model_name, title=histogram_name,
                        )
                    }
                )

    # Set all building blocks for each reaction
    set_all_building_blocks(reactions=reactions, building_blocks=building_block_smiles)

    # Set allowed building blocks for each reaction
    print("Loading and setting allowed building blocks for each reaction...")
    load_and_set_allowed_reaction_building_blocks(
        reactions=reactions,
        chemical_spaces=chemical_spaces,
        reaction_to_building_blocks_paths=reaction_to_building_blocks_paths,
    )

    # Set up device for model
    if use_gpu:
        device = torch.device("cuda")
    else:
        device = torch.device("cpu")

    # Define model scoring function
    print("Loading models and creating model scorer...")
    scorer = MoleculeScorer(
        model_paths=model_paths,
        model_types=model_types,
        fingerprint_types=fingerprint_types,
        model_weights=model_weights,
        device=device,
        smiles_to_scores=building_block_smiles_to_scores,
    )

    # Set up RL model if applicable
    if search_type == "rl":
        # Seed PyTorch for reproducibility
        torch.manual_seed(rng_seed)

        # Set up RL model args
        rl_model_args = {
            "prediction_type": rl_prediction_type,
            "model_weights": model_weights,
            "model_paths": model_paths if rl_pretrained else None,
            "num_workers": num_workers,
            "num_epochs": rl_train_epochs,
            "device": device,
            "extended_evaluation": rl_extended_evaluation,
        }

        # Select RL model class and update RL model args
        if rl_model_type == "mlp_rdkit":
            rl_model_class = RLModelMLP
        elif rl_model_type.startswith("chemprop"):
            rl_model_class = RLModelChemprop
            rl_model_args["use_rdkit_features"] = rl_model_type == "chemprop_rdkit"
        else:
            raise ValueError(f"Invalid RL model type: {rl_model_type}")

        # Create RL model
        rl_model = rl_model_class(**rl_model_args)

        # Print RL model architecture
        print(f"RL model architecture: {rl_model.models}")
    else:
        rl_model = None

    # Set up Generator
    print("Setting up generator...")
    generator = Generator(
        search_type=search_type,
        chemical_space_to_building_block_smiles_to_id=chemical_space_to_building_block_smiles_to_id,
        max_reactions=max_reactions,
        scorer=scorer,
        model_weights=model_weights,
        success_comparators=success_comparators,
        explore_weight=explore_weight,
        num_expand_nodes=num_expand_nodes,
        rl_base_temperature=rl_base_temperature,
        rl_temperature_similarity_target=rl_temperature_similarity_target,
        rl_train_frequency=rl_train_frequency,
        optimization=optimization,
        reactions=reactions,
        rng_seed=rng_seed,
        no_building_block_diversity=no_building_block_diversity,
        store_nodes=store_nodes,
        verbose=verbose,
        rl_model=rl_model,
        replicate=replicate,
        wandb_log=wandb_log,
    )

    # Search for molecules
    print(f"Generating molecules for {n_rollout:,} rollouts...")
    start_time = datetime.now()
    molecules_save_path = save_dir / "molecules.csv"

    for rollout_start in range(0, n_rollout, save_frequency):
        # Set rollout end
        rollout_end = min(rollout_start + save_frequency, n_rollout)

        # Generate molecules
        print(f"Running rollouts {rollout_start:,} through {rollout_end - 1:,}...")
        generator.generate(n_rollout=rollout_end - rollout_start)

        # Get full molecule nodes
        nodes = generator.get_full_molecule_nodes()

        # Save molecules
        print("Saving molecules...")
        save_generated_molecules(
            nodes=nodes,
            chemical_space_to_building_block_id_to_smiles=chemical_space_to_building_block_id_to_smiles,
            model_names=model_weights.model_names,
            save_path=molecules_save_path,
        )

    # Compute, print, and save stats
    stats = {
        "generation_time": datetime.now() - start_time,
        "num_nonzero_reaction_molecules": len(nodes),
        "approx_num_nodes_searched": generator.approx_num_nodes_searched,
    }

    print(f'Generation time = {stats["generation_time"]}')
    print(
        f'Number of full molecule, nonzero reaction nodes = {stats["num_nonzero_reaction_molecules"]:,}'
    )
    print(
        f'Approximate total number of nodes searched = {stats["approx_num_nodes_searched"]:,}'
    )

    if store_nodes:
        stats["num_nodes_searched"] = generator.num_nodes_searched
        print(f'Total number of nodes searched = {stats["num_nodes_searched"]:,}')

    pd.DataFrame(data=[stats]).to_csv(save_dir / "generation_stats.csv", index=False)

    # Log rollout stats
    if wandb_log:
        wandb.log(
            {
                "Generation Time": stats["generation_time"].total_seconds(),
                "Number of Generated Molecules": stats[
                    "num_nonzero_reaction_molecules"
                ],
                "Approximate Number of Nodes Searched": stats[
                    "approx_num_nodes_searched"
                ],
            }
            | (
                {"Exact Number of Nodes Searched": stats["num_nodes_searched"]}
                if store_nodes
                else {}
            )
        )

        # Build table of scores of generated molecules
        node_scores = [node.individual_scores for node in nodes]
        wandb_gen_table = wandb.Table(
            data=node_scores, columns=list(model_weights.model_names)
        )

        # Log histograms of scores of generated molecules
        for model_name in model_weights.model_names:
            histogram_name = f"Generated Molecule {model_name} Scores"
            wandb.log(
                {
                    histogram_name: wandb.plot.histogram(
                        wandb_gen_table, model_name, title=histogram_name,
                    )
                }
            )

    # Save RL model
    # TODO: resolve pickling errors
    # if rl_model is not None:
    #     rl_model.save(save_dir / 'rl_model.pt')


def generate_command_line() -> None:
    """Run generate function from command line."""
    tapify(generate)
