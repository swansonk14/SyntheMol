"""Generate molecules combinatorially using a search guided by a molecular property predictor."""
from copy import deepcopy
from datetime import datetime
from pathlib import Path
from typing import Literal

import pandas as pd
import torch
import wandb
from tap import tapify

from synthemol.constants import (
    BUILDING_BLOCKS_PATH,
    CHEMICAL_SPACES,
    FINGERPRINT_TYPES,
    FEATURES_SIZE_MAPPING,
    ID_COL,
    SCORE_TYPES,
    OLD_REACTION_ORDER,
    OLD_REACTIONS,
    OPTIMIZATION_TYPES,
    REACTION_TO_BUILDING_BLOCKS_PATH,
    RL_MODEL_TYPES,
    RL_PREDICTION_TYPES,
    SCORE_COL,
    SMILES_COL,
)
from synthemol.generate.score_weights import ScoreWeights
from synthemol.generate.rl_models import RLModelChemprop, RLModelMLP
from synthemol.reactions import (
    CHEMICAL_SPACE_TO_REACTIONS,
    load_and_set_allowed_reaction_building_blocks,
    set_all_building_blocks,
)
from synthemol.generate.generator import Generator
from synthemol.generate.scorer import MoleculeScorer
from synthemol.generate.utils import parse_comparator_string, save_generated_molecules


def generate(
    search_type: Literal["mcts", "rl"],
    save_dir: Path,
    score_types: list[SCORE_TYPES],
    score_model_paths: list[str] | None = None,
    score_fingerprint_types: list[str] | None = None,
    score_names: list[str] | None = None,
    base_score_weights: list[float] | None = None,
    success_thresholds: list[str] | None = None,
    chemical_spaces: tuple[CHEMICAL_SPACES, ...] = ("real",),
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
    rolling_average_weight: float = 0.98,
    rl_model_type: RL_MODEL_TYPES = "chemprop",
    rl_model_fingerprint_type: str | None = None,
    rl_model_paths: list[Path] | None = None,
    rl_prediction_types: tuple[RL_PREDICTION_TYPES] = ("classification",),
    rl_base_temperature: float = 0.1,
    rl_temperature_similarity_target: float = 0.6,
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
    replicate_mcts: bool = False,
    replicate_rl: bool = False,
    wandb_log: bool = False,
    wandb_project_name: str = "synthemol",
    wandb_run_name: str | None = None,
    h2o_solvents: bool = False,
) -> None:
    """Generate molecules combinatorially using a search guided by a molecular property predictor.

    :param search_type: The type of search to perform. 'mcts' = Monte Carlo tree search. 'rl' = Reinforcement learning.
    :param save_dir: Path to directory where the generated molecules will be saved.
    :param score_types: List of types of scores to score molecules.
    :param score_model_paths: For score types that are model-based ("random_forest" and "chemprop"), the corresponding
        score model path should be a path to a directory of model checkpoints (ensemble)
        or to a specific PKL or PT file containing a trained model with a single output.
        For score types that are not model-based, the corresponding score model path must be "None".
        If all score types are not model-based, this argument can be None.
    :param score_fingerprint_types: For score types that are model-based and require fingerprints as input, the corresponding
        fingerprint type should be the type of fingerprint (e.g., "rdkit").
        For model-based scores that don't require fingerprints or non-model-based scores,
        the corresponding fingerprint type must be "None".
        If all score types do not require fingerprints, this argument can be None.
    :param score_names: List of names for each score. If None, scores will be named "Score 1", "Score 2", etc.
    :param base_score_weights: Initial weights for each score for defining the reward function.
        If None, defaults to equal weights for each score.
    :param success_thresholds: The threshold for each score for defining success of the form ">= 0.5".
        If provided, the score weights will be dynamically set to maximize joint success
        across all scores.
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
    :param rolling_average_weight: The weight to use for the rolling average of similarity (dynamic temperature)
        and success (dynamic score weights).
    :param rl_model_type: The type of RL model to use. 'mlp' = MLP model. 'chemprop' = Chemprop model. 
    :param rl_model_fingerprint_types: The corresponding fingerprint type for the RL model. MLP models require the fingerprint type to be defined, 
        but Chemprop RL models can have a fingerprint type of None. 
    :param rl_model_paths: List of paths with each path pointing to a PT file containing a trained model that will be
        used as the initial weights for the RL models. If None, RL models are trained from scratch.
    :param rl_prediction_types: The types of predictions made by the RL models, which determines the loss functions.
        'classification' = binary classification. 'regression' = regression.
    :param rl_base_temperature: The initial temperature parameter for the softmax function used to select building blocks.
        Higher temperature means more exploration. If rl_temperature_similarity_target is provided,
        the temperature is adjusted based on generated molecule diversity.
    :param rl_temperature_similarity_target: Adjusts the temperature to obtain the maximally scoring molecules
        that are at most this similar to previously generated molecules. Starts with the temperature provided
        by rl_base_temperature. If -1, the temperature is not adjusted.
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
    :param replicate_mcts: This is necessary to replicate the results from the  MCTS paper
        but otherwise should not be used since it limits the potential choices of building blocks.
    :param replicate_rl: This is necessary to replicate the result from the RL paper
        but otherwise should not be used since it does not include reaction fixes.
    :param wandb_log: Whether to log results to Weights & Biases.
    :param wandb_project_name: The name of the Weights & Biases project to log results to.
    :param wandb_run_name: The name of the Weights & Biases run to log results to.
    :param h2o_solvents: Whether to concatenate H2O solvent features with the molecule features during prediction.
    """
    # Convert score_model_paths to Path/None
    # TODO: change tapify to allow list[Path | Literal["None"] | None]

    if rl_model_type == "mlp" and rl_model_fingerprint_type is None:
        raise ValueError("MLP RL models must have a fingerprint type that is not None")

    if score_model_paths is not None:
        score_model_paths: list[Path | None] = [
            Path(score_model_path) if score_model_path not in ("None", None) else None
            for score_model_path in score_model_paths
        ]

    # Convert score_fingerprint_types to FINGERPRINT_TYPES/None
    # TODO: change tapify to allow list[FINGERPRINT_TYPES | Literal["None"] | None]
    if score_fingerprint_types is not None:
        score_fingerprint_types: list[FINGERPRINT_TYPES | None] = [
            score_fingerprint_type
            if score_fingerprint_type not in ("None", None)
            else None
            for score_fingerprint_type in score_fingerprint_types
        ]

    # Change type of building blocks score columns for compatibility with Pandas
    building_blocks_score_columns = list(building_blocks_score_columns)

    # Get number of scores
    num_scores = len(score_types)

    # Set up default score weights as equal weighting
    if base_score_weights is None:
        base_score_weights = [1 / num_scores] * num_scores

    # Check lengths of model arguments match
    for arg in [
        score_types,
        score_model_paths,
        rl_model_paths,
        score_fingerprint_types,
        score_names,
        base_score_weights,
        building_blocks_score_columns,
        success_thresholds,
    ]:
        if arg is not None and len(arg) != num_scores:
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

    # Set up score weights, with option for dynamic weights if success_thresholds is provided
    if success_thresholds is not None:
        success_comparators = tuple(
            parse_comparator_string(success_threshold)
            for success_threshold in success_thresholds
        )
        score_weights = ScoreWeights(
            base_score_weights=base_score_weights,
            immutable=False,
            score_names=score_names,
        )
    else:
        success_comparators = None
        score_weights = ScoreWeights(
            base_score_weights=base_score_weights,
            immutable=True,
            score_names=score_names,
        )

    # Create save directory
    save_dir.mkdir(parents=True, exist_ok=True)

    # Get (unique) reactions
    reactions = tuple(
        dict.fromkeys(
            reaction
            for chemical_space in sorted(set(chemical_spaces))
            for reaction in CHEMICAL_SPACE_TO_REACTIONS[chemical_space]
        )
    )

    # Optionally make changes to reactions to replicate previous experiments
    if replicate_mcts or replicate_rl:
        # Only keep old reactions (and copy them to avoid modifying original reactions)
        reactions = [
            deepcopy(reaction)
            for reaction in reactions
            if reaction.reaction_id in OLD_REACTIONS
        ]

        # Undo changes to old reactions
        # TODO: OH to O[H]?
        for reaction in reactions:
            reaction.id = reaction.reaction_id

            if reaction.id == 240790:
                reaction.post_reactions = None

            if reaction.id == 271948:
                reaction.post_reactions = None
                reaction.reactants = reaction.reactants[::-1]

    print(f"Using {len(reactions):,} reactions")

    # Load building blocks
    print("Loading building blocks...")

    # Optionally make changes to building blocks to replicate previous experiments
    if replicate_mcts:
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

        # Reorder old reactions
        reactions = tuple(
            sorted(
                reactions,
                key=lambda reaction: OLD_REACTION_ORDER.index(reaction.reaction_id),
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
            wandb_run_name = (
                f"{search_type}"
                f"{f'_{rl_model_type}' if search_type == 'rl' else ''}"
                f"_{'_'.join(score_types)}"
                f"_{'_'.join(chemical_spaces)}"
            )

        # Initialize Weights & Biases logging
        wandb.init(
            project=wandb_project_name,
            name=wandb_run_name,
            config={
                "search_type": search_type,
                "save_dir": save_dir,
                "score_types": score_types,
                "score_model_paths": score_model_paths,
                "score_fingerprint_types": score_fingerprint_types,
                "score_names": score_names,
                "base_score_weights": base_score_weights,
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
                "rolling_average_weight": rolling_average_weight,
                "rl_model_type": rl_model_type,
                "rl_fingerprint_type": rl_model_fingerprint_type,
                "rl_model_paths": rl_model_paths,
                "rl_prediction_types": rl_prediction_types,
                "rl_base_temperature": rl_base_temperature,
                "rl_temperature_similarity_target": rl_temperature_similarity_target,
                "rl_train_frequency": rl_train_frequency,
                "rl_train_epochs": rl_train_epochs,
                "rl_extended_evaluation": rl_extended_evaluation,
                "h2o_solvents": h2o_solvents,
                "optimization": optimization,
                "rng_seed": rng_seed,
                "no_building_block_diversity": no_building_block_diversity,
                "store_nodes": store_nodes,
                "verbose": verbose,
                "replicate": replicate_mcts,
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
                columns=list(score_weights.score_names),
            )

            # Log building block score histograms
            for score_name in score_weights.score_names:
                histogram_name = f"{chemical_space} Building Block {score_name} Scores"
                wandb.log(
                    {
                        histogram_name: wandb.plot.histogram(
                            wandb_bb_table, score_name, title=histogram_name,
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
    device = torch.device("cuda" if use_gpu else "cpu")

    # Define scorer
    print("Creating scorer...")
    scorer = MoleculeScorer(
        score_types=score_types,
        score_weights=score_weights,
        model_paths=score_model_paths,
        fingerprint_types=score_fingerprint_types,
        h2o_solvents=h2o_solvents,
        device=device,
        smiles_to_scores=building_block_smiles_to_scores,
    )

    # Set up RL model if applicable
    if search_type == "rl":
        # Seed PyTorch for reproducibility
        torch.manual_seed(rng_seed)

        # Set up RL model args
        rl_model_args = {
            "prediction_types": rl_prediction_types,
            "score_weights": score_weights,
            "model_paths": rl_model_paths,
            "num_workers": num_workers,
            "num_epochs": rl_train_epochs,
            "device": device,
            "extended_evaluation": rl_extended_evaluation,
            "features_type": rl_model_fingerprint_type,
            "features_size": FEATURES_SIZE_MAPPING[rl_model_fingerprint_type],
            "h2o_solvents": h2o_solvents
        }

        # Select RL model class and update RL model args
        if rl_model_type == "mlp":
            rl_model_class = RLModelMLP
        elif rl_model_type == "chemprop":
            rl_model_class = RLModelChemprop
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
        score_weights=score_weights,
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
        rolling_average_weight=rolling_average_weight,
        replicate=replicate_mcts,
        wandb_log=wandb_log,
    )

    # Search for molecules
    print(f"Generating molecules for {n_rollout:,} rollouts...")
    start_time = datetime.now()
    molecules_save_path = save_dir / "molecules.csv"
    nodes = []

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
            score_names=score_weights.score_names,
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
            data=node_scores, columns=list(score_weights.score_names)
        )

        # Log histograms of scores of generated molecules
        for score_name in score_weights.score_names:
            histogram_name = f"Generated Molecule {score_name} Scores"
            wandb.log(
                {
                    histogram_name: wandb.plot.histogram(
                        wandb_gen_table, score_name, title=histogram_name,
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
