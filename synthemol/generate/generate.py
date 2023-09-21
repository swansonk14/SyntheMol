"""Generate molecules combinatorially using a Monte Carlo tree search guided by a molecular property predictor."""
from datetime import datetime
from pathlib import Path
from typing import Literal

import pandas as pd
import torch
import wandb
from chemprop.nn_utils import initialize_weights
from tap import tapify

from synthemol.constants import (
    BUILDING_BLOCKS_PATH,
    FINGERPRINT_TYPES,
    MODEL_TYPES,
    OPTIMIZATION_TYPES,
    REACTION_TO_BUILDING_BLOCKS_PATH,
    REAL_BUILDING_BLOCK_ID_COL,
    RL_MODEL_TYPES,
    RL_PREDICTION_TYPES,
    SCORE_COL,
    SMILES_COL
)
from synthemol.models import RLModelChemprop, RLModelRDKit
from synthemol.reactions import (
    Reaction,
    REACTIONS,
    load_and_set_allowed_reaction_building_blocks,
    set_all_building_blocks
)
from synthemol.generate.generator import Generator
from synthemol.generate.utils import create_model_scoring_fn, save_generated_molecules


def generate(
        search_type: Literal['mcts', 'rl'],
        save_dir: Path,
        model_types: list[MODEL_TYPES],
        model_paths: list[Path],
        fingerprint_types: list[FINGERPRINT_TYPES],
        model_weights: tuple[float, ...] = (1.0,),
        building_blocks_path: Path = BUILDING_BLOCKS_PATH,
        reaction_to_building_blocks_path: Path | None = REACTION_TO_BUILDING_BLOCKS_PATH,
        building_blocks_id_column: str = REAL_BUILDING_BLOCK_ID_COL,
        building_blocks_score_column: str = SCORE_COL,
        building_blocks_smiles_column: str = SMILES_COL,
        reactions: tuple[Reaction, ...] = REACTIONS,
        max_reactions: int = 1,
        n_rollout: int = 10,
        explore_weight: float = 10.0,
        num_expand_nodes: int | None = None,
        rl_model_type: RL_MODEL_TYPES = 'rdkit',
        rl_prediction_type: RL_PREDICTION_TYPES = 'classification',
        rl_temperature: float = 0.1,
        rl_temperature_similarity_target: float | None = None,
        rl_train_frequency: int = 10,
        rl_train_epochs: int = 5,
        num_workers: int = 0,
        use_gpu: bool = False,
        optimization: OPTIMIZATION_TYPES = 'maximize',
        rng_seed: int = 0,
        no_building_block_diversity: bool = False,
        store_nodes: bool = False,
        save_frequency: int = 1000,
        verbose: bool = False,
        replicate: bool = False,
        wandb_log: bool = False,
        wandb_project_name: str = 'synthemol',
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
    :param model_weights: List of weights for each model/ensemble in model_paths for defining the reward function.
    :param building_blocks_path: Path to CSV file containing molecular building blocks.
    :param reaction_to_building_blocks_path: Path to PKL file containing mapping from REAL reactions to allowed building blocks.
    :param building_blocks_id_column: Name of the column containing IDs for each building block.
    :param building_blocks_score_column: Name of column containing scores for each building block.
    :param building_blocks_smiles_column: Name of the column containing SMILES for each building block.
    :param reactions: A tuple of reactions that combine molecular building blocks.
    :param max_reactions: Maximum number of reactions that can be performed to expand building blocks into molecules.
    :param n_rollout: The number of times to run the generation process.
    :param explore_weight: The hyperparameter that encourages exploration.
    :param num_expand_nodes: The number of child nodes to include when expanding a given node. If None, all child nodes will be included.
    :param rl_model_type: The type of RL model to use. 'rdkit' = MLP RDKIT model. 'chemprop' = pretrained Chemprop model.
    :param rl_prediction_type: The type of prediction made by the RL model, which determines the loss function.
                               'classification' = binary classification. 'regression' = regression.
    :param rl_temperature: The temperature parameter for the softmax function used to select building blocks.
                           Higher temperature means more exploration. If rl_temperature_similarity_target is provided,
                           the temperature is adjusted based on generated molecule diversity.
    :param rl_temperature_similarity_target: If provided, adjusts the temperature to obtain the maximally scoring molecules
                                             that are at most this similar to previously generated molecules. Starts with
                                             the temperature provided by rl_temperature.
    :param rl_train_frequency: The number of rollouts between each training step of the RL model.
    :param rl_train_epochs: The number of epochs to train the RL model for each training step.
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
    # Check lengths of model arguments match
    if len({len(arg) for arg in [model_types, model_paths, fingerprint_types, model_weights]}) != 1:
        raise ValueError('model_types, model_paths, fingerprint_types, and model_weights must all have the same length.')

    # Create save directory
    save_dir.mkdir(parents=True, exist_ok=True)

    # Load building blocks
    print('Loading building blocks...')

    # Optionally alter building blocks loading to precisely replicate previous experiments
    if replicate:
        # Change score loading dtype to ensure numerical precision
        building_block_data = pd.read_csv(building_blocks_path, dtype={building_blocks_score_column: str})
        building_block_data[building_blocks_score_column] = building_block_data[building_blocks_score_column].astype(float)

        # Reorder reactions
        old_reactions_order = [275592, 22, 11, 527, 2430, 2708, 240690, 2230, 2718, 40, 1458, 271948, 27]
        reactions = tuple(sorted(reactions, key=lambda reaction: old_reactions_order.index(reaction.id)))

        # Deduplicate building blocks by SMILES
        building_block_data.drop_duplicates(subset=building_blocks_smiles_column, inplace=True)
    # Otherwise, load the building blocks normally
    else:
        building_block_data = pd.read_csv(building_blocks_path)

    print(f'Loaded {len(building_block_data):,} building blocks')

    # Ensure unique building block IDs
    if building_block_data[building_blocks_id_column].nunique() != len(building_block_data):
        raise ValueError('Building block IDs are not unique.')

    # Map building blocks SMILES to IDs, IDs to SMILES, and SMILES to scores
    building_block_smiles_to_id = dict(zip(
        building_block_data[building_blocks_smiles_column],
        building_block_data[building_blocks_id_column]
    ))
    building_block_id_to_smiles = dict(zip(
        building_block_data[building_blocks_id_column],
        building_block_data[building_blocks_smiles_column]
    ))
    building_block_smiles_to_score = dict(zip(
        building_block_data[building_blocks_smiles_column],
        building_block_data[building_blocks_score_column]
    ))

    print(f'Found {len(building_block_smiles_to_id):,} unique building blocks')

    # Optionally, set up Weights & Biases logging and log building block stats
    if wandb_log:
        # Set up Weights & Biases run name
        if wandb_run_name is None:
            wandb_run_name = f'{search_type}' + (f'_{rl_model_type}' if search_type == 'rl' else '')
            for model_type, fingerprint_type, model_weight in zip(
                    model_types,
                    fingerprint_types,
                    model_weights
            ):
                wandb_run_name += (f'_{model_weight}_{model_type}' +
                                   (f'_{fingerprint_type}' if fingerprint_type != 'none' else ''))

        # Initialize Weights & Biases logging
        wandb.init(
            project=wandb_project_name,
            name=wandb_run_name,
            config={
                'search_type': search_type,
                'save_dir': save_dir,
                'model_paths': model_paths,
                'model_types': model_types,
                'fingerprint_types': fingerprint_types,
                'model_weights': model_weights,
                'building_blocks_path': building_blocks_path,
                'reaction_to_building_blocks_path': reaction_to_building_blocks_path,
                'building_blocks_id_column': building_blocks_id_column,
                'building_blocks_score_column': building_blocks_score_column,
                'building_blocks_smiles_column': building_blocks_smiles_column,
                'max_reactions': max_reactions,
                'n_rollout': n_rollout,
                'explore_weight': explore_weight,
                'num_expand_nodes': num_expand_nodes,
                'rl_model_type': rl_model_type,
                'rl_prediction_type': rl_prediction_type,
                'rl_temperature': rl_temperature,
                'rl_temperature_similarity_target': rl_temperature_similarity_target,
                'rl_train_frequency': rl_train_frequency,
                'rl_train_epochs': rl_train_epochs,
                'optimization': optimization,
                'rng_seed': rng_seed,
                'no_building_block_diversity': no_building_block_diversity,
                'store_nodes': store_nodes,
                'verbose': verbose,
                'replicate': replicate
            }
        )

        # Log building block count
        wandb.log({'Building Block Count': len(building_block_smiles_to_id)})

        # Log building block score histogram
        building_block_scores = [[score] for score in building_block_smiles_to_score.values()]
        table = wandb.Table(data=building_block_scores, columns=['Score'])
        wandb.log({'building_block_scores': wandb.plot.histogram(table, 'Score', title='Building Block Scores')})

    # Set all building blocks for each reaction
    set_all_building_blocks(
        reactions=reactions,
        building_blocks=set(building_block_smiles_to_id)
    )

    # Optionally, set allowed building blocks for each reaction
    if reaction_to_building_blocks_path is not None:
        print('Loading and setting allowed building blocks for each reaction...')
        load_and_set_allowed_reaction_building_blocks(
            reactions=reactions,
            reaction_to_reactant_to_building_blocks_path=reaction_to_building_blocks_path
        )

    # Set up device for model
    if use_gpu:
        device = torch.device('cuda')
    else:
        device = torch.device('cpu')

    # Define model scoring function
    print('Loading models and creating model scoring function...')
    model_scoring_fn = create_model_scoring_fn(
        model_paths=model_paths,
        model_types=model_types,
        fingerprint_types=fingerprint_types,
        model_weights=model_weights,
        device=device,
        smiles_to_score=building_block_smiles_to_score
    )

    # Set up RL model if applicable
    if search_type == 'rl':
        torch.manual_seed(rng_seed)

        if rl_model_type == 'rdkit':
            rl_model = RLModelRDKit(
                num_workers=num_workers,
                num_epochs=rl_train_epochs,
                device=device
            )
        elif rl_model_type.startswith('chemprop'):
            if len(model_paths) > 1 and rl_model_type == 'chemprop_pretrained':
                raise ValueError('Cannot use pretrained RL Chemprop model with multiple model paths.')

            # TODO: Fix this so that chemprop_scratch works without a pretrained Chemprop model
            if model_types[0] != 'chemprop':
                raise ValueError('For RL Chemprop, the first model in model_paths must be a Chemprop model.')

            rl_model = RLModelChemprop(
                model_path=model_paths[0],
                num_workers=num_workers,
                device=device
            )

            if rl_model_type == 'chemprop_scratch':
                initialize_weights(rl_model.model)
        else:
            raise ValueError(f'Invalid RL model type: {rl_model_type}')
    else:
        rl_model = None

    # Set up Generator
    print('Setting up generator...')
    generator = Generator(
        search_type=search_type,
        building_block_smiles_to_id=building_block_smiles_to_id,
        max_reactions=max_reactions,
        scoring_fn=model_scoring_fn,
        explore_weight=explore_weight,
        num_expand_nodes=num_expand_nodes,
        rl_temperature=rl_temperature,
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
        wandb_log=wandb_log
    )

    # Search for molecules
    print(f'Generating molecules for {n_rollout:,} rollouts...')
    start_time = datetime.now()
    molecules_save_path = save_dir / 'molecules.csv'

    for rollout_start in range(0, n_rollout, save_frequency):
        # Set rollout end
        rollout_end = min(rollout_start + save_frequency, n_rollout)

        # Generate molecules
        print(f'Running rollouts {rollout_start:,} through {rollout_end - 1:,}...')
        generator.generate(n_rollout=rollout_end - rollout_start)

        # Get full molecule nodes
        nodes = generator.get_full_molecule_nodes()

        # Save molecules
        print('Saving molecules...')
        save_generated_molecules(
            nodes=nodes,
            building_block_id_to_smiles=building_block_id_to_smiles,
            save_path=molecules_save_path
        )

    # Compute, print, and save stats
    stats = {
        'generation_time': datetime.now() - start_time,
        'num_nonzero_reaction_molecules': len(nodes),
        'approx_num_nodes_searched': generator.approx_num_nodes_searched
    }

    print(f'Generation time = {stats["generation_time"]}')
    print(f'Number of full molecule, nonzero reaction nodes = {stats["num_nonzero_reaction_molecules"]:,}')
    print(f'Approximate total number of nodes searched = {stats["approx_num_nodes_searched"]:,}')

    if store_nodes:
        stats['num_nodes_searched'] = generator.num_nodes_searched
        print(f'Total number of nodes searched = {stats["num_nodes_searched"]:,}')

    pd.DataFrame(data=[stats]).to_csv(save_dir / 'generation_stats.csv', index=False)

    # Log rollout stats
    if wandb_log:
        wandb.log({
            'Generation Time': stats['generation_time'].total_seconds(),
            'Number of Generated Molecules': stats['num_nonzero_reaction_molecules'],
            'Approximate Number of Nodes Searched': stats['approx_num_nodes_searched']
        } | ({'Exact Number of Nodes Searched': stats['num_nodes_searched']} if store_nodes else {}))

        node_scores = [[node.P] for node in nodes]
        table = wandb.Table(data=node_scores, columns=['Score'])
        wandb.log({'generation_scores': wandb.plot.histogram(table, 'Score', title='Generated Molecule Scores')})

    # Save RL model if applicable
    # TODO: fix RL model saving
    if rl_model is not None:
        rl_model.save(save_dir / 'rl_model.pt')


def generate_command_line() -> None:
    """Run generate function from command line."""
    tapify(generate)
