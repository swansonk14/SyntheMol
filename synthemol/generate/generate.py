"""Generate molecules combinatorially using a Monte Carlo tree search guided by a molecular property predictor."""
from datetime import datetime
from pathlib import Path

import pandas as pd
from tap import tapify

from synthemol.constants import (
    BUILDING_BLOCKS_PATH,
    FINGERPRINT_TYPES,
    MODEL_TYPES,
    OPTIMIZATION_TYPES,
    REACTION_TO_BUILDING_BLOCKS_PATH,
    REAL_BUILDING_BLOCK_ID_COL,
    SCORE_COL,
    SMILES_COL
)
from synthemol.reactions import (
    Reaction,
    REACTIONS,
    load_and_set_allowed_reaction_building_blocks,
    set_all_building_blocks
)
from synthemol.generate.generator import Generator
from synthemol.generate.utils import create_model_scoring_fn, save_generated_molecules


def generate(
        model_path: Path,
        model_type: MODEL_TYPES,
        save_dir: Path,
        building_blocks_path: Path = BUILDING_BLOCKS_PATH,
        fingerprint_type: FINGERPRINT_TYPES | None = None,
        reaction_to_building_blocks_path: Path | None = REACTION_TO_BUILDING_BLOCKS_PATH,
        building_blocks_id_column: str = REAL_BUILDING_BLOCK_ID_COL,
        building_blocks_score_column: str = SCORE_COL,
        building_blocks_smiles_column: str = SMILES_COL,
        reactions: tuple[Reaction] = REACTIONS,
        max_reactions: int = 1,
        n_rollout: int = 10,
        explore_weight: float = 10.0,
        num_expand_nodes: int | None = None,
        optimization: OPTIMIZATION_TYPES = 'maximize',
        rng_seed: int = 0,
        no_building_block_diversity: bool = False,
        store_nodes: bool = False,
        verbose: bool = False,
        replicate: bool = False
) -> None:
    """Generate molecules combinatorially using a Monte Carlo tree search guided by a molecular property predictor.

    :param model_path: Path to a directory of model checkpoints or to a specific PKL or PT file containing a trained model.
    :param model_type: Type of model to train.
    :param building_blocks_path: Path to CSV file containing molecular building blocks.
    :param save_dir: Path to directory where the generated molecules will be saved.
    :param fingerprint_type: Type of fingerprints to use as input features.
    :param reaction_to_building_blocks_path: Path to PKL file containing mapping from REAL reactions to allowed building blocks.
    :param building_blocks_id_column: Name of the column containing IDs for each building block.
    :param building_blocks_score_column: Name of column containing scores for each building block.
    :param building_blocks_smiles_column: Name of the column containing SMILES for each building block.
    :param reactions: A tuple of reactions that combine molecular building blocks.
    :param max_reactions: Maximum number of reactions that can be performed to expand building blocks into molecules.
    :param n_rollout: The number of times to run the generation process.
    :param explore_weight: The hyperparameter that encourages exploration.
    :param num_expand_nodes: The number of child nodes to include when expanding a given node. If None, all child nodes will be included.
    :param optimization: Whether to maximize or minimize the score.
    :param rng_seed: Seed for random number generators.
    :param no_building_block_diversity: Whether to turn off the score modification that encourages diverse building blocks.
    :param store_nodes: Whether to store in memory all the nodes of the search tree.
                        This doubles the speed of the search but significantly increases
                        the memory usage (e.g., 450 GB for 20,000 rollouts instead of 600 MB).
    :param replicate: This is necessary to replicate the results from the paper, but otherwise should not be used
                      since it limits the potential choices of building blocks.
    :param verbose: Whether to print out additional information during generation.
    """
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

    # Define model scoring function
    print('Loading models and creating model scoring function...')
    model_scoring_fn = create_model_scoring_fn(
        model_path=model_path,
        model_type=model_type,
        fingerprint_type=fingerprint_type,
        smiles_to_score=building_block_smiles_to_score
    )

    # Set up Generator
    print('Setting up generator...')
    generator = Generator(
        building_block_smiles_to_id=building_block_smiles_to_id,
        max_reactions=max_reactions,
        scoring_fn=model_scoring_fn,
        explore_weight=explore_weight,
        num_expand_nodes=num_expand_nodes,
        optimization=optimization,
        reactions=reactions,
        rng_seed=rng_seed,
        no_building_block_diversity=no_building_block_diversity,
        store_nodes=store_nodes,
        verbose=verbose,
        replicate=replicate
    )

    # Search for molecules
    print('Generating molecules...')
    start_time = datetime.now()
    nodes = generator.generate(n_rollout=n_rollout)

    # Compute, print, and save stats
    stats = {
        'mcts_time': datetime.now() - start_time,
        'num_nonzero_reaction_molecules': len(nodes),
        'approx_num_nodes_searched': generator.approx_num_nodes_searched
    }

    print(f'MCTS time = {stats["mcts_time"]}')
    print(f'Number of full molecule, nonzero reaction nodes = {stats["num_nonzero_reaction_molecules"]:,}')
    print(f'Approximate total number of nodes searched = {stats["approx_num_nodes_searched"]:,}')

    if store_nodes:
        stats['num_nodes_searched'] = generator.num_nodes_searched
        print(f'Total number of nodes searched = {stats["num_nodes_searched"]:,}')

    pd.DataFrame(data=[stats]).to_csv(save_dir / 'mcts_stats.csv', index=False)

    # Save generated molecules
    print('Saving molecules...')
    save_generated_molecules(
        nodes=nodes,
        building_block_id_to_smiles=building_block_id_to_smiles,
        save_path=save_dir / 'molecules.csv'
    )


def generate_command_line() -> None:
    """Run generate function from command line."""
    tapify(generate)
