"""Assess the novelty, scores, and diversity of generated molecules."""
from collections import Counter
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from chemfunc import compute_top_similarities

from SyntheMol.constants import SCORE_COL, SMILES_COL


def plot_scores(
        scores: list[float],
        save_dir: Path,
        score_name: str = SCORE_COL
) -> None:
    """Plot score distribution.

    :param scores: A list of scores.
    :param save_dir: The directory where the plot will be saved.
    :param score_name: The name of the score.
    """
    # Plot score distribution
    plt.clf()
    plt.hist(scores, bins=100)
    plt.xlabel(score_name)
    plt.ylabel('Count')
    plt.title(f'{score_name} Distribution')
    plt.savefig(save_dir / f'{score_name}.pdf', bbox_inches='tight')

    # Save score distribution
    fig_data = pd.DataFrame({score_name: scores})
    fig_data.to_csv(save_dir / f'{score_name}.csv', index=False)


def plot_similarity(
        smiles: list[str],
        similarity_type: str,
        save_dir: Path,
        reference_smiles: list[str] | None = None,
        reference_name: str | None = None
) -> None:
    """Plot similarity distribution within a list of SMILES or between that list and a reference list.

    :param smiles: A list of SMILES.
    :param similarity_type: The type of similarity.
    :param save_dir: The directory where the plot will be saved.
    :param reference_smiles: A list of reference SMILES to compare against.
    :param reference_name: The name of the reference list of SMILES.
    """
    # Compute similarities
    max_similarities = compute_top_similarities(
        similarity_type=similarity_type,
        mols=smiles,
        reference_mols=reference_smiles
    )

    # Get reference name
    if reference_name is None:
        reference_name = 'Internal'

    # Get save name
    save_name = save_dir / f'{reference_name.lower()}_{similarity_type}_similarity'

    # Plot diversity distribution
    plt.clf()
    plt.hist(max_similarities, bins=100)
    plt.xlabel(f'Maximum {reference_name} {similarity_type.title()} Similarity')
    plt.ylabel('Count')
    plt.title(f'Maximum {reference_name} {similarity_type.title()} Similarity Distribution')
    plt.savefig(f'{save_name}.pdf', bbox_inches='tight')

    # Save diversity distribution
    fig_data = pd.DataFrame({f'max_{similarity_type}_similarity': max_similarities})
    fig_data.to_csv(f'{save_name}.csv', index=False)


def plot_reactions_numbers(
        num_reactions: list[int],
        save_dir: Path
) -> None:
    """Plot the frequency of each number of reactions.

    :param num_reactions: A list of numbers of reactions per molecule.
    :param save_dir: The directory where the plot will be saved.
    """
    # Get reaction counts
    reaction_counts = Counter(num_reactions)
    min_reactions, max_reactions = min(reaction_counts), max(reaction_counts)
    reaction_nums = range(min_reactions, max_reactions + 1)
    reaction_counts = [reaction_counts[num_reactions] for num_reactions in reaction_nums]

    # Plot reaction counts
    plt.clf()
    plt.bar(reaction_nums, reaction_counts)
    plt.xticks(reaction_nums)
    plt.xlabel('Number of Reactions')
    plt.ylabel('Count')
    plt.title('Number of Reactions')
    plt.savefig(save_dir / 'reaction_numbers.pdf', bbox_inches='tight')

    # Save reaction counts
    fig_data = pd.DataFrame({'reaction_nums': reaction_nums, 'reaction_counts': reaction_counts})
    fig_data.to_csv(save_dir / 'reaction_numbers.csv', index=False)


def plot_reaction_usage(
        data: pd.DataFrame,
        save_dir: Path
) -> None:
    """Plot the frequency with which each reaction is used (unique reactions per molecule).

    :param data: DataFrame containing reaction usage.
    :param save_dir: The directory where the plot will be saved.
    """
    # Get reaction usage
    reaction_columns = [column for column in data.columns if column.startswith('reaction_') and column.endswith('_id')]
    reaction_data = data[reaction_columns]
    reaction_counter = Counter(
        reaction
        for _, reaction_row in reaction_data.iterrows()
        for reaction in {int(reaction) for reaction in reaction_row.dropna()}
    )

    reactions = sorted(reaction_counter)
    reaction_counts = [reaction_counter[reaction] for reaction in reactions]
    xticks = np.arange(len(reaction_counts))

    # Plot reaction usage
    plt.clf()
    plt.bar(xticks, reaction_counts)
    plt.xticks(ticks=xticks, labels=reactions, rotation=45)
    plt.xlabel('Reaction')
    plt.ylabel('Count (# molecules containing the reaction)')
    plt.title('Reaction Counts')
    plt.savefig(save_dir / 'reaction_counts.pdf', bbox_inches='tight')

    # Save reaction usage
    fig_data = pd.DataFrame({'reaction': reactions, 'count': reaction_counts})
    fig_data.to_csv(save_dir / 'reaction_counts.csv', index=False)


def plot_building_block_usage(
        data: pd.DataFrame,
        save_dir: Path
) -> None:
    """Plot the frequency with which each building block is used (unique building blocks per molecule).

    :param data: DataFrame containing building block usage.
    :param save_dir: The directory where the plot will be saved.
    """
    # Get building block usage
    building_block_columns = [
        column
        for column in data.columns
        if column.startswith('building_block_') and column.endswith('_id')
    ]
    building_block_data = data[building_block_columns]
    building_block_counter = Counter(
        building_block
        for _, building_block_row in building_block_data.iterrows()
        for building_block in {int(building_block) for building_block in building_block_row.dropna()}
        if building_block != -1
    )

    building_blocks_with_counts = building_block_counter.most_common()
    building_blocks, building_block_counts = zip(*building_blocks_with_counts)

    # Plot building block usage
    plt.clf()
    plt.scatter(np.arange(len(building_block_counts)), building_block_counts)

    plt.xlabel('Sorted Building Block Index')
    plt.ylabel('Count (# molecules containing the building block)')
    plt.title('Building Block Counts')
    plt.savefig(save_dir / 'building_block_counts.pdf', bbox_inches='tight')

    # Save building block usage
    fig_data = pd.DataFrame({'building_block': building_blocks, 'count': building_block_counts})
    fig_data.to_csv(save_dir / 'building_block_counts.csv', index=False)


def plot_generated_molecule_analysis(
        data_path: Path,
        save_dir: Path,
        reference_paths: list[Path] | None = None,
        score_column: str = SCORE_COL,
        smiles_column: str = SMILES_COL,
        reference_smiles_column: str = SMILES_COL
) -> None:
    """Assess the novelty, scores, and diversity of generated molecules.

    :param data_path: Path to CSV file containing the generated molecules.
    :param save_dir: Directory where the plots will be saved.
    :param reference_paths: Optional list of paths to CSV files containing reference molecules for computing similarity.
    :param score_column: Name of the column containing the scores.
    :param smiles_column: Name of the column containing the SMILES strings.
    :param reference_smiles_column: Name of the column containing the SMILES strings in the reference files.
    """
    # Load generated molecules
    data = pd.read_csv(data_path)

    # Count molecules
    print(f'Number of molecules = {len(data):,}')

    # Create save directory
    save_dir.mkdir(parents=True, exist_ok=True)

    # Score distribution
    plot_scores(
        scores=data[score_column],
        save_dir=save_dir
    )

    # Similarity within generated molecules
    plot_similarity(
        smiles=data[smiles_column],
        similarity_type='tanimoto',
        save_dir=save_dir
    )

    if reference_paths is not None:
        for reference_path in reference_paths:
            # Load reference molecules
            reference_molecules = pd.read_csv(reference_path)

            print(f'Number of reference molecules in {reference_path.stem} = {len(reference_molecules):,}')

            # Similarity between generated molecules and reference molecules
            plot_similarity(
                smiles=data[smiles_column],
                similarity_type='tversky',
                save_dir=save_dir,
                reference_smiles=reference_molecules[reference_smiles_column],
                reference_name=reference_path.stem
            )

    # Number of reactions
    plot_reactions_numbers(
        num_reactions=data['num_reactions'],
        save_dir=save_dir
    )

    # Usage of reactions
    plot_reaction_usage(
        data=data,
        save_dir=save_dir
    )

    # Usage of building blocks
    plot_building_block_usage(
        data=data,
        save_dir=save_dir
    )


if __name__ == '__main__':
    from tap import tapify

    tapify(plot_generated_molecule_analysis)
