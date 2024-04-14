"""Selects a set of diverse, novel hit molecules from a set of molecules."""
from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd
from chemfunc import compute_pairwise_tanimoto_similarities

from synthemol.constants import SCORE_COL, SMILES_COL


def compute_average_maximum_similarity(similarities: np.ndarray) -> float:
    """Computes the average maximum pairwise similarity within a set of molecules.

    :param similarities: A 2D numpy array of pairwise molecule similarities.
    :return: Average maximum pairwise similarity.
    """
    np.fill_diagonal(similarities, -np.inf)
    max_similarities = np.max(similarities, axis=1)
    average_max_similarity = float(np.mean(max_similarities))

    return average_max_similarity


def get_approximate_maximum_independent_set(
    similarities: np.ndarray, threshold: float, num_tries: int = 10
) -> set[int]:
    """Gets a maximal independent set approximating the maximum independent set of molecules within a similarity threshold.

    Note: This does not provide a formal approximation since there are no theoretical guarantees on using
        maximal independent sets to approximate maximum independent sets.

    :param similarities: A 2D numpy array of pairwise molecule similarities.
    :param threshold: Similarity threshold, which is the maximum similarity allowed between molecules in the set.
    :param num_tries: Number of times to run the approximation (i.e., number of maximal independent sets to try).
    :return: The approximate maximum independent set (set of indices of molecules).
    """
    # Create similarities graph
    similarities_matrix = similarities > threshold
    similarities_graph = nx.from_numpy_array(similarities_matrix)

    # Compute approximate maximum independent set
    max_independent_set = set()
    for seed in range(num_tries):
        independent_set = nx.maximal_independent_set(similarities_graph, seed=seed)
        if len(independent_set) > len(max_independent_set):
            max_independent_set = independent_set

    return max_independent_set


def select_molecules(
    data_path: Path,
    save_molecules_path: Path,
    save_analysis_path: Path | None = None,
    smiles_column: str = SMILES_COL,
    score_columns: tuple[str, ...] = (SCORE_COL,),
    novelty_columns: tuple[str, ...] = (
        "train_hits_tversky_nearest_neighbor_similarity",
        "chembl_tversky_nearest_neighbor_similarity",
    ),
    score_comparators: tuple[str, ...] = (">= 0.5",),
    novelty_thresholds: tuple[float, ...] = (0.6, 0.6),
    similarity_threshold: float = 0.6,
    max_rollout: int | None = None,
    select_num: int | None = None,
    sort_column: str | None = None,
    descending: bool = False,
) -> None:
    """Selects a set of diverse, novel hit molecules from a set of molecules

    :param data_path: Path to CSV file containing molecules.
    :param save_molecules_path: Path to a CSV file where the selected molecules will be saved.
    :param save_analysis_path: Optional path to CSV file where the analysis will be saved.
    :param smiles_column: Name of the column containing SMILES.
    :param score_columns: Name of the columns containing scores.
    :param novelty_columns: Name of the columns containing similarity scores compared to known hits (for novelty).
    :param score_comparators: Comparators to use for calculating the hits and similarity (one per score column).
        The comparators should be of the form ">= 0.5".
    :param novelty_thresholds: Thresholds to use for filtering by novelty (one per novelty column).
        Molecules with similarity scores at or below the threshold are considered novel.
    :param similarity_threshold: Threshold to use for calculating the maximum independent set (diverse molecules).
        Molecules with pairwise similarity scores at or below the threshold are considered diverse.
    :param max_rollout: Maximum rollout number to include in the analysis.
    :param select_num: Optional number of molecules to select (otherwise keeps entire maximum independent set).
    :param sort_column: Optional name of the column to sort by before selecting molecules.
        If None, selects molecules in the order of the input file.
    :param descending: Whether to sort in descending order.
    """
    # Check arguments
    if len(score_columns) != len(score_comparators):
        raise ValueError(
            "Number of score columns must match the number of score comparators."
        )

    if len(novelty_columns) != len(novelty_thresholds):
        raise ValueError(
            "Number of novelty columns must match the number of novelty thresholds."
        )

    # Load molecules
    molecules = pd.read_csv(data_path)

    print(f"Number of molecules = {len(molecules):,}")

    # Filter by rollout
    if max_rollout is not None:
        molecules = molecules[molecules["rollout_num"] <= max_rollout]

        print(
            f"Number of molecules through rollout {max_rollout:,} = {len(molecules):,}"
        )

    num_molecules = len(molecules)

    # Filter to hits (i.e., molecules matching/exceeding all score thresholds)
    hits = molecules.query(
        " & ".join(
            f"`{score_column}` {score_comparator}"
            for score_column, score_comparator in zip(score_columns, score_comparators)
        )
    )

    num_hits = len(hits)
    percent_hits = 100 * num_hits / num_molecules
    print(f"Number of hits = {num_hits:,} ({percent_hits:.2f}%)")

    # Filter to novel hits
    novel_hits = hits.query(
        " & ".join(
            f"`{novelty_column}` <= {novelty_threshold}"
            for novelty_column, novelty_threshold in zip(
                novelty_columns, novelty_thresholds
            )
        )
    )

    num_novel_hits = len(novel_hits)
    percent_novel_hits = 100 * num_novel_hits / num_molecules
    print(f"Number of novel hits = {num_novel_hits:,} ({percent_novel_hits:.2f}%)")

    # If novel hits, compute diversity and select independent set
    if len(novel_hits) > 0:
        # Calculate pairwise similarities
        novel_hit_molecules = novel_hits[smiles_column].tolist()
        pairwise_similarities = compute_pairwise_tanimoto_similarities(
            novel_hit_molecules
        )

        # Compute average maximum pairwise similarity
        average_maximum_similarity = compute_average_maximum_similarity(
            similarities=pairwise_similarities
        )

        print(
            f"Novel hits average maximum similarity = {average_maximum_similarity:.2f}"
        )

        # Compute maximum independent set
        maximum_independent_set_indices = get_approximate_maximum_independent_set(
            similarities=pairwise_similarities, threshold=similarity_threshold
        )
        maximum_independent_set_size = len(maximum_independent_set_indices)

        print(
            f"Novel hits maximum independent set size = {maximum_independent_set_size:,}"
        )

        # Select molecules from maximum independent set
        selected = novel_hits.iloc[sorted(maximum_independent_set_indices)]

    # Otherwise, skip diversity analysis
    else:
        average_maximum_similarity = float("nan")
        maximum_independent_set_size = 0
        selected = pd.DataFrame()

    # Sort and select molecules
    if len(selected) > 0:
        # Sort molecules
        if sort_column is not None:
            print(f"Sorting molecules by {sort_column}")
            selected = selected.sort_values(sort_column, ascending=not descending)

        # Select molecules
        if select_num is not None and len(selected) > select_num:
            print(f"Selecting {select_num:,} molecules")
            selected = selected.iloc[:select_num]

    # Save selected molecules
    save_molecules_path.parent.mkdir(parents=True, exist_ok=True)
    selected.to_csv(save_molecules_path, index=False)

    # Save analysis
    if save_analysis_path is not None:
        save_analysis_path.parent.mkdir(parents=True, exist_ok=True)

        pd.DataFrame(
            [
                {
                    "num_molecules": num_molecules,
                    "num_hits": num_hits,
                    "percent_hits": percent_hits,
                    "num_novel_hits": num_novel_hits,
                    "percent_novel_hits": percent_novel_hits,
                    "novel_hits_average_maximum_similarity": average_maximum_similarity,
                    "novel_hits_maximum_independent_set_size": maximum_independent_set_size,
                }
            ]
        ).to_csv(save_analysis_path, index=False)


if __name__ == "__main__":
    from tap import tapify

    tapify(select_molecules)
