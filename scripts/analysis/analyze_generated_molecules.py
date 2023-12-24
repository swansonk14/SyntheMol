"""Analyzes a set of generated molecules for hits, novelty, and diversity."""
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
    :param threshold: Similarity threshold.
    :param num_tries: Number of times to run the approximation (i.e., number of maximal independent sets to try).
    :return: The approximate maximum independent set (set of indices of molecules).
    """
    # Create similarities graph
    similarities_matrix = similarities >= threshold
    similarities_graph = nx.from_numpy_array(similarities_matrix)

    # Compute approximate maximum independent set
    max_independent_set = set()
    for seed in range(num_tries):
        independent_set = nx.maximal_independent_set(similarities_graph, seed=seed)
        if len(independent_set) > len(max_independent_set):
            max_independent_set = independent_set

    return max_independent_set


def analyze_generated_molecules(
    data_path: Path,
    save_analysis_path: Path,
    save_molecules_path: Path,
    smiles_column: str = SMILES_COL,
    score_columns: tuple[str, ...] = (SCORE_COL,),
    novelty_columns: tuple[str, ...] = (
        "train_hits_tversky_nearest_neighbor_similarity",
        "chembl_tversky_nearest_neighbor_similarity",
    ),
    score_thresholds: tuple[float, ...] = (0.5,),
    novelty_thresholds: tuple[float, ...] = (0.6, 0.6),
    similarity_threshold: float = 0.6,
    max_rollout: int | None = None,
) -> None:
    """Analyzes a set of generated molecules for hits, novelty, and diversity.

    :param data_path: Path to CSV file containing generated molecules.
    :param save_analysis_path: Path to CSV file where the analysis will be saved.
    :param save_molecules_path: Path to a CSV file where the selected molecules will be saved.
    :param smiles_column: Name of the column containing SMILES.
    :param score_columns: Name of the columns containing scores.
    :param novelty_columns: Name of the columns containing similarity scores compared to known hits (for novelty).
    :param score_thresholds: Thresholds to use for calculating the hits and similarity (one per score column).
    :param novelty_thresholds: Thresholds to use for filtering by novelty (one per novelty column).
    :param similarity_threshold: Threshold to use for calculating the maximum independent set (diverse molecules).
    :param max_rollout: Maximum rollout number to include in the analysis.
    """
    # Create save directories
    save_analysis_path.parent.mkdir(parents=True, exist_ok=True)
    save_molecules_path.parent.mkdir(parents=True, exist_ok=True)

    # Load generated molecules
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
            f"`{score_column}` >= {score_threshold}"
            for score_column, score_threshold in zip(score_columns, score_thresholds)
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

    # Create DataFrame with results
    results = pd.DataFrame(
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
    )

    # Save results
    results.to_csv(save_analysis_path, index=False)
    selected.to_csv(save_molecules_path, index=False)


if __name__ == "__main__":
    from tap import tapify

    tapify(analyze_generated_molecules)
