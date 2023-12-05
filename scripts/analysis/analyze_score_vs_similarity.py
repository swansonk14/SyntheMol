"""Analyzes the score and similarity of a set of generated molecules."""
from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd
from chemfunc import compute_pairwise_tanimoto_similarities
from tqdm import tqdm

from synthemol.constants import SCORE_COL, SMILES_COL


def compute_average_maximum_similarity(similarities: np.ndarray) -> float:
    """Computes the average maximum similarity between a set of molecules.

    :param similarities: A 2D numpy array of pairwise molecule similarities.
    :return: Average maximum similarity.
    """
    np.fill_diagonal(similarities, -np.inf)
    max_similarities = np.max(similarities, axis=1)
    average_max_similarity = float(np.mean(max_similarities))

    return average_max_similarity


def compute_approximate_maximum_independent_set(
    similarities: np.ndarray, threshold: float, num_tries: int = 10
) -> int:
    """Computes the approximate size of the maximum independent set of molecules within a similarity threshold.

    Note: This does not provide a formal approximation since there are no theoretical guarantees on using
          maximal independent sets to approximate maximum independent sets.

    :param similarities: A 2D numpy array of pairwise molecule similarities.
    :param threshold: Similarity threshold.
    :param num_tries: Number of times to run the approximation (i.e., number of maximal independent sets to try).
    :return: The approximate size of the maximum independent set.
    """
    # Create similarities graph
    similarities_matrix = similarities >= threshold
    similarities_graph = nx.from_numpy_array(similarities_matrix)

    # Compute approximate maximum independent set
    max_independent_set_size = 0
    for seed in range(num_tries):
        independent_set = nx.maximal_independent_set(similarities_graph, seed=seed)
        max_independent_set_size = max(max_independent_set_size, len(independent_set))

    return max_independent_set_size


def analyze_score_vs_similarity(
    data_path: Path,
    save_path: Path,
    smiles_column: str = SMILES_COL,
    score_column: str = SCORE_COL,
    novelty_column: str = "train_hits_tversky_nearest_neighbor_similarity",
    novelty_threshold: float = 0.5,
    similarity_threshold: float = 0.5,
    score_thresholds: tuple[float, ...] = (0.9, 0.95, 0.98, 0.99),
    max_rollout: int | None = None,
) -> None:
    """Analyzes the score and similarity of a set of generated molecules.

    :param data_path: Path to CSV file containing generated molecules.
    :param save_path: Path to CSV file where the analysis will be saved.
    :param smiles_column: Name of the column containing SMILES.
    :param score_column: Name of the column containing scores.
    :param novelty_column: Name of the column containing similarity scores compared to known hits (for novelty).
    :param novelty_threshold: Threshold to use for filtering by novelty.
    :param similarity_threshold: Threshold to use for calculating the maximum independent set.
    :param score_thresholds: Thresholds to use for calculating the hits and similarity.
    :param max_rollout: Maximum rollout number to include in the analysis.
    """
    # Load data
    data = pd.read_csv(data_path)

    print(f"Number of molecules = {len(data):,}")

    # Filter by rollout
    if max_rollout is not None:
        data = data[data["rollout_num"] <= max_rollout]
        print(f"Number of molecules after filtering by rollout = {len(data):,}")

    # For each threshold, calculate the percent of hits and similarity among the hits
    results = []
    for score_threshold in tqdm(score_thresholds, desc="score thresholds"):
        # Select molecules above threshold as hits
        hits = data[data[score_column] >= score_threshold]

        if len(hits) == 0:
            results.append(
                {
                    "score_threshold": score_threshold,
                    "num_hits": 0,
                    "percent_hits": 0,
                    "similarity": 0,
                    "max_independent_set": 0,
                    "num_hits_novel": 0,
                    "percent_hits_novel": 0,
                    "similarity_novel": 0,
                    "max_independent_set_novel": 0,
                }
            )
            continue

        # Calculate statistics
        hit_molecules = hits[smiles_column].tolist()
        pairwise_similarities = compute_pairwise_tanimoto_similarities(hit_molecules)

        result = {
            "score_threshold": score_threshold,
            "num_hits": len(hits),
            "percent_hits": len(hits) / len(data),
            "similarity": compute_average_maximum_similarity(
                similarities=pairwise_similarities
            ),
            "max_independent_set": compute_approximate_maximum_independent_set(
                similarities=pairwise_similarities, threshold=similarity_threshold
            ),
        }

        # Compute same results for novel molecules
        hits_novel = hits[hits[novelty_column] <= novelty_threshold]

        if len(hits_novel) == 0:
            result["num_hits_novel"] = 0
            result["percent_hits_novel"] = 0
            result["similarity_novel"] = 0
            result["max_independent_set_novel"] = 0
            results.append(result)
            continue

        hit_molecules_novel = hits_novel[smiles_column].tolist()
        pairwise_similarities_novel = compute_pairwise_tanimoto_similarities(
            hit_molecules_novel
        )

        result["num_hits_novel"] = len(hits_novel)
        result["percent_hits_novel"] = len(hits_novel) / len(data)
        result["similarity_novel"] = compute_average_maximum_similarity(
            similarities=pairwise_similarities_novel
        )
        result[
            "max_independent_set_novel"
        ] = compute_approximate_maximum_independent_set(
            similarities=pairwise_similarities_novel, threshold=similarity_threshold
        )

        results.append(result)

    # Create DataFrame with results
    results = pd.DataFrame(results)

    # Save results
    save_path.parent.mkdir(parents=True, exist_ok=True)
    results.to_csv(save_path, index=False)


if __name__ == "__main__":
    from tap import tapify

    tapify(analyze_score_vs_similarity)
