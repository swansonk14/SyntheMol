"""Analyzes the score and similarity of a set of generated molecules."""
from pathlib import Path

import numpy as np
import pandas as pd
from chemfunc import compute_pairwise_tanimoto_similarities
from tqdm import tqdm

from synthemol.constants import SCORE_COL, SMILES_COL


def compute_average_maximum_similarity(molecules: list[str]) -> float:
    """Computes the average maximum similarity between a set of molecules.

    :param molecules: List of SMILES strings.
    :return: Average maximum similarity.
    """
    pairwise_similarities = compute_pairwise_tanimoto_similarities(molecules)
    np.fill_diagonal(pairwise_similarities, -np.inf)
    max_similarities = np.max(pairwise_similarities, axis=1)
    average_max_similarity = float(np.mean(max_similarities))

    return average_max_similarity


def analyze_score_vs_similarity(
    data_path: Path,
    save_path: Path,
    smiles_column: str = SMILES_COL,
    score_column: str = SCORE_COL,
    novelty_column: str = "train_hits_tversky_nearest_neighbor_similarity",
    novelty_threshold: float = 0.5,
    score_thresholds: tuple[float, ...] = (0.5, 0.75, 0.9, 0.95, 0.98, 0.99),
) -> None:
    """Analyzes the score and similarity of a set of generated molecules.

    :param data_path: Path to CSV file containing generated molecules.
    :param save_path: Path to CSV file where the analysis will be saved.
    :param smiles_column: Name of the column containing SMILES.
    :param score_column: Name of the column containing scores.
    :param novelty_column: Name of the column containing similarity scores compared to known hits (for novelty).
    :param novelty_threshold: Threshold to use for filtering by novelty.
    :param score_thresholds: Thresholds to use for calculating the hits and similarity.
    """
    # Load data
    data = pd.read_csv(data_path)

    print(f"Number of molecules = {len(data):,}")

    # Mark novelty
    data["novel"] = data[novelty_column] <= novelty_threshold

    # For each threshold, calculate the percent of hits and similarity among the hits
    num_hits = []
    percent_hits = []
    similarity = []
    num_hits_novel = []
    percent_hits_novel = []
    similarity_novel = []
    for score_threshold in tqdm(score_thresholds, desc="score thresholds"):
        # Select molecules above threshold as hits
        hits = data[data[score_column] >= score_threshold]

        if len(hits) == 0:
            num_hits.append(0)
            percent_hits.append(0)
            similarity.append(0)
            num_hits_novel.append(0)
            percent_hits_novel.append(0)
            similarity_novel.append(0)
            continue

        # Calculate statistics
        num_hits.append(len(hits))
        percent_hits.append(len(hits) / len(data))
        hit_molecules = hits[smiles_column].tolist()
        similarity.append(compute_average_maximum_similarity(hit_molecules))

        # Compute same results for novel molecules
        hits_novel = hits[hits["novel"]]

        if len(hits_novel) == 0:
            num_hits_novel.append(0)
            percent_hits_novel.append(0)
            similarity_novel.append(0)
            continue

        num_hits_novel.append(len(hits_novel))
        percent_hits_novel.append(len(hits_novel) / len(data))
        hit_molecules_novel = hits_novel[smiles_column].tolist()
        similarity_novel.append(compute_average_maximum_similarity(hit_molecules_novel))

    # Create DataFrame with results
    results = pd.DataFrame(
        {
            "score_threshold": score_thresholds,
            "num_hits": num_hits,
            "percent_hits": percent_hits,
            "similarity": similarity,
            "num_hits_novel": num_hits_novel,
            "percent_hits_novel": percent_hits_novel,
            "similarity_novel": similarity_novel,
        }
    )

    # Save results
    save_path.parent.mkdir(parents=True, exist_ok=True)
    results.to_csv(save_path, index=False)


if __name__ == "__main__":
    from tap import tapify

    tapify(analyze_score_vs_similarity)
