"""Adds labels to test predictions from Chemprop (and corrects SMILES format)."""
from pathlib import Path

import pandas as pd
from tqdm import tqdm

from synthemol.constants import SMILES_COL


def add_labels_to_test_preds(
    true_path: Path,
    preds_path: Path,
    value_column: str,
    smiles_column: str = SMILES_COL,
) -> None:
    """Adds labels to test predictions from Chemprop (and corrects SMILES format).

    :param true_path: Path to a CSV file containing the true labels.
    :param preds_path: Path to a CSV file containing the predictions or to a directory containing test_preds.csv files.
    :param value_column: Name of the column containing the predicted/true values.
    :param smiles_column: Name of the column containing SMILES.
    """
    # Load true labels
    true = pd.read_csv(true_path, index_col=smiles_column)

    # Get preds paths
    if preds_path.is_dir():
        preds_paths = sorted(preds_path.rglob("test_preds.csv"))
    else:
        preds_paths = [preds_path]

    # Add labels to each preds file
    for preds_path in tqdm(preds_paths):
        # Load preds
        preds = pd.read_csv(preds_path, index_col=smiles_column)

        # Correct SMILES format if needed (remove beginning "['" and ending "']" and fix backslashes)
        if preds.index[0].startswith("['"):
            preds.index = preds.index.str[2:-2]
            preds.index = preds.index.str.replace("\\\\", "\\")

        # Add labels
        preds[f"{value_column}_true"] = true[value_column]

        # Ensure every SMILES got a label
        if preds[f"{value_column}_true"].isna().any():
            breakpoint()
            raise ValueError(
                f"Not all SMILES in {preds_path} were found in {true_path}."
            )

        # Save preds
        preds.to_csv(preds_path)


if __name__ == "__main__":
    from tap import tapify

    tapify(add_labels_to_test_preds)
