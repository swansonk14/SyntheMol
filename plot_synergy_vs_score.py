"""Plot a molecule's synergy vs prediction score."""
from pathlib import Path
from typing import Literal

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from rdkit import Chem
from tap import Tap


class Args(Tap):
    synergy_path: Path  # Path to Excel file containing molecule synergy.
    preds_path: Path  # Path to CSV file containing MCTS generated molecules with prediction scores.
    prediction_set: Literal['rf', 'chemprop', 'chemprop_rdkit', 'random']  # Name of the prediction set whose molecules will be plotted.
    save_path: Path  # Path to PDF/PNG file where the plot will be saved.
    synergy_smiles_column: str = 'Smile'  # Column name containing SMILES strings in the synergy data.
    synergy_column: str = 'FICi(16)'  # Name of the column containing synergy.
    synergy_threshold: float = 0.5  # Threshold for synergy values to be plotted.
    preds_smiles_column: str = 'smiles'  # Column name containing SMILES strings in the prediction data.
    preds_score_column: str = 'score'  # Column name containing prediction scores in the prediction data.

    def process_args(self) -> None:
        self.save_path.parent.mkdir(parents=True, exist_ok=True)


def plot_synergy_vs_score(args: Args) -> None:
    """Plot a molecule's synergy vs prediction score."""
    # Load data
    synergy = pd.read_excel(args.synergy_path)
    preds = pd.read_csv(args.preds_path)

    print(f'Size of synergy data = {len(synergy):,}')
    print(f'Size of prediction data = {len(preds):,}')

    # Canonicalize synergy SMILES
    synergy[args.preds_smiles_column] = [Chem.MolToSmiles(Chem.MolFromSmiles(smiles)) for smiles in synergy['Smile']]

    # Limit synergy data to molecules in the prediction set
    synergy = synergy[synergy[args.preds_smiles_column].isin(preds[args.preds_smiles_column])]

    # Get all tested SMILES
    tested_smiles = set(synergy[args.preds_smiles_column])

    print(f'Number of unique tested molecules = {len(tested_smiles):,}')

    # Map SMILES to prediction scores
    smiles_to_score = {
        smiles: score for smiles, score in zip(preds[args.preds_smiles_column], preds[args.preds_score_column])
    }

    # Set prediction scores in synergy data
    synergy[args.preds_score_column] = [smiles_to_score[smiles] for smiles in synergy[args.preds_smiles_column]]

    # Sort synergy data by prediction score
    synergy = synergy.sort_values(args.preds_score_column).reset_index(drop=True)

    # Create mask for molecules with activity
    activity_mask = np.array([
        isinstance(syn, float) and not np.isnan(syn)
        for syn in synergy[args.synergy_column]
    ])

    print(f'Number of molecules with activity = {sum(activity_mask):,}')

    # Create mask for molecules with synergy
    synergy_mask = np.array([
        isinstance(syn, float) and not np.isnan(syn) and syn <= args.synergy_threshold
        for syn in synergy[args.synergy_column]
    ])

    print(f'Number of molecules with synergy = {sum(synergy_mask):,}')

    # Plot synergy vs score
    plt.clf()
    plt.scatter(synergy[~synergy_mask].index, synergy[~synergy_mask][args.preds_score_column],
                color='blue', label='No Activity')
    plt.scatter(synergy[activity_mask].index, synergy[activity_mask][args.preds_score_column],
                color='orange', label='Activity')
    plt.scatter(synergy[synergy_mask].index, synergy[synergy_mask][args.preds_score_column],
                color='red', label='Activity + Synergy')
    plt.xlabel('Molecule Index')
    plt.ylabel('Prediction Score')
    plt.title(f'{args.prediction_set} Prediction Score vs Activity/Synergy with SPR-741')
    plt.legend()

    args.save_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(args.save_path, bbox_inches='tight')


if __name__ == '__main__':
    plot_synergy_vs_score(Args().parse_args())
