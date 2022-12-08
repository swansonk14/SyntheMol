"""Plot a molecule's synergy vs prediction score."""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from rdkit import Chem
from tap import Tap

from constants import PREDICTION_SETS


class Args(Tap):
    data_path: Path  # Path to Excel file containing molecule synergy.
    save_dir: Path  # Path to directory where the plots will be saved.

    synergy_sheet_name: str = 'Synergy'  # Name of Excel sheet containing synergy data.
    synergy_smiles_column: str = 'Smile'  # Column name containing SMILES strings in the synergy data.
    synergy_prediction_set_column: str = 'prediction_set'  # Column name containing prediction set in the synergy data.
    synergy_column: str = 'FICi(16)'  # Name of the column containing synergy.
    synergy_threshold: float = 0.5  # Threshold for synergy values to be plotted.
    preds_sheet_name: str = 'Hits'  # Name of Excel sheet containing prediction data.
    preds_smiles_column: str = 'smiles'  # Column name containing SMILES strings in the prediction data.
    preds_score_column: str = 'score'  # Column name containing prediction scores in the prediction data.

    def process_args(self) -> None:
        self.save_dir.mkdir(parents=True, exist_ok=True)


def plot_synergy_vs_score(args: Args) -> None:
    """Plot a molecule's synergy vs prediction score."""
    # Load data
    synergy = pd.read_excel(args.data_path, sheet_name=args.synergy_sheet_name)
    preds = pd.read_excel(args.data_path, sheet_name=args.preds_sheet_name)

    print(f'Size of synergy data = {len(synergy):,}')
    print(f'Size of prediction data = {len(preds):,}')

    # Canonicalize synergy SMILES
    synergy[args.preds_smiles_column] = [Chem.MolToSmiles(Chem.MolFromSmiles(smiles)) for smiles in synergy['Smile']]

    # Get all tested SMILES
    tested_smiles = set(synergy[args.preds_smiles_column])

    print(f'Number of unique molecules = {len(tested_smiles):,}')

    # Map SMILES to prediction scores
    smiles_to_score = {
        smiles: score for smiles, score in zip(preds[args.preds_smiles_column], preds[args.preds_score_column])
    }

    # Set prediction scores in synergy data
    synergy[args.preds_score_column] = [smiles_to_score[smiles] for smiles in synergy[args.preds_smiles_column]]

    # Sort synergy data by prediction score
    synergy = synergy.sort_values(args.preds_score_column).reset_index(drop=True)

    # Plot synergy vs score
    for prediction_set in PREDICTION_SETS:
        # Limit synergy data to prediction set
        prediction_set_synergy = synergy[synergy[args.synergy_prediction_set_column] == prediction_set].reset_index(drop=True)

        print(f'\n{prediction_set}')
        print(f'Number of molecules = {len(prediction_set_synergy):,}')

        # Create mask for molecules with activity
        activity_mask = np.array([
            isinstance(syn, float) and not np.isnan(syn)
            for syn in prediction_set_synergy[args.synergy_column]
        ])

        print(f'Number of molecules with activity = {sum(activity_mask):,}')

        # Create mask for molecules with synergy
        synergy_mask = np.array([
            isinstance(syn, float) and not np.isnan(syn) and syn <= args.synergy_threshold
            for syn in prediction_set_synergy[args.synergy_column]
        ])

        print(f'Number of molecules with synergy = {sum(synergy_mask):,}')

        plt.clf()
        plt.scatter(prediction_set_synergy[~synergy_mask].index,
                    prediction_set_synergy[~synergy_mask][args.preds_score_column],
                    color='blue', label='No Activity')
        plt.scatter(prediction_set_synergy[activity_mask].index,
                    prediction_set_synergy[activity_mask][args.preds_score_column],
                    color='orange', label='Activity')
        plt.scatter(prediction_set_synergy[synergy_mask].index,
                    prediction_set_synergy[synergy_mask][args.preds_score_column],
                    color='red', label='Activity + Synergy')
        plt.xlabel('Molecule Index')
        plt.ylabel('Prediction Score')
        plt.title(f'{prediction_set} Prediction Score vs Activity/Synergy with SPR-741')
        plt.legend()

        plt.savefig(args.save_dir / f'{prediction_set}_synergy.pdf', bbox_inches='tight')


if __name__ == '__main__':
    plot_synergy_vs_score(Args().parse_args())
