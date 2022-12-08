"""Plots each compound's inhibition value vs prediction score across concentrations."""
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from tap import Tap

from constants import PREDICTION_SETS


class Args(Tap):
    data_path: Path  # Path to Excel file containing inhibition data.
    sheet_name: str  # Name of Excel sheet containing inhibition data.
    save_dir: Path  # Path to directory where plots will be saved.

    preds_sheet_name: str = 'Hits'  # Name of Excel sheet containing prediction data.
    preds_smiles_column: str = 'smiles'  # Column name containing SMILES strings in the prediction data.
    preds_score_column: str = 'score'  # Column name containing prediction scores in the prediction data.
    preds_compound_id_column: str = 'compound_id'  # Column name containing compound IDs in the prediction data.
    preds_prediction_set_column: str = 'prediction_set'  # Column name containing prediction sets in the prediction data.

    def process_args(self) -> None:
        self.save_dir.mkdir(parents=True, exist_ok=True)


def plot_inhibition_vs_score(args: Args) -> None:
    """Plots each compound's inhibition value vs prediction score across concentrations."""
    # Load data
    inhibition = pd.read_excel(args.data_path, sheet_name=args.sheet_name)
    preds = pd.read_excel(args.data_path, sheet_name=args.preds_sheet_name)

    # Compute the average of every pair of adjacent rows (i.e., mean of two replicates)
    inhibition = inhibition.groupby(inhibition.index // 2).mean().reset_index(drop=True)

    # Ensure preds are in order starting at 1
    assert list(range(1, len(preds) + 1)) == list(preds[args.preds_compound_id_column])

    # For each prediction set and for each concentration, plot the inhibition value vs prediction score
    for prediction_set in PREDICTION_SETS:
        # Create save directory for prediction set
        save_dir = args.save_dir / prediction_set
        save_dir.mkdir(parents=True, exist_ok=True)

        # Get mask for prediction set
        prediction_set_mask = preds[args.preds_prediction_set_column] == prediction_set

        for concentration in inhibition.columns:
            # Compute spearman correlation between prediction scores and inhibition values
            spearman_corr = inhibition[concentration][prediction_set_mask].corr(
                preds[args.preds_score_column][prediction_set_mask],
                method='spearman'
            )

            # Plot inhibition vs score
            plt.clf()
            plt.scatter(preds[args.preds_score_column][prediction_set_mask],
                        inhibition[concentration][prediction_set_mask])
            plt.xlabel('Prediction score')
            plt.ylabel('Mean Inhibition')
            plt.text(0.98, 0.98, f'Spearman corr. = {spearman_corr:.3f}',
                     horizontalalignment='right', verticalalignment='top',
                     transform=plt.gca().transAxes)
            plt.title(f'{prediction_set} Inhibition vs Prediction Score at {concentration} ug/mL')
            plt.savefig(save_dir / f'inhibition_vs_score_{concentration}.pdf', bbox_inches='tight')


if __name__ == '__main__':
    plot_inhibition_vs_score(Args().parse_args())
