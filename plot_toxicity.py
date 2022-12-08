"""Plot predicted toxicity for each prediction set of molecules."""
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from tap import Tap


class Args(Tap):
    data_path: Path  # Path to Excel file containing molecule toxicity.
    save_dir: Path  # Path to directory where the plots will be saved.

    toxicity_sheet_name: str = 'Hits'  # Name of Excel sheet containing toxicity data.
    toxicity_columns: list[str] = ['FDA_APPROVED', 'CT_TOX']  # Column names containing toxicity data.
    prediction_set_column: str = 'prediction_set'  # Column name containing prediction set in the synergy data.

    def process_args(self) -> None:
        self.save_dir.mkdir(parents=True, exist_ok=True)


def plot_synergy_vs_score(args: Args) -> None:
    """Plot a molecule's synergy vs prediction score."""
    # Load data
    toxicity = pd.read_excel(args.data_path, sheet_name=args.toxicity_sheet_name)

    print(f'Size of toxicity data = {len(toxicity):,}')

    # Get prediction sets
    prediction_sets = sorted(toxicity[args.prediction_set_column].unique())

    # Plot toxicity for each toxicity metric
    for toxicity_column in args.toxicity_columns:
        plt.clf()

        # Sort toxicity data by prediction set and then toxicity score within prediction set
        toxicity = toxicity.sort_values([args.prediction_set_column, toxicity_column]).reset_index(drop=True)

        for prediction_set in prediction_sets:
            prediction_set_toxicity = toxicity[toxicity[args.prediction_set_column] == prediction_set]
            plt.scatter(prediction_set_toxicity.index, prediction_set_toxicity[toxicity_column],
                        label=prediction_set)

        plt.xlabel('Molecule Index')
        plt.ylabel(f'Predicted {toxicity_column}')
        plt.title(f'Predicted {toxicity_column} by Prediction Set')
        plt.legend()
        plt.savefig(args.save_dir / f'{toxicity_column}.pdf', bbox_inches='tight')


if __name__ == '__main__':
    plot_synergy_vs_score(Args().parse_args())
