"""Evaluate the uncertainty of an ensemble for selecting molecules."""
from pathlib import Path
from typing import Literal

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score
from tap import Tap
from tqdm import trange


class Args(Tap):
    chemprop_data_path: Path  # Path to CSV file containing chemprop ensemble predictions on a test set.
    rf_data_path: Path  # Path to CSV file containing random forest ensemble predictions on a test set.
    save_dir: Path  # Path to directory where plots will be saved.
    activity_column: str = 'activity'  # Name of the column containing true activity values.
    preds_column_suffix: str = 'preds'  # Suffix of the names of the columns containing predictions.
    smiles_column: str = 'smiles'  # Name of the column containing SMILES.
    top_ks: list[int] = [10, 25, 50, 75, 100, 150, 200, 500]  # Number of top scoring molecules to select for computing hit ratios.
    uncertainty_basis: Literal['score', 'percentile'] = 'percentile'  # The basis on which to compute the uncertainty.
    """Either the raw model score or the percentile rank of the score compared to all other scores in the dataset."""
    model_comparison: Literal['intra', 'inter'] = 'intra'  # Whether to compare between models of the same type or different types.

    def process_args(self) -> None:
        self.save_dir.mkdir(parents=True, exist_ok=True)


def evaluate_ensemble_uncertainty(args: Args) -> None:
    """Evaluate the uncertainty of an ensemble for selecting molecules."""
    # Load data
    chemprop_data = pd.read_csv(args.chemprop_data_path)
    rf_data = pd.read_csv(args.rf_data_path)

    # Ensure both data sets use same set of molecules and activities in the same order
    assert chemprop_data[args.smiles_column].equals(rf_data[args.smiles_column])
    assert chemprop_data[args.activity_column].equals(rf_data[args.activity_column])

    # Get activities
    activities = chemprop_data[args.activity_column].to_numpy()

    # Get number of molecules
    num_molecules = len(activities)

    print(f'Data size = {num_molecules:,}')

    # Get predictions
    chemprop_preds_columns = [column for column in chemprop_data.columns if column.endswith(args.preds_column_suffix)]
    chemprop_preds = chemprop_data[chemprop_preds_columns].to_numpy().transpose()  # (num_models, num_molecules)
    chemprop_ensemble_pred = chemprop_preds.mean(axis=0)

    rf_preds_columns = [column for column in rf_data.columns if column.endswith(args.preds_column_suffix)]
    rf_preds = rf_data[rf_preds_columns].to_numpy().transpose()  # (num_models, num_molecules)
    rf_ensemble_pred = rf_preds.mean(axis=0)

    # Compute hit ratio
    hit_ratio = np.mean(activities)
    print(f'Full data hit ratio = {hit_ratio:.3f}')
    print()

    # Evaluate models and uncertainty measures
    for model_name, model_preds, ensemble_pred in [('chemprop', chemprop_preds, chemprop_ensemble_pred),
                                                   ('RF', rf_preds, rf_ensemble_pred)]:
        # Compute metrics for individual models and ensembles
        model_roc_aucs = [roc_auc_score(activities, preds) for preds in model_preds]
        model_prc_aucs = [average_precision_score(activities, preds) for preds in model_preds]
        ensemble_roc_auc = roc_auc_score(activities, ensemble_pred)
        ensemble_prc_auc = average_precision_score(activities, ensemble_pred)

        ensemble_top_k_hit_ratios = {}
        for top_k in args.top_ks:
            top_k_indices = ensemble_pred.argsort()[-top_k:]
            top_k_activities = activities[top_k_indices]
            ensemble_top_k_hit_ratios[top_k] = np.mean(top_k_activities)

        num_models = len(model_preds)

        print(f'{model_name} {num_models} models')
        print(f'ROC-AUC = {np.mean(model_roc_aucs):.3f} +/- {np.std(model_roc_aucs):.3f} '
              f'vs {ensemble_roc_auc:.3f} ensemble')
        print(f'PRC-AUC = {np.mean(model_prc_aucs):.3f} +/- {np.std(model_prc_aucs):.3f} '
              f'vs {ensemble_prc_auc:.3f} ensemble')
        for top_k, ensemble_top_k_hit_ratio in ensemble_top_k_hit_ratios.items():
            print(f'Ensemble top {top_k} hit ratio = {ensemble_top_k_hit_ratio:.3f}')
        print()

        # Compute model uncertainty
        if args.model_comparison == 'intra':
            if args.uncertainty_basis == 'score':
                uncertainty = np.std(model_preds, axis=0)
            elif args.uncertainty_basis == 'percentile':
                model_percentiles = model_preds.argsort(axis=1).argsort(axis=1) / (num_molecules - 1)
                uncertainty = np.std(model_percentiles, axis=0)
            else:
                raise ValueError(f'Uncertainty basis "{args.uncertainty_basis}" is not supported.')
        elif args.model_comparison == 'inter':
            if args.uncertainty_basis == 'score':
                diffs = np.abs(chemprop_preds - rf_preds)
                uncertainty = np.mean(diffs, axis=0)
            elif args.uncertainty_basis == 'percentile':
                chemprop_percentiles = chemprop_preds.argsort(axis=1).argsort(axis=1) / (num_molecules - 1)
                rf_percentiles = rf_preds.argsort(axis=1).argsort(axis=1) / (num_molecules - 1)
                diffs = np.abs(chemprop_percentiles - rf_percentiles)
                uncertainty = np.mean(diffs, axis=0)
            else:
                raise ValueError(f'Uncertainty basis "{args.uncertainty_basis}" is not supported.')
        else:
            raise ValueError(f'Model comparison "{args.model_comparison}" is not supported.')

        uncertainty_argsort = np.argsort(-uncertainty)  # sort from highest to lowest uncertainty
        ensemble_pred_uncertainty_sorted = ensemble_pred[uncertainty_argsort]
        activities_uncertainty_sorted = activities[uncertainty_argsort]

        # Plot uncertainty vs score
        hit_mask = activities == 1
        non_hit_mask = activities == 0

        ensemble_percentiles = ensemble_pred.argsort().argsort() / (num_molecules - 1)

        plt.clf()
        plt.scatter(ensemble_percentiles[non_hit_mask], uncertainty[non_hit_mask], color='blue', label='non-hit')
        plt.scatter(ensemble_percentiles[hit_mask], uncertainty[hit_mask], color='red', label='hit')
        plt.xlabel('Ensemble Prediction Score Percentile')
        plt.ylabel(f'{args.model_comparison.title()}-Model {args.uncertainty_basis.title()} Uncertainty')
        plt.title(f'{model_name} Uncertainty vs Prediction Score Percentile')
        plt.legend()
        plt.gcf().set_size_inches(12, 9.5)
        plt.savefig(args.save_dir / f'{model_name}_uncertainty_vs_score.pdf', bbox_inches='tight')

        # Evaluate at different uncertainty levels
        roc_aucs, prc_aucs, hit_ratios = [], [], []
        top_k_hit_ratio_dict = {top_k: [] for top_k in args.top_ks}

        for i in trange(num_molecules):
            # Filter ensemble preds and activities by uncertainty level
            filtered_ensemble_pred = ensemble_pred_uncertainty_sorted[i:]
            filtered_activities = activities_uncertainty_sorted[i:]

            # Compute hit ratio
            hit_ratios.append(np.mean(filtered_activities))

            # Compute top k hit ratios
            for top_k in args.top_ks:
                if len(filtered_ensemble_pred) >= top_k:
                    filtered_top_k_indices = filtered_ensemble_pred.argsort()[-top_k:]
                    filtered_top_k_activities = filtered_activities[filtered_top_k_indices]
                    top_k_hit_ratio_dict[top_k].append(np.mean(filtered_top_k_activities))

            # Compute AUCs
            if len(set(filtered_activities)) == 1:
                break

            roc_aucs.append(roc_auc_score(filtered_activities, filtered_ensemble_pred))
            prc_aucs.append(average_precision_score(filtered_activities, filtered_ensemble_pred))

        # Plot AUCs
        plt.clf()
        plt.plot(np.arange(len(roc_aucs)), roc_aucs, color='blue', label='ROC-AUC')
        plt.plot(np.arange(len(prc_aucs)), prc_aucs, color='red', label='PRC-AUC')
        plt.hlines(y=ensemble_roc_auc, xmin=0, xmax=len(roc_aucs), linestyles='--', color='blue', label='Ensemble ROC-AUC')
        plt.hlines(y=ensemble_prc_auc, xmin=0, xmax=len(prc_aucs), linestyles='--', color='red', label='Ensemble PRC-AUC')
        plt.xlabel('Number of molecules removed by uncertainty')
        plt.ylabel('AUC')
        plt.title(f'{model_name}: AUCs by uncertainty for {num_models} model ensemble')
        plt.legend()
        plt.gcf().set_size_inches(12, 9.5)
        plt.savefig(args.save_dir / f'{model_name}_aucs.pdf', bbox_inches='tight')

        # Plot hit ratios
        plt.clf()

        for top_k, top_k_hit_ratios in top_k_hit_ratio_dict.items():
            p = plt.plot(np.arange(len(top_k_hit_ratios)), top_k_hit_ratios, label=f'Top {top_k} hit ratio')
            plt.hlines(y=ensemble_top_k_hit_ratios[top_k], xmin=0, xmax=len(top_k_hit_ratios), linestyles='--', color=p[0].get_color())

        plt.plot(np.arange(len(hit_ratios)), hit_ratios, color='blue', label='Hit ratio')
        plt.hlines(y=hit_ratio, xmin=0, xmax=len(hit_ratios), linestyles='--', color='blue')

        plt.xlabel('Number of molecules removed by uncertainty')
        plt.ylabel('Hit ratio')
        plt.title(f'{model_name}: hit ratios by uncertainty for {num_models} model ensemble')
        plt.legend(loc='upper left')
        plt.gcf().set_size_inches(12, 9.5)
        plt.savefig(args.save_dir / f'{model_name}_hit_ratios.pdf', bbox_inches='tight')


if __name__ == '__main__':
    evaluate_ensemble_uncertainty(Args().parse_args())
