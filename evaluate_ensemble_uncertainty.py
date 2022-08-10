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
        # TODO: bigger figure size
        plt.savefig(args.save_dir / f'{model_name}_aucs.pdf', bbox_inches='tight')

        # Plot hit ratios
        plt.clf()
        plt.plot(np.arange(len(hit_ratios)), hit_ratios, color='blue', label='Hit ratio')
        plt.hlines(y=hit_ratio, xmin=0, xmax=len(hit_ratios), linestyles='--', color='blue', label='Full data hit ratio')

        for top_k, top_k_hit_ratios in top_k_hit_ratio_dict.items():
            p = plt.plot(np.arange(len(top_k_hit_ratios)), top_k_hit_ratios, label=f'Top {top_k} hit ratio')
            plt.hlines(y=ensemble_top_k_hit_ratios[top_k], xmin=0, xmax=len(top_k_hit_ratios), linestyles='--', color=p[0].get_color(), label=f'Full data top {top_k} hit ratio')

        plt.xlabel('Number of molecules removed by uncertainty')
        plt.ylabel('Hit ratio')
        plt.title(f'{model_name}: hit ratios by uncertainty for {num_models} model ensemble')
        plt.legend(loc='upper left')
        plt.savefig(args.save_dir / f'{model_name}_hit_ratios.pdf', bbox_inches='tight')





# def evaluate_ensemble_confidence(args: Args) -> None:
#     """Evaluate the confidence of an ensemble for selecting molecules."""
#     # Load data
#     chemprop_data = pd.read_csv(args.chemprop_data_path)
#     rf_data = pd.read_csv(args.rf_data_path)
#
#     chemprop_preds_columns = [column for column in chemprop_data.columns if column.endswith(args.preds_column_suffix)]
#     rf_preds_columns = [column for column in rf_data.columns if column.endswith(args.preds_column_suffix)]
#
#     chemprop_preds = chemprop_data[chemprop_preds_columns].to_numpy()
#     rf_preds = rf_data[rf_preds_columns].to_numpy()
#
#     chemprop_percentiles = chemprop_preds.argsort(axis=0).argsort(axis=0) / (len(chemprop_data) - 1)
#     rf_percentiles = rf_preds.argsort(axis=0).argsort(axis=0) / (len(rf_data) - 1)
#
#     model_num = 0
#
#     hit_mask = chemprop_data[args.activity_column] == 1
#
#     # agree = np.abs(chemprop_percentiles[:, model_num] - rf_percentiles[:, model_num]) <= 0.05
#     agree = np.abs(chemprop_percentiles[:, model_num] - chemprop_percentiles[:, model_num + 2]) <= 0.05
#
#     print(sum(chemprop_data[args.activity_column] == 1) / len(chemprop_data))
#     print(sum(chemprop_data[args.activity_column][agree] == 1) / sum(agree))
#
#     print('Chemprop')
#     print(f'ROC-AUC = {roc_auc_score(chemprop_data[args.activity_column][agree], chemprop_percentiles[:, model_num][agree]):.3f}')
#     print(f'PRC-AUC = {average_precision_score(chemprop_data[args.activity_column][agree], chemprop_percentiles[:, model_num][agree]):.3f}')
#     print()
#
#
#     print('RF')
#     print(f'ROC-AUC = {roc_auc_score(chemprop_data[args.activity_column][agree], rf_percentiles[:, model_num][agree]):.3f}')
#     print(f'PRC-AUC = {average_precision_score(chemprop_data[args.activity_column][agree], rf_percentiles[:, model_num][agree]):.3f}')
#     print()
#
#     print('Ensemble')
#     print(f'ROC-AUC = {roc_auc_score(chemprop_data[args.activity_column][agree], (chemprop_percentiles[:, model_num] + rf_percentiles[:, model_num])[agree]):.3f}')
#     print(f'PRC-AUC = {average_precision_score(chemprop_data[args.activity_column][agree], (chemprop_percentiles[:, model_num] + rf_percentiles[:, model_num])[agree]):.3f}')
#     print()
#
#     exit()
#
#     import matplotlib.pyplot as plt
#     plt.clf()
#     # plt.scatter(chemprop_percentiles[:, model_num][~hit_mask], np.abs(chemprop_percentiles[:, model_num] - rf_percentiles[:, model_num])[~hit_mask], label='nonhit', color='blue')
#     # plt.scatter(chemprop_percentiles[:, model_num][hit_mask], np.abs(chemprop_percentiles[:, model_num] - rf_percentiles[:, model_num])[hit_mask], label='hit', color='red')
#     plt.scatter(chemprop_percentiles[:, model_num][~hit_mask], rf_percentiles[:, model_num][~hit_mask], label='nonhit', color='blue')
#     plt.scatter(chemprop_percentiles[:, model_num][hit_mask], rf_percentiles[:, model_num][hit_mask], label='hit', color='red')
#     plt.xlabel('Chemprop Percentile')
#     # plt.ylabel('Absolute Value Percentile Difference')
#     plt.ylabel('RF Percentile')
#     plt.legend()
#     plt.show()
#
#     exit()
#
#
#     # Get preds
#     preds_columns = [column for column in data.columns if column.endswith(args.preds_column_suffix)]
#     all_preds = data[preds_columns].to_numpy()
#
#     # Convert preds scores into percentiles for comparability
#     ranking = all_preds.argsort(axis=0).argsort(axis=0)
#     all_preds = ranking / (len(data) - 1)  # percentiles
#
#     # Compute variance
#     preds_std = np.std(all_preds, axis=1)
#
#     # Use first model's predictions
#     pred_column = preds_columns[0]
#
#     # Compute performance of first model as a function of std percentile among ensemble
#     # percentiles = np.arange(0, 101, 10)
#     percentiles = [0, 90, 95, 100]
#     percentile_stds = np.percentile(preds_std, percentiles)
#
#     hit_mask = data[args.activity_column] == 1
#     hit_stds = preds_std[hit_mask]
#     nonhit_stds = preds_std[~hit_mask]
#     print(np.mean(hit_stds))
#     print(np.mean(nonhit_stds))
#     print()
#
#     preds = all_preds.mean(axis=1)
#     print(f'ROC-AUC = {roc_auc_score(data[args.activity_column], preds):.3f}')
#     print(f'PRC-AUC = {average_precision_score(data[args.activity_column], preds):.3f}')
#     print(f'Positive enrichment = {sum(data[args.activity_column] == 1) / len(data):.3f}')
#     print()
#
#     high_score_mask = preds >= 0.98
#     # for std in np.arange(0.35, 0.0, -0.05):
#     for std in np.arange(0.05, 0.0, -0.005):
#         std_mask = preds_std <= std
#         mask = high_score_mask & std_mask
#         print(f'std = {std:.3f}')
#         activities = data[args.activity_column][mask]
#         print(activities.value_counts())
#         print(f'ROC-AUC = {roc_auc_score(data[args.activity_column][mask], preds[mask]):.3f}')
#         print(f'PRC-AUC = {average_precision_score(data[args.activity_column][mask], preds[mask]):.3f}')
#         print(f'Positive enrichment = {sum(activities == 1) / len(activities):.3f}')
#         print()
#
#     # exit()
#
#     import matplotlib.pyplot as plt
#     for model_num in range(len(preds_columns)):
#         plt.clf()
#         # plt.scatter(all_preds.mean(axis=1), preds_std, color='blue', s=5)
#         # plt.scatter(all_preds[:, model_num], preds_std, color='blue', s=5)
#         plt.scatter(all_preds.mean(axis=1)[~hit_mask], nonhit_stds, label='nonhit', color='blue', s=5)
#         plt.scatter(all_preds.mean(axis=1)[hit_mask], hit_stds, label='hit', color='red', s=5)
#         # plt.hist(nonhit_stds, density=True, bins=100, label='nonhit', alpha=0.5)
#         # plt.hist(hit_stds, density=True, bins=100, label='hit', alpha=0.5)
#         # plt.xlabel('Prediction Score')
#         plt.xlabel('Ranking')
#         plt.ylabel('Standard Deviation')
#         plt.legend()
#         plt.show()
#         exit()
#     exit()
#
#     # for model_num in range(len(preds_columns)):
#     #     print(roc_auc_score(data[args.activity_column], all_preds[:, model_num]))
#     #     print(average_precision_score(data[args.activity_column], all_preds[:, model_num]))
#     #     print()
#     #
#     # print(roc_auc_score(data[args.activity_column], preds_std))
#     # print(average_precision_score(data[args.activity_column], preds_std))
#
#
#     for i in range(len(percentiles) - 1):
#         print(f'Std percentiles {percentiles[i]} to {percentiles[i + 1]}')
#
#         selected_data = data[(preds_std > percentile_stds[i]) & (preds_std <= percentile_stds[i + 1])]
#         print(selected_data[args.activity_column].value_counts())
#
#         activity_classes = selected_data[args.activity_column].unique()
#         if len(activity_classes) == 1:
#             continue
#
#         roc_auc = roc_auc_score(selected_data[args.activity_column], selected_data[pred_column])
#         prc_auc = average_precision_score(selected_data[args.activity_column], selected_data[pred_column])
#
#         print(f'ROC-AUC = {roc_auc:.3}')
#         print(f'PRC-AUC = {prc_auc:.3}')
#         print()


if __name__ == '__main__':
    evaluate_ensemble_uncertainty(Args().parse_args())
