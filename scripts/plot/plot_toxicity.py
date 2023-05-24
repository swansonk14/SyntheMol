"""Plot predicted toxicity for a set of generated molecules compared to predictions on the test set of molecules."""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import percentileofscore
from sklearn.metrics import precision_recall_curve


def plot_toxicity(
        test_dir: Path,
        generated_path: Path,
        save_dir: Path,
        test_tox_column: str = 'CT_TOX',
        test_pred_column: str = 'prediction',
        generated_pred_column: str = 'CT_TOX'
) -> None:
    """Plot predicted toxicity for a set of generated molecules compared to predictions on the test set of molecules.

    :param test_dir: Path to a directory containing test set predictions from cross-validation.
    :param generated_path: Path to a CSV file containing predictions on a set of molecules.
    :param save_dir: Path to directory where the plot will be saved.
    :param test_tox_column: Name of the column in test_dir files containing toxicity labels.
    :param test_pred_column: Name of the column in test_dir files containing toxicity predictions.
    :param generated_pred_column: Name of the column in generated_path containing toxicity predictions.
    """
    # Load test preds
    test = []
    while True:
        test_preds_path = test_dir / f'model_{len(test)}_test_preds.csv'
        if not test_preds_path.exists():
            break
        test.append(pd.read_csv(test_preds_path))

    print(f'Found {len(test):,} test set predictions')

    test = pd.concat(test)

    print(f'Size of test set predictions = {len(test):,}')

    # Load generated preds and sort by prediction score
    generated = pd.read_csv(generated_path)

    print(f'Size of generated predictions = {len(generated):,}')

    # Compute threshold for max F1 on test preds
    precision, recall, thresholds = precision_recall_curve(
        test[test_tox_column],
        test[test_pred_column]
    )

    f1_scores = 2 * recall * precision / (recall + precision)
    f1_scores[np.isnan(f1_scores)] = 0

    best_f1 = np.max(f1_scores)
    best_f1_index = np.argmax(f1_scores)
    best_precision = precision[best_f1_index]
    best_recall = recall[best_f1_index]
    best_threshold = thresholds[best_f1_index]

    print(f'Best F1 threshold = {best_threshold:.3f}')
    print(f'F1 = {best_f1:.3f}')
    print(f'Precision = {best_precision:.3f}')
    print(f'Recall = {best_recall:.3f}')

    # Get toxic and non-toxic molecules from test set
    toxic = test[test[test_tox_column] == 1]
    nontoxic = test[test[test_tox_column] == 0]

    # Compute percentiles of generated preds within test preds
    toxic_percentiles = [
        percentileofscore(toxic[test_pred_column], pred)
        for pred in generated[generated_pred_column]
    ]
    nontoxic_percentiles = [
        percentileofscore(nontoxic[test_pred_column], pred)
        for pred in generated[generated_pred_column]
    ]

    for i in range(len(generated)):
        print(f'Generated molecule {i}')
        print(f'Prediction = {generated[generated_pred_column].iloc[i]:.3f}')
        print(f'Toxic percentile = {toxic_percentiles[i]:.2f}')
        print(f'Non-toxic percentile = {nontoxic_percentiles[i]:.2f}')
        print()

    # Create save directory
    save_dir.mkdir(parents=True, exist_ok=True)

    # Plot toxicity of test preds and generated preds
    plt.hist(nontoxic[test_pred_column], bins=100, density=True, color='blue', label='non-toxic', alpha=0.5)
    plt.hist(toxic[test_pred_column], bins=100, density=True, color='red', label='toxic', alpha=0.5)
    plt.vlines(generated[generated_pred_column], plt.ylim()[0], 0.1 * plt.ylim()[1], color='black', label='generated')
    plt.legend()
    plt.xlabel('Predicted toxicity')
    plt.ylabel('Density')
    plt.title('Predicted toxicity vs test set predictions')
    plt.savefig(save_dir / 'toxicity.pdf', bbox_inches='tight')

    # Save plot data
    fig_data = {
        'generated': generated[generated_pred_column],
        'test_toxic': toxic[test_pred_column],
        'test_nontoxic': nontoxic[test_pred_column],
    }
    max_len = max(len(values) for values in fig_data.values())
    fig_data = pd.DataFrame({
        key: np.pad(values, (0, max_len - len(values)), constant_values=np.nan)
        for key, values in fig_data.items()
    })
    fig_data.to_csv(save_dir / 'toxicity.csv', index=False)


if __name__ == '__main__':
    from tap import tapify

    tapify(plot_toxicity)
