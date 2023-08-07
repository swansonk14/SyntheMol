import pandas as pd
from pathlib import Path
from tap import tapify 


def calc_average(preds_path: Path):
    test_preds = pd.read_csv(preds_path)
    test_preds['smiles'] = test_preds['smiles'].str.strip('[]')
    test_preds[['smiles', 'solvent_smiles']] = test_preds['smiles'].str.split(',', expand=True)
    test_preds['smiles'] = test_preds['smiles'].str.strip(' \'')
    test_preds['solvent_smiles'] = test_preds['solvent_smiles'].str.strip(' \'')
    group_df = test_preds[['smiles', 'PLQY']].groupby(by='smiles').mean()
    group_df.to_csv('./average_scores.csv')

if __name__ == '__main__':
    tapify(calc_average)
    