"""Maps generated molecules to REAL IDs in the format expected by Enamine."""
from pathlib import Path

import pandas as pd
from rdkit import Chem
from rdkit.Chem.PandasTools import WriteSDF
from tap import Tap

from real_reactions import REACTION_ID_TO_REAL_IDS


class Args(Tap):
    data_path: Path  # Path to CSV file containing generated molecules with reaction and reagent IDs.
    smiles_save_path: Path  # Path to CSV file where molecule SMILES and REAL IDs will be saved.
    sdf_save_path: Path  # Path to SDF file where molecules and REAL IDs will be saved.

    def process_args(self) -> None:
        self.smiles_save_path.parent.mkdir(parents=True, exist_ok=True)
        self.sdf_save_path.parent.mkdir(parents=True, exist_ok=True)


def map_generated_molecules_to_real_ids(args: Args) -> None:
    """Maps generated molecules to REAL IDs in the format expected by Enamine."""
    # Load data
    data = pd.read_csv(args.data_path)
    print(f'Data size = {len(data):,}')

    # Get only one-reaction molecules
    data = data[data['num_reactions'] == 1]
    print(f'Number of one-reaction molecules = {len(data):,}')

    # Get reagent columns
    reagent_columns = sorted(
        column for column in data.columns if column.startswith('reagent_1_') and column.endswith('_id')
    )

    # Create new DataFrame with SMILES and REAL IDs
    smiles, mols, real_ids = [], [], []
    for _, row in data.iterrows():
        reagent_ids = '____'.join(str(int(reagent_id)) for reagent_id in row[reagent_columns].dropna())
        mol = Chem.MolFromSmiles(row['smiles'])

        for reaction_id in sorted(REACTION_ID_TO_REAL_IDS[row['reaction_1_id']]):
            smiles.append(row['smiles'])
            mols.append(mol)
            real_ids.append(f'm_{reaction_id}____{reagent_ids}')

    real_data = pd.DataFrame(data={
        'real_id': real_ids,
        'smiles': smiles,
        'mol': mols
    })
    print(f'Number of unique REAL IDs = {len(real_data):,}')

    # Save data as SDF
    with open(args.sdf_save_path, 'w') as f:
        WriteSDF(real_data, f, molColName='mol', idName='real_id')

    # Save data as CSV
    del real_data['mol']
    real_data.to_csv(args.smiles_save_path, index=False)


if __name__ == '__main__':
    map_generated_molecules_to_real_ids(Args().parse_args())
