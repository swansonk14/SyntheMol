"""Maps generated molecules to REAL IDs in the format expected by Enamine."""
from pathlib import Path

import pandas as pd
from rdkit import Chem
from rdkit.Chem.PandasTools import WriteSDF
from tap import Tap


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
    real_data = pd.DataFrame(data=[
        {
            'real_id': f'm_{row["reaction_1_id"]}____'
                       f'{"____".join(str(int(reagent_id)) for reagent_id in row[reagent_columns].dropna())}',
            'smiles': row['smiles'],
            'mol': Chem.MolFromSmiles(row['smiles'])
        } for _, row in data.iterrows()
    ])

    # Save data as SDF
    with open(args.sdf_save_path, 'w') as f:
        WriteSDF(real_data, f, molColName='mol', idName='real_id')

    # Save data as CSV
    del real_data['mol']
    real_data.to_csv(args.smiles_save_path, index=False)


if __name__ == '__main__':
    map_generated_molecules_to_real_ids(Args().parse_args())
