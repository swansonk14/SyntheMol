"""Maps generated molecules to REAL IDs in the format expected by Enamine."""
from pathlib import Path

import pandas as pd
from rdkit import Chem
from rdkit.Chem.PandasTools import WriteSDF
from tap import Tap


class Args(Tap):
    data_path: Path  # Path to CSV file containing generated molecules with reaction and reagent IDs.
    save_dir: Path  # Path to directory where CSV and SDF files with REAL IDs will be saved.

    def process_args(self) -> None:
        self.save_dir.mkdir(parents=True, exist_ok=True)


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

    # Compute REAL IDs without prefixes
    real_ids = [
        f'{row["reaction_1_id"]}____{"____".join(str(int(reagent_id)) for reagent_id in row[reagent_columns].dropna())}'
        for _, row in data.iterrows()
    ]

    # Compute mols
    mols = [Chem.MolFromSmiles(smiles) for smiles in data['smiles']]

    # Loop through prefixes
    for prefix in ['m', 's']:
        # Create new DataFrame with molecules and REAL IDs
        real_data = pd.DataFrame(data={
            'real_id': [f'{prefix}_{real_id}' for real_id in real_ids],
            'smiles': data['smiles'],
            'mol': mols
        })

        # Save data as SDF
        with open(args.save_dir / f'type_{prefix}.sdf', 'w') as f:
            WriteSDF(real_data, f, molColName='mol', idName='real_id')

        # Save data as CSV
        del real_data['mol']
        real_data.to_csv(args.save_dir / f'type_{prefix}.csv', index=False)


if __name__ == '__main__':
    map_generated_molecules_to_real_ids(Args().parse_args())
