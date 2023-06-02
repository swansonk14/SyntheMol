"""Maps generated molecules to REAL IDs in the format expected by Enamine."""
from pathlib import Path

import pandas as pd
from rdkit import Chem
from rdkit.Chem.PandasTools import WriteSDF


def map_generated_molecules_to_real_ids(
        data_path: Path,
        save_dir: Path
) -> None:
    """Maps generated molecules to REAL IDs in the format expected by Enamine.

    Note: Currently only works with one-reaction molecules.

    :param data_path: Path to CSV file containing generated molecules with reaction and building block IDs.
    :param save_dir: Path to directory where CSV and SDF files with REAL IDs will be saved.
    """
    # Load data
    data = pd.read_csv(data_path)
    print(f'Data size = {len(data):,}')

    # Get only one-reaction molecules
    data = data[data['num_reactions'] == 1]
    print(f'Number of one-reaction molecules = {len(data):,}')

    # Get building block columns
    building_block_columns = sorted(
        column for column in data.columns if column.startswith('building_block_1_') and column.endswith('_id')
    )

    # Compute REAL IDs without prefixes
    real_ids = [
        f'{row["reaction_1_id"]}____{"____".join(str(int(building_block_id)) for building_block_id in row[building_block_columns].dropna())}'
        for _, row in data.iterrows()
    ]

    # Compute mols
    mols = [Chem.MolFromSmiles(smiles) for smiles in data['smiles']]

    # Create save directory
    save_dir.mkdir(parents=True, exist_ok=True)

    # Loop through prefixes
    for prefix in ['m', 's']:
        # Create new DataFrame with molecules and REAL IDs
        real_data = pd.DataFrame(data={
            'real_id': [f'{prefix}_{real_id}' for real_id in real_ids],
            'smiles': data['smiles'],
            'mol': mols
        })

        # Save data as SDF
        with open(save_dir / f'type_{prefix}.sdf', 'w') as f:
            WriteSDF(real_data, f, molColName='mol', idName='real_id')

        # Save data as CSV
        del real_data['mol']
        real_data.to_csv(save_dir / f'type_{prefix}.csv', index=False)


if __name__ == '__main__':
    from tap import tapify

    tapify(map_generated_molecules_to_real_ids)
