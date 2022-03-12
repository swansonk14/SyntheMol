"""Converts molecules in SDF format to a CSV with SMILES."""
from pathlib import Path

import pandas as pd
from rdkit import Chem
from tap import Tap
from tqdm import tqdm


class Args(Tap):
    data_path: Path  # Path to an SDF file.
    save_path: Path  # Path to a CSV file where SMILES strings will be saved.

    def process_args(self) -> None:
        self.save_path.parent.mkdir(parents=True, exist_ok=True)


def sdf_to_smiles(args: Args) -> None:
    """Converts molecules in SDF format to a CSV with SMILES."""
    # Load SDF file
    num_skipped = 0
    mols = []

    with Chem.SDMolSupplier(str(args.data_path)) as suppl:
        for mol in tqdm(suppl, desc='Loading SDF'):
            if mol is None:
                num_skipped += 1
            else:
                mols.append(mol)

    print(f'Number of molecules = {len(mols):,}')
    print(f'Number skipped = {num_skipped:,}')

    # Convert Mols to (unique) SMILES
    smiles = sorted({Chem.MolToSmiles(mol) for mol in tqdm(mols, desc='Mol to SMILES')})

    print(f'Number of unique SMILES = {len(smiles):,}')

    # Save data as CSV
    data = pd.DataFrame(data={'smiles': smiles})
    data.to_csv(args.save_path, index=False)


if __name__ == '__main__':
    sdf_to_smiles(Args().parse_args())
