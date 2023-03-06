"""Export reaction SMARTS and reaction IDs in a CSV file."""
from pathlib import Path

import pandas as pd
from tap import tapify

from SyntheMol.reactions import REAL_REACTIONS


def export_reaction_smarts(
        save_path: Path
) -> None:
    """Export reaction SMARTS and reaction IDs in a CSV file.

    :param save_path: Path to a CSV file where reaction SMARTS and reaction IDs will be saved.
    """
    # Gather SMARTS and reaction IDs from REAL reactions
    smarts, reaction_ids = [], []
    for reaction in REAL_REACTIONS:
        smarts.append(
            f'{".".join(f"({reactant.smarts_with_atom_mapping})" for reactant in reaction.reactants)}'
            f'>>({reaction.product.smarts_with_atom_mapping})'
        )
        reaction_ids.append(reaction.id)

    # Save reaction SMARTS and reaction IDs in a CSV file
    save_path.parent.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame({'smarts': smarts, 'reaction_id': reaction_ids})
    df.to_csv(save_path, index=False)


if __name__ == '__main__':
    tapify(export_reaction_smarts)
