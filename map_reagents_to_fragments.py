"""Maps REAL reagents to fragments that match those reagents."""
import json
from pathlib import Path

import pandas as pd
from tap import Tap
from tqdm import tqdm

from reactions import convert_to_mol
from real_reactions import REAL_REACTIONS
from synnet_reactions import SYNNET_REACTIONS


class Args(Tap):
    fragment_path: Path  # Path to CSV file containing REAL fragments.
    smiles_column: str = 'smiles'  # Name of the column containing SMILES.
    save_path: Path  # Path to JSON file where a dictionary mapping reagents to fragments will be saved.
    synnet_rxn: bool = False  # Whether to include SynNet reactions in addition to REAL reactions.

    def process_args(self) -> None:
        self.save_path.parent.mkdir(parents=True, exist_ok=True)


def map_reagents_to_fragments(fragments: list[str], synnet_rxn: bool = False) -> dict[str, list[str]]:
    """Maps REAL reagents to fragments that match those reagents.

    :param fragments: A list of molecular fragments (SMILES).
    :param synnet_rxn: Whether to include SynNet reactions in addition to REAL reactions.
    :return: A dictionary mapping REAL reagent SMARTS to a sorted list of fragment SMILES that match the reagent.
    """
    if synnet_rxn:
        reactions = REAL_REACTIONS + SYNNET_REACTIONS
    else:
        reactions = REAL_REACTIONS

    # Get fragment SMILES
    fragments = sorted(set(fragments))

    # Convert fragment SMILES to mols
    fragment_mols = [
        convert_to_mol(fragment, add_hs=True)
        for fragment in tqdm(fragments, desc='SMILES to mol')
    ]

    # Map reagents to fragments
    reagent_to_fragments = {}

    for reaction in tqdm(reactions, desc='Looping reactions'):
        for reagent in tqdm(reaction.reagents, desc='Looping reagents', leave=False):
            reagent_str = str(reagent)

            if reagent_str in reagent_to_fragments:
                continue

            reagent_to_fragments[reagent_str] = [
                fragment
                for fragment, fragment_mol in tqdm(zip(fragments, fragment_mols),
                                                   total=len(fragments), desc='Matching fragments', leave=False)
                if reagent.has_substruct_match(fragment_mol)
            ]

    return reagent_to_fragments


def run_map_reagents_to_fragments(args: Args) -> None:
    """Maps REAL reagents to fragments that match those reagents."""
    # Load data
    fragment_data = pd.read_csv(args.fragment_path)

    # Map reagents to fragments
    reagent_to_fragments = map_reagents_to_fragments(
        fragments=fragment_data[args.smiles_column],
        synnet_rxn=args.synnet_rxn
    )

    # Print statistics
    for reagent, fragments in reagent_to_fragments.items():
        print(f'Reagent {reagent} matches {len(fragments):,} fragments')

    # Save data
    with open(args.save_path, 'w') as f:
        json.dump(reagent_to_fragments, f, indent=4, sort_keys=True)


if __name__ == '__main__':
    run_map_reagents_to_fragments(Args().parse_args())
