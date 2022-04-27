"""Generate random molecules combinatorially."""
import json
from itertools import product
from pathlib import Path
from typing import Optional, Union

import numpy as np
import pandas as pd
from rdkit import Chem
from tap import Tap
from tqdm import trange

from real_reactions import convert_to_mol, Reaction, REAL_REACTIONS


class Args(Tap):
    fragment_path: Path  # Path to CSV file containing molecular building blocks.
    smiles_column: str = 'smiles'  # Name of the column containing SMILES.
    reagent_to_fragments_path: Path  # Path to a JSON file containing a dictionary mapping from reagents to fragments.
    num_molecules: int = 10  # Number of molecules to generate.
    max_num_reactions: int = 3  # Maximum number of reactions that can be performed to expand fragments into molecules.
    save_path: Path  # Path to CSV file where generated molecules will be saved.

    def process_args(self) -> None:
        self.save_path.parent.mkdir(parents=True, exist_ok=True)


LogEntry = dict[str, Union[int, list[int], list[str]]]


# TODO: all documentation
def get_reagent_matches_per_mol(reaction: Reaction, mols: list[Chem.Mol]) -> list[list[int]]:
    return [
        [
            reagent_index
            for reagent_index, reagent in enumerate(reaction.reagents)
            if reagent.has_substruct_match(mol)
        ]
        for mol in mols
    ]


def sample_next_fragment(fragments: list[str], reagent_to_fragments: dict[str, set[str]]) -> Optional[str]:
    mols = [convert_to_mol(fragment, add_hs=True) for fragment in fragments]

    # For each reaction these fragments can participate in, get all the unfilled reagents
    unfilled_reagents = set()
    for reaction in REAL_REACTIONS:
        reagent_indices = set(range(reaction.num_reagents))

        # Skip reaction if there's no room to add more reagents
        if len(mols) >= reaction.num_reagents:
            continue

        # For each mol, get a list of indices of reagents it matches
        reagent_matches_per_mol = get_reagent_matches_per_mol(reaction=reaction, mols=mols)

        # Loop through products of reagent indices that the mols match to
        # and for each product, if it matches to all separate reagents,
        # then include the missing reagents in the set of unfilled reagents
        for matched_reagent_indices in product(*reagent_matches_per_mol):
            matched_reagent_indices = set(matched_reagent_indices)

            if len(matched_reagent_indices) == len(mols):
                unfilled_reagents |= {reaction.reagents[index] for index in reagent_indices - matched_reagent_indices}

    if len(unfilled_reagents) == 0:
        return None

    # Get all the fragments that match the other reagents
    available_fragments = sorted(set.union(*[
        reagent_to_fragments[reagent.smarts] for reagent in unfilled_reagents
    ]))
    selected_fragment = np.random.choice(available_fragments)

    return selected_fragment


def find_reactions_for_fragments(fragments: list[str]) -> list[tuple[Reaction, dict[str, int]]]:
    mols = [convert_to_mol(fragment, add_hs=True) for fragment in fragments]
    matching_reactions = []

    for reaction in REAL_REACTIONS:
        if len(mols) != reaction.num_reagents:
            continue

        # For each mol, get a list of indices of reagents it matches
        reagent_matches_per_mol = get_reagent_matches_per_mol(reaction=reaction, mols=mols)

        # Include every assignment of fragments to reagents that fills all the reagents
        for matched_reagent_indices in product(*reagent_matches_per_mol):
            if len(set(matched_reagent_indices)) == reaction.num_reagents:
                fragment_to_reagent_index = dict(zip(fragments, matched_reagent_indices))
                matching_reactions.append((reaction, fragment_to_reagent_index))

    return matching_reactions


def run_random_reaction(fragment: str,
                        fragment_to_index: dict[str, int],
                        reagent_to_fragments: dict[str, set[str]]) -> tuple[str, dict[str, int]]:
    # Create fragments list
    fragments = [fragment]

    # Select second fragment
    second_fragment = sample_next_fragment(fragments=fragments, reagent_to_fragments=reagent_to_fragments)

    if second_fragment is None:
        raise ValueError('Cannot find a second fragment.')

    # Add second fragment
    fragments.append(second_fragment)

    # Possibly select a third fragment
    # TODO: change this probability
    if np.random.rand() < 0.5:
        third_fragment = sample_next_fragment(fragments=fragments, reagent_to_fragments=reagent_to_fragments)

        if third_fragment is not None:
            fragments.append(third_fragment)

    # Find all reactions that match the current fragments
    matching_reactions = find_reactions_for_fragments(fragments=fragments)

    if len(matching_reactions) == 0:
        raise ValueError('Cannot finding a matching reaction.')

    # Sample reaction
    reaction_index = np.random.choice(len(matching_reactions))
    reaction, fragment_to_reagent_index = matching_reactions[reaction_index]

    # Put fragments in the right order for the reaction
    fragments.sort(key=lambda frag: fragment_to_reagent_index[frag])

    # Run reaction
    products = reaction.run_reactants(fragments)

    if len(products) == 0:
        raise ValueError('Reaction failed to produce products.')

    assert all(len(product) == 1 for product in products)

    # Convert product mols to SMILES (and remove Hs)
    products = sorted({Chem.MolToSmiles(Chem.RemoveHs(product[0])) for product in products})

    # Sample a product molecule
    molecule = np.random.choice(products)

    # Create reaction log
    reaction_log = {
        'reaction_id': reaction.id,
        'reagent_ids': [fragment_to_index.get(fragment, -1) for fragment in fragments],
        'reagent_smiles': [fragment for fragment in fragments]
    }

    return molecule, reaction_log


def generate_random_molecule(fragments: list[str],
                             fragment_to_index: dict[str, int],
                             reagent_to_fragments: dict[str, set[str]],
                             max_num_reactions: int) -> tuple[str, list[LogEntry]]:
    # Select first fragment
    fragment = np.random.choice(fragments)

    # Start construction log
    construction_log = []

    # Loop over reactions to expand the fragment
    for reaction_index in range(max_num_reactions):
        # Run a reaction to expand the fragment
        fragment, reaction_log = run_random_reaction(
            fragment=fragment,
            fragment_to_index=fragment_to_index,
            reagent_to_fragments=reagent_to_fragments
        )

        # Add entry to construction log
        construction_log.append(reaction_log)

        # Random chance of stopping early
        # TODO: change this probability
        if np.random.rand() < 0.33:
            break

    return fragment, construction_log


def save_molecules(molecules: list[str],
                   construction_logs: list[list[LogEntry]],
                   save_path: Path) -> None:
    # Convert construction logs from lists to dictionaries
    construction_dicts = []
    max_reaction_num = 0
    reaction_num_to_max_reagent_num = {}

    for construction_log in construction_logs:
        construction_dict = {}
        max_reaction_num = max(max_reaction_num, len(construction_log))

        for reaction_index, reaction_log in enumerate(construction_log):
            reaction_num = reaction_index + 1
            construction_dict[f'reaction_{reaction_num}_id'] = reaction_log['reaction_id']

            reaction_num_to_max_reagent_num[reaction_num] = max(
                reaction_num_to_max_reagent_num.get(reaction_num, 0),
                len(reaction_log['reagent_ids'])
            )

            for reagent_index, (reagent_id, reagent_smiles) in enumerate(zip(reaction_log['reagent_ids'],
                                                                             reaction_log['reagent_smiles'])):
                reagent_num = reagent_index + 1
                construction_dict[f'reagent_{reaction_num}_{reagent_num}_id'] = reagent_id
                construction_dict[f'reagent_{reaction_num}_{reagent_num}_smiles'] = reagent_smiles

        construction_dicts.append(construction_dict)

    # Specify column order for CSV file
    columns = ['smiles']

    for reaction_num in range(1, max_reaction_num + 1):
        columns.append(f'reaction_{reaction_num}')

        for reagent_num in range(1, reaction_num_to_max_reagent_num[reaction_num] + 1):
            columns.append(f'reagent_{reaction_num}_{reagent_num}_id')
            columns.append(f'reagent_{reaction_num}_{reagent_num}_smiles')

    # Save data
    # TODO: save construction logs
    # TODO: specify column orders
    data = pd.DataFrame(data=[
        {
            'smiles': molecule,
            **construction_dict
        }
        for molecule, construction_dict in zip(molecules, construction_dicts)
    ])
    data.to_csv(save_path, index=False)


def generate_random_molecules(args: Args) -> None:
    """Generate random molecules combinatorially."""
    # Load fragments
    fragments = pd.read_csv(args.fragment_path)[args.smiles_column]

    # Map fragments to indices
    fragment_to_index = {}

    for index, fragment in enumerate(fragments):
        if fragment not in fragment_to_index:
            fragment_to_index[fragment] = index

    # Load mapping from reagents to fragments
    with open(args.reagent_to_fragments_path) as f:
        reagent_to_fragments = {
            reagent: set(fragments)
            for reagent, fragments in json.load(f).items()
        }

    # Set random seed
    np.random.seed(0)

    # Get usable fragments
    usable_fragments = sorted(set.union(*reagent_to_fragments.values()))

    # Generate random molecules
    molecules, construction_logs = zip(*[
        generate_random_molecule(
            fragments=usable_fragments,
            fragment_to_index=fragment_to_index,
            reagent_to_fragments=reagent_to_fragments,
            max_num_reactions=args.max_num_reactions
        )
        for _ in trange(args.num_molecules)
    ])

    # Save generated molecules
    save_molecules(
        molecules=molecules,
        construction_logs=construction_logs,
        save_path=args.save_path
    )


if __name__ == '__main__':
    generate_random_molecules(Args().parse_args())
