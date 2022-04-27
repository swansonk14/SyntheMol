"""Generate random molecules combinatorially."""
import json
from itertools import product
from pathlib import Path
from typing import Any, Optional

import numpy as np
import pandas as pd
from rdkit import Chem
from tap import Tap
from tqdm import trange

from real_reactions import Reaction, REAL_REACTIONS


class Args(Tap):
    fragment_path: Path  # Path to CSV file containing molecular building blocks.
    smiles_column: str = 'smiles'  # Name of the column containing SMILES.
    reagent_to_fragments_path: Path  # Path to a JSON file containing a dictionary mapping from reagents to fragments.
    num_molecules: int = 10  # Number of molecules to generate.
    max_num_reactions: int = 3  # Maximum number of reactions that can be performed to expand fragments into molecules.
    save_path: Path  # Path to CSV file where generated molecules will be saved.

    def process_args(self) -> None:
        self.save_path.parent.mkdir(parents=True, exist_ok=True)


def get_reagent_matches_per_mol(reaction: Reaction, mols: list[Chem.Mol]) -> list[list[int]]:
    return [
        [
            reagent_index
            for reagent_index, reagent in enumerate(reaction.reagents)
            if reagent.has_substruct_match(mol)
        ]
        for mol in mols
    ]


def sample_next_fragment(fragments: list[str], reagent_to_fragments: dict[str, list[str]]) -> Optional[str]:
    mols = [Chem.AddHs(Chem.MolFromSmiles(fragment)) for fragment in fragments]

    # Get all the other reagents that match a reaction this fragment can participate in
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
        for indices in product(*reagent_matches_per_mol):
            index_set = set(indices)

            if len(index_set) == len(mols):
                unfilled_reagents |= {reaction.reagents[i] for i in reagent_indices - index_set}

    if len(unfilled_reagents) == 0:
        return None

    # Get all the fragments that match the other reagents
    available_fragments = sorted(set.union(*[
        set(reagent_to_fragments[reagent.smarts]) for reagent in unfilled_reagents
    ]))
    selected_fragment = np.random.choice(available_fragments)

    return selected_fragment


def find_reactions_for_fragments(fragments: list[str]) -> list[tuple[Reaction, tuple[int, ...]]]:
    mols = [Chem.AddHs(Chem.MolFromSmiles(fragment)) for fragment in fragments]
    matching_reactions = []

    for reaction in REAL_REACTIONS:
        if len(mols) != reaction.num_reagents:
            continue

        # For each mol, get a list of indices of reagents it matches
        reagent_matches_per_mol = get_reagent_matches_per_mol(reaction=reaction, mols=mols)

        # If any assignment of fragments to reagents fills all the reagents, then include the reaction
        for reagent_indices in product(*reagent_matches_per_mol):
            if len(set(reagent_indices)) == reaction.num_reagents:
                matching_reactions.append((reaction, reagent_indices))

    return matching_reactions


def run_random_reaction(fragment: str, reagent_to_fragments: dict[str, list[str]]) -> str:
    """Given a fragment, do the following:
    1) Find all the reactants that match this fragment.
    2) Find all the reactions that contain those reactants.
    3) Find all the other reactants for each of those reactions.
    4) Find all the fragments that match those other reactants.
    5) Randomly select a fragment from that list.
    6) Repeat if more than two reactants.
    7) Select reaction among all the reactions that match those reactants.
    8) Run that reactions to generate the molecule.
    """
    # Create fragments list
    fragments = [fragment]

    # Select second fragment
    second_fragment = sample_next_fragment(fragments=fragments, reagent_to_fragments=reagent_to_fragments)

    if second_fragment is None:
        raise ValueError('Cannot find a second fragment.')

    fragments.append(second_fragment)

    # Possibly select a third fragment
    # TODO: change this probability
    if np.random.rand() < 0.5:
        third_fragment = sample_next_fragment(fragments=fragments, reagent_to_fragments=reagent_to_fragments)

        if third_fragment is not None:
            fragments.append(third_fragment)

    # Find all reactions that match the current fragments
    reactions_with_reagent_indices = find_reactions_for_fragments(fragments=fragments)

    if len(reactions_with_reagent_indices) == 0:
        raise ValueError('Cannot finding a matching reaction.')

    # Sample reaction
    reaction_indices = np.arange(len(reactions_with_reagent_indices))
    reaction_index = np.random.choice(reaction_indices)
    reaction, reagent_indices = reactions_with_reagent_indices[reaction_index]

    # Put fragments in the right order for the reaction
    fragment_to_index = dict(zip(fragments, reagent_indices))
    fragments.sort(key=lambda frag: fragment_to_index[frag])

    # Run reaction
    products = reaction.run_reactants(fragments)

    if len(products) == 0:
        raise ValueError('Reaction failed to produce products.')

    assert all(len(product) == 1 for product in products)

    # Convert product mols to SMILES (and remove Hs)
    products = sorted({Chem.MolToSmiles(Chem.RemoveHs(product[0])) for product in products})

    # Sample a product molecule
    molecule = np.random.choice(products)

    return molecule


def generate_random_molecule(fragments: list[str],
                             reagent_to_fragments: dict[str, list[str]],
                             max_num_reactions: int) -> str:
    fragment = np.random.choice(fragments)

    for _ in range(max_num_reactions):
        # Run a reaction to expand the fragment
        fragment = run_random_reaction(fragment=fragment, reagent_to_fragments=reagent_to_fragments)

        # Random chance of stopping early
        # TODO: change this probability
        if np.random.rand() < 0.33:
            break

    return fragment


def generate_random_molecules(args: Args) -> None:
    """Generate random molecules combinatorially."""
    # Load fragments
    fragments = sorted(pd.read_csv(args.fragment_path)[args.smiles_column])

    # Load mapping from reagents to fragments
    with open(args.reagent_to_fragments_path) as f:
        reagent_to_fragments = json.load(f)

    # Set random seed
    np.random.seed(0)

    # Generate random molecules
    # TODO: keep track of building blocks and reactions that are used
    molecules = [
        generate_random_molecule(
            fragments=fragments,
            reagent_to_fragments=reagent_to_fragments,
            max_num_reactions=args.max_num_reactions
        )
        for _ in trange(args.num_molecules)
    ]

    # Save data
    data = pd.DataFrame(data={'smiles': molecules})
    data.to_csv(args.save_path, index=False)


if __name__ == '__main__':
    generate_random_molecules(Args().parse_args())
