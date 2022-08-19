"""Counts the number of feasible molecules using reaction SMARTS from the top REAL reactions."""
import json
from collections import defaultdict
from pathlib import Path
from typing import Any, Optional

import numpy as np
import pandas as pd
from rdkit import Chem
from tap import Tap
from tqdm import tqdm

from map_reagents_to_fragments import map_reagents_to_fragments
from real_reactions import REAL_REACTIONS


class Args(Tap):
    fragment_path: Path  # Path to CSV file containing building blocks.
    reagent_to_fragments_path: Optional[Path] = None  # Path to a JSON file containing a pre-computed dictionary mapping from reagents to fragments.
    num_reactions: int = 1  # Number of sequential reactions to perform.
    num_products_per_reaction: int = 100  # Number of product molecules to sample per reaction.
    smiles_column: str = 'smiles'  # Name of the column in molecule_path containing SMILES.


RNG = np.random.default_rng(seed=0)


def random_choice(array: list[Any]) -> Any:
    return array[RNG.integers(0, len(array))]


def count_feasible_molecules(args: Args) -> None:
    """Counts the number of feasible molecules using reaction SMARTS from the top REAL reactions."""
    # Load data
    fragment_data = pd.read_csv(args.fragment_path)

    # Map reagents to fragments
    if args.reagent_to_fragments_path is not None:
        with open(args.reagent_to_fragments_path) as f:
            reagent_to_fragments = json.load(f)
    else:
        reagent_to_fragments = map_reagents_to_fragments(fragments=fragment_data[args.smiles_column])

    # Map reagents to current set of molecules (which starts as just the fragments)
    reagent_to_molecules = {
        reagent: [
            {
                'smiles': fragment,
                'weight': 1
            }
            for fragment in fragments
        ]
        for reagent, fragments in reagent_to_fragments.items()
    }

    # Count feasible molecules for each reaction for each reaction step
    for reaction_step in range(args.num_reactions):
        print(f'Reaction step {reaction_step + 1}')
        all_products = []

        # Count feasible molecules for each reaction and sample products
        for reaction in tqdm(REAL_REACTIONS, desc='Looping reactions'):
            reagent_indices = set(range(reaction.num_reagents))
            num_possible_products = 0
            reaction_products = []

            # Avoid double counting when running first reaction
            molecule_index_range = range(1 if reaction_step == 0 else reaction.num_reagents)

            # Loop over reagents, treating each one as the source of the current molecule in turn
            for molecule_reagent_index in molecule_index_range:
                # Get molecule SMILES and weights
                reagent_str = str(reaction.reagents[molecule_reagent_index])
                molecule_smiles = [molecule['smiles'] for molecule in reagent_to_molecules[reagent_str]]
                molecule_weights = np.array([molecule['weight'] for molecule in reagent_to_molecules[reagent_str]])
                normed_molecule_weights = molecule_weights / np.sum(molecule_weights)

                # Count and sample SMILES for the molecule
                reagent_index_to_data = {
                    molecule_reagent_index: {
                        'count': sum(molecule_weights),
                        'sampled_smiles': RNG.choice(
                            molecule_smiles, size=args.num_products_per_reaction, p=normed_molecule_weights
                        )
                    }
                }

                # Count and sample SMILES for the remaining fragments
                for fragment_reagent_index in reagent_indices - {molecule_reagent_index}:
                    # Get fragment SMILES
                    reagent_str = str(reaction.reagents[fragment_reagent_index])
                    fragment_smiles = reagent_to_fragments[reagent_str]

                    # Count and sample SMILES for the fragment
                    reagent_index_to_data[fragment_reagent_index] = {
                        'count': len(fragment_smiles),
                        'sampled_smiles': RNG.choice(fragment_smiles, size=args.num_products_per_reaction)
                    }

                # Count number of feasible molecules
                # TODO: clean up the diff
                num_possible_products += reaction.count_feasible_products(*[
                    reagent_index_to_data[i]['count'] for i in range(reaction.num_reagents)
                ], diff=reaction_step > 0 and not (reaction.id == 1 and molecule_reagent_index == 0))

                # Create products from sampled reagents
                for i in range(args.num_products_per_reaction):
                    products = reaction.run_reactants([
                        reagent_index_to_data[index]['sampled_smiles'][i] for index in reagent_indices
                    ])

                    if len(products) == 0:
                        raise ValueError('Reaction failed to produce products.')

                    assert all(len(product) == 1 for product in products)

                    # Convert mols to SMILES, remove Hs, and deduplicate (preserving order for random reproducibility)
                    products = list(dict.fromkeys(Chem.MolToSmiles(Chem.RemoveHs(product[0])) for product in products))

                    # Sample a product molecule
                    product = random_choice(products)
                    reaction_products.append(product)

            print(f'Reaction {reaction.id} can produce {num_possible_products:.2e} molecules')

            # Add weights to products
            weight = num_possible_products / len(reaction_products)

            all_products += [
                {
                    'smiles': product,
                    'weight': weight
                }
                for product in reaction_products
            ]

        if reaction_step < args.num_reactions - 1:
            # Update mapping from reagent to molecules
            product_smiles = [product['smiles'] for product in all_products]
            product_weight = [product['weight'] for product in all_products]

            product_to_weight = defaultdict(int)
            for smiles, weight in zip(product_smiles, product_weight):
                product_to_weight[smiles] += weight

            reagent_to_molecules = map_reagents_to_fragments(fragments=product_smiles)

            reagent_to_molecules = {
                reagent: [
                    {
                        'smiles': molecule,
                        'weight': product_to_weight[molecule]
                    }
                    for molecule in molecules
                ]
                for reagent, molecules in reagent_to_molecules.items()
            }


if __name__ == '__main__':
    count_feasible_molecules(Args().parse_args())
