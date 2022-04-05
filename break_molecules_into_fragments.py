"""Break molecules into non-overlapping fragments with minimal remaining atoms/bonds."""
from itertools import combinations, product
from pathlib import Path

import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.FragmentMatcher import FragmentMatcher
from tap import Tap
from tqdm import tqdm


class Args(Tap):
    fragment_path: Path  # Path to CSV file containing molecular fragments.
    molecule_path: Path  # Path to CSV file containing full molecules.
    smiles_column: str = 'smiles'  # Name of the column in molecule_path containing SMILES.
    save_path: Path  # Path to CSV file where fragment indicators for each molecule will be saved.

    def process_args(self) -> None:
        self.save_path.parent.mkdir(parents=True, exist_ok=True)


def break_molecules_into_fragments(args: Args) -> None:
    """Break molecules into non-overlapping fragments with minimal remaining atoms/bonds."""
    # Load fragments
    # TODO: canonicalize
    fragments = pd.read_csv(args.fragment_path)
    fragments = sorted(set(fragments[args.smiles_column]))

    # Map from fragment SMILES to fragment matcher
    frag_matchers = {}
    for fragment in tqdm(fragments):
        frag_matcher = FragmentMatcher()
        frag_matcher.Init(fragment)
        frag_matchers[fragment] = frag_matcher

    # Load molecules
    molecules = pd.read_csv(args.molecule_path, sep='\t')

    # Convert SMILES to RDKit molecules
    molecules['mol'] = [Chem.MolFromSmiles(smiles) for smiles in tqdm(molecules[args.smiles_column])]

    # Filter out bad molecules
    molecules = molecules[~molecules['mol'].isna()]

    # Find non-overlapping fragments in molecules
    for mol in tqdm(molecules['mol']):
        # Get fragment matches in the molecule
        frags_to_match = {
            f'{fragment}_{i}': set(match)
            for fragment in fragments
            for i, match in enumerate(frag_matchers[fragment].GetMatches(mol))
        }
        frags_in_molecule = sorted(frags_to_match)
        print(f'Number of fragments = {len(frags_in_molecule)}')

        # Find non-overlapping sets of two fragments in a molecule
        non_overlapping_pairs = {}
        for frag_1, frag_2 in combinations(frags_in_molecule, r=2):
            atoms_1, atoms_2 = frags_to_match[frag_1], frags_to_match[frag_2]
            union_atoms = atoms_1 | atoms_2

            # If non-overlapping fragments, combine
            if len(union_atoms) == len(atoms_1) + len(atoms_2):
                frag_pair = tuple(sorted([frag_1, frag_2]))
                non_overlapping_pairs[frag_pair] = union_atoms

        frag_pairs_in_molecule = sorted(non_overlapping_pairs)
        print(f'Number of non-overlapping pairs = {len(frag_pairs_in_molecule)}')

        # Find non-overlapping sets of three fragments in a molecule
        non_overlapping_triplets = {}
        for frag, frag_pair in product(frags_in_molecule, frag_pairs_in_molecule):
            atoms, pair_atoms = frags_to_match[frag], non_overlapping_pairs[frag_pair]
            union_atoms = atoms | pair_atoms

            # if non-overlapping fragment and fragment pair, combine
            if len(union_atoms) == len(atoms) + len(pair_atoms):
                frag_triplet = tuple(sorted([frag, *frag_pair]))
                non_overlapping_triplets[frag_triplet] = union_atoms

        frag_triplets_in_molecule = sorted(non_overlapping_triplets)
        print(f'Number of non-overlapping triplets = {len(frag_triplets_in_molecule)}')

        frag_triplets_in_molecule = sorted(frag_triplets_in_molecule,
                                           key=lambda frag_triplet: len(non_overlapping_triplets[frag_triplet]),
                                           reverse=True)

        print(frag_triplets_in_molecule[2:4])
        triplet_coverage = [len(non_overlapping_triplets[frag_triplet]) for frag_triplet in frag_triplets_in_molecule]
        print(triplet_coverage[:10])
        exit()

        for frag_triplet in frag_triplets_in_molecule[:4]:
            print([frags_to_match[frag] for frag in frag_triplet])
            img = Draw.MolsToGridImage(
                [mol] * 4,
                molsPerRow=4,
                subImgSize=(600, 600),
                highlightAtomLists=[frags_to_match[frag] for frag in frag_triplet] + [set(range(mol.GetNumAtoms())) - non_overlapping_triplets[frag_triplet]]
            )
            img.show()
            img.close()
        exit()

        triplet_coverage = sorted((len(atoms) for atoms in non_overlapping_triplets.values()), reverse=True)
        print(mol.GetNumAtoms())
        print(triplet_coverage[:10])
        max_cov = triplet_coverage[0]
        print(sum(cov == max_cov for cov in triplet_coverage))

    raise NotImplementedError

    # Delete RDKit molecules
    del molecules['mol']

    # Save data
    molecules.to_csv(args.save_path, index=False)


if __name__ == '__main__':
    break_molecules_into_fragments(Args().parse_args())
