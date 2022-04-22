"""SMARTS representations of the REAL reactions."""
from rdkit import Chem


class SidechainChecker:
    matchers = {
        'alkyl': {
            'function': lambda atom: atom.GetAtomicNum() == 6 and not atom.GetIsAromatic(),
            'type': 'all'
        },
        'aryl': {
            'function': lambda atom: atom.GetIsAromatic(),
            'type': 'any'
        },
        'cycle': {
            'function': lambda atom: atom.IsInRing(),
            'type': 'any'
        }
    }

    def __init__(self, atoms_to_check: dict[int, str]) -> None:
        self.atoms_to_check = atoms_to_check  # Maps from index in query SMARTS of a side chain to match type

    def __call__(self, mol: Chem.Mol, matched_atoms: list[int]) -> bool:
        visited_atoms = set(matched_atoms)

        # Loop over the atoms that start the side chains to check
        for query_index, match_type in self.atoms_to_check.items():
            # Get the index of the atom in the molecule
            atom_index = matched_atoms[query_index]
            atom = mol.GetAtomWithIdx(atom_index)

            # Do a breadth-first search from that atom, checking all of its neighbors that aren't in the query
            stack = [atom]
            matcher_function = self.matchers[match_type]['function']
            matcher_type = self.matchers[match_type]['type']
            matched = matcher_type == 'all'

            while stack:
                atom = stack.pop(0)

                if matcher_type == 'all' and not matcher_function(atom):
                    return False
                elif matcher_function(atom):
                    matched = True
                    break

                visited_atoms.add(atom.GetIdx())

                for neighbor_atom in atom.GetNeighbors():
                    if neighbor_atom.GetIdx() not in visited_atoms:
                        stack.append(neighbor_atom)

            if not matched:
                return False

        return True


# Dictionary mapping our reaction ID to SMARTS for the reagents and products
REAL_REACTIONS = {
    1: {
        'reagents': [
            {
                'smarts': 'CC(C)(C)OC(=O)[N:1]([*:2])C.C[NH1:3][*:4]',  # TODO: fix C matching
                'checker': None  # TODO: figure out carbon chain length (note: chain could be non-aromatic ring)
            },
            {
                'smarts': '[OH1][C:5]([*:6])=[O:7]',
                'checker': None
            },
            {
                'smarts': '[OH1][C:8]([*:9])=[O:10]',
                'checker': None
            }
        ],
        'product': {  # TODO: is a checker needed for the product and/or reaction smarts?
            # TODO: is atom mapping needed
            'smarts': 'C[N:3]([*:4])[C:5](=[O:7])[*:6].C[N:1]([*:2])[C:8](=[O:10])[*:9]',
            'checker': None  # TODO: figure out carbon chain
        }
    },
    2: {
        'reagents': [
            '[*:1][NH1:2][*:3]',
            '[OH1][C:4]([*:5])=[O:6]'
        ],
        'product': '[*:5][C:4](=[O:6])[N:2]([*:1])[*:3]'
    },
    3: {  # TODO: Ar?
        'reagents': [

        ],
        'product': ''
    },
    4: {
        'reagents': [
            '[*:1][NH2:2]',
            '[*:3][NH2:4]',
        ],
        'product': 'O=C([NH1:2][*:1])[NH1:4][*:3]'
    },
    5: {  # TODO: Alk, Ar, cycle?
        'reagents': [

        ],
        'product': ''
    },
    6: {
        'reagents': [
            {
                'smarts': '[NH2:1][*:2]',
                'checker': SidechainChecker(atoms_to_check={1: 'aryl'})
            },
            {
                'smarts': '[NH2:3][*:4]',
                'checker': SidechainChecker(atoms_to_check={1: 'alkyl'})
            }
        ],
        'product': {
            'smarts': 'O=C([NH1:1][*:2])C(=O)[NH1:3][*:4]',
            'checker': SidechainChecker(atoms_to_check={3: 'aryl', 7: 'alkyl'})
        }
    }
}

# Create reaction SMARTS from reagent and product SMARTS
for reaction in REAL_REACTIONS.values():
    # TODO: do we need parentheses for the reagents to force them to be separate molecules?
    # reaction = reagent_1.reagent_2...reagent_n>>product
    reaction['reaction'] = f'{".".join(reaction["reagent"])}>>{reaction["product"]}'

# Dictionary mapping our reaction ID to a set of corresponding REAL reaction IDs
REACTION_ID_TO_REAL_ID = {
    1: {275592},
    2: {11, 22, 527, 240690},
    3: {2430},
    4: {2708},
    5: {2230},
    6: {2718}
}
