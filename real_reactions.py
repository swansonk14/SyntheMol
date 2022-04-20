"""SMARTS representations of the REAL reactions."""
# Dictionary mapping our reaction ID to SMARTS for the reagents and products
REAL_REACTIONS = {
    1: {
        'reagents': [
            'CC(C)(C)OC(=O)[N:1]([*:2])C[NH1:3][*:4]'  # TODO: figure out carbon chain length
            '[OH1][C:5]([*:6])=[O:7]',
            '[OH1][C:8]([*:9])=[O:10]'
        ],
        'products': 'C[N:3]([*:4])[C:5](=[O:7])[*:6].C[N:1]([*:2])[C:8](=[O:10])[*:9]'  # TODO: figure out carbon chain
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
    5: {2230}
}
