"""SMARTS representations of the REAL reactions."""
# Dictionary mapping our reaction ID to SMARTS for the reagents and products
REAL_REACTIONS = {
    1: {
        'reagents': [
            ''
            '[OH1][C:4]([*:5])=[O:6]',
            '[OH1][C:7]([*:8])=[O:9]'
        ],
        'products': ''
    },
    2: {
        'reagents': [
            '[*:1][NH1:2][*:3]',
            '[OH1][C:4]([*:5])=[O:6]'
        ],
        'product': '[*:5][C:4](=[O:6])[N:2]([*:1])[*:3]'
    }
}

# Create reaction SMARTS from reagent and product SMARTS
for reaction in REAL_REACTIONS.values():
    # TODO: do we need parentheses for the reagents to force them to be separate molecules?
    # reaction = reagent_1.reagent_2...reagent_n>>product
    reaction['reaction'] = f'{".".join(reaction["reagent"])}>>{reaction["product"]}'

# Dictionary mapping our reaction ID to a set of corresponding REAL reaction IDs
# TODO: check if our reactions map to even more REAL reactions with same reactants/products
REACTION_ID_TO_REAL_ID = {
    1: {275592},
    2: {11, 22, 527, 240690}
}
