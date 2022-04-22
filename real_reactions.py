"""SMARTS representations of the REAL reactions."""
# TODO: other reaction IDs that match
# TODO: requires Chem.AddHs()
REAL_REACTIONS = [
    {
        'reagents': [
            'CC(C)(C)OC(=O)[N:1]([*:2])C.C[NH1:3][*:4]',  # TODO: fix C matching
            '[OH1][C:5]([*:6])=[O:7]',
            '[OH1][C:8]([*:9])=[O:10]'
        ],
        'product': 'C[N:3]([*:4])[C:5](=[O:7])[*:6].C[N:1]([*:2])[C:8](=[O:10])[*:9]',  # TODO: fix C matching
        'real_ids': {275592}
    },
    {
        'reagents': [
            '[*:1][NH1:2][*:3]',
            '[OH1][C:4]([*:5])=[O:6]'
        ],
        'product': '[*:5][C:4](=[O:6])[N:2]([*:1])[*:3]',
        'real_ids': {11, 22, 527, 188690, 240690, 269946, 270006, 270112, 272270, 272430, 273450, 273652, 280130}
    },
    {
        'reagents': [
            '[*:1][NH1:2][*:3]',
            '[*:4][NH1:5][*:6]'
        ],
        'product': 'O=C([N:2]([*:1])[*:3])[N:5]([*:4])[*:6]',
        'real_ids': {512, 2430, 2554, 2708, 272164, 273390, 273392, 273452, 273454, 273458, 273464, 273574}
    },
    {
        'reagents': [
            '[*:1][NH1:2][*:3]',
            '[F,Cl,Br,I:4][*:5]'
        ],
        'product': '[*:1][N:2]([*:2])[*:3]',
        'real_ids': {7, 27, 34, 38, 44, 61, 2230, 269956, 269982, 270122, 270166, 270344, 272692, 272710, 273654}
    },
    {
        'reagents': [
            '[*:1][O,S:2]',
            '[F,Cl,Br,I:3][*:4]'
        ],
        'product': '[*:1][O,S:2][*:4]',
        'real_ids': {7, 34, 272692, 272710, 273654}
    },
    {
        'reagents': [
            '[*:1][NH1:2][*:3]',
            '[*:4][NH1:5][*:6]',
        ],
        'product': 'O=C(C(=O)[N:2]([*:1])[*:3])[N:5]([*:4])[*:6]',
        'real_ids': {2718, 271948, 271949, 273460, 273462, 273492, 273494, 273496, 273498}
    },
    {
        'reagents': [
            '[*:1][NH1:2][*:3]',
            '[O:4]=[S:5](=[O:6])([F,Cl,Br,I:7])[*:8]'
        ],
        'product': '[O:4]=[S:5](=[O:6])([*:8])[N:2]([*:1])[*:3]',
        'real_ids': {20, 40, 196680, 232682, 270084, 270188, 271082, 273578, 274078}
    },
    {
        'reagents': [
            '[O:1]=[C:2]([O:3])[*:4]'
            '[F,Cl,Br,I:5][*:6]'
        ],
        'product': '[O:1]=[C:2]([*:4])[O:3][*:6]',
        'real_ids': {1458}
    }
]

# Create reaction SMARTS from reagent and product SMARTS
for reaction in REAL_REACTIONS:
    # TODO: do we need parentheses for the reagents to force them to be separate molecules?
    # reaction = reagent_1.reagent_2...reagent_n>>product
    reaction['reaction'] = f'{".".join(reaction["reagent"])}>>{reaction["product"]}'
