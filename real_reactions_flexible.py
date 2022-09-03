"""SMARTS representations of the **generalized** most common REAL reactions.
Reference: https://docs.google.com/document/d/1jEFqvKT4bFe4Il90Pg53AHCtQQf5eA46RST4ZWQJhyI/edit?usp=sharing
"""
from reactions import (
    CarbonPathChecker,
    count_three_reagents_with_two_same,
    count_two_different_reagents,
    count_two_same_reagents,
    QueryMol,
    Reaction
)


REAL_REACTIONS = [
    Reaction(
        reagents=[
            QueryMol(
                smarts='CC(C)(C)OC(=O)[N:1]([*:2])[*:3].[*:4][N:5]([H])[*:6]',
                checker_class=CarbonPathChecker
            ),
            QueryMol('[OH1][C:7]([*:8])=[O:9]'),
            QueryMol('[OH1][C:10]([*:11])=[O:12]')
        ],
        product=QueryMol(
            smarts='[*:4][N:5]([*:6])[C:7](=[O:9])[*:8].[*:3][N:1]([*:2])[C:10](=[O:12])[*:11]',
            checker_class=CarbonPathChecker
        ),
        reaction_id=1,
        real_ids={275592},
        counting_fn=count_three_reagents_with_two_same
    ),
    Reaction(
        reagents=[
            QueryMol('[*:1][N:2]([H])[*:3]'),
            QueryMol('[F,Cl,Br,I][*:4]')
        ],
        product=QueryMol('[*:1][N:2]([*:3])[*:4]'),
        reaction_id=2,
        real_ids={7, 27, 34, 38, 44, 61, 2230, 269956, 269982, 270122, 270166, 270344, 272692, 272710, 273654},
        counting_fn=count_two_different_reagents
    ),
    Reaction(
        reagents=[
            QueryMol('[*:1][OH1,SH1:2]'),
            QueryMol('[F,Cl,Br,I][*:3]')
        ],
        product=QueryMol('[*:1][O,S:2][*:3]'),
        reaction_id=3,
        real_ids={7, 34, 272692, 272710, 273654},
        counting_fn=count_two_different_reagents
    ),
    Reaction(
        reagents=[
            QueryMol('[*:1][N:2]([H])[*:3]'),
            QueryMol('[*:4][N:5]([H])[*:6]')
        ],
        product=QueryMol('O=C([N:2]([*:1])[*:3])[N:5]([*:4])[*:6]'),
        reaction_id=4,
        real_ids={512, 2430, 2554, 2708, 272164, 273390, 273392, 273452, 273454, 273458, 273464, 273574},
        counting_fn=count_two_same_reagents
    ),
    Reaction(
        reagents=[
            QueryMol('[*:1][N:2]([H])[*:3]'),
            QueryMol('[*:4][N:5]([H])[*:6]')
        ],
        product=QueryMol('O=C(C(=O)[N:2]([*:1])[*:3])[N:5]([*:4])[*:6]'),
        reaction_id=5,
        real_ids={2718, 271948, 271949, 273460, 273462, 273492, 273494, 273496, 273498},
        counting_fn=count_two_same_reagents
    ),
    Reaction(
        reagents=[
            QueryMol('[*:1][N:2]([H])[*:3]'),
            QueryMol('[OH1][C:4]([*:5])=[O:6]')
        ],
        product=QueryMol('[*:5][C:4](=[O:6])[N:2]([*:1])[*:3]'),
        reaction_id=6,
        real_ids={11, 22, 527, 188690, 240690, 269946, 270006, 270112, 272270, 272430, 273450, 273652, 280130},
        counting_fn=count_two_different_reagents
    ),
    Reaction(
        reagents=[
            QueryMol('[O:1]=[C:2]([OH1:3])[*:4]'),
            QueryMol('[F,Cl,Br,I][*:5]')
        ],
        product=QueryMol('[O:1]=[C:2]([*:4])[O:3][*:5]'),
        reaction_id=7,
        real_ids={1458},
        counting_fn=count_two_different_reagents
    ),
    Reaction(
        reagents=[
            QueryMol('[*:1][N:2]([H])[*:3]'),
            QueryMol('[O:4]=[S:5](=[O:6])([F,Cl,Br,I])[*:7]')
        ],
        product=QueryMol('[O:4]=[S:5](=[O:6])([*:7])[N:2]([*:1])[*:3]'),
        reaction_id=8,
        real_ids={20, 40, 196680, 232682, 270084, 270188, 271082, 273578, 274078},
        counting_fn=count_two_different_reagents
    )
]

# Check that reaction IDs are unique and count from 1 to the number of reactions
reaction_ids = {reaction.id for reaction in REAL_REACTIONS}

if reaction_ids != set(range(1, len(REAL_REACTIONS) + 1)):
    raise ValueError('REAL reaction IDs are not unique and/or do not count from 1 to the number of reactions.')

REACTION_ID_TO_REAL_IDS = {
    reaction.id: reaction.real_ids
    for reaction in REAL_REACTIONS
}
