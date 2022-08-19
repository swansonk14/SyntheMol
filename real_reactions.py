"""SMARTS representations of the REAL reactions.
Reference: https://docs.google.com/document/d/1LDgRXf4P-uOXQEmgJPgVhuOK32I2u0FXSB8tzDenN6U/edit?usp=sharing
"""
from reactions import (
    alkyl_checker,
    aryl_checker,
    CarbonPathChecker,
    count_three_reagents_with_two_same,
    count_two_different_reagents,
    count_two_same_reagents,
    cycle_checker,
    h_checker,
    QueryMol,
    Reaction,
    RGroupChecker
)

# TODO: for all these reactions, they may match more REAL reactions than the ones noted
# TODO: do product QueryMols need checkers?
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
        product=QueryMol('[*:4][N:5]([*:6])[C:7](=[O:9])[*:8].[*:3][N:1]([*:2])[C:10](=[O:12])[*:11]'),
        reaction_id=1,
        real_ids={275592},
        counting_fn=count_three_reagents_with_two_same
    ),
    Reaction(
        reagents=[
            QueryMol('[*:1][N:2]([H])[*:3]'),
            QueryMol('[OH1][C:4]([*:5])=[O:6]')
        ],
        product=QueryMol('[*:5][C:4](=[O:6])[N:2]([*:1])[*:3]'),
        reaction_id=2,
        # Note: 240690 requires Boc-protected amino group in side chain and all use different catalysts
        real_ids={22, 11, 527, 240690},
        counting_fn=count_two_different_reagents
    ),
    Reaction(
        reagents=[
            QueryMol(
                smarts='[*:1][N:2]([H])[H:3]',
                checker_class=RGroupChecker,
                checker_kwargs={'checkers': {aryl_checker}}
            ),
            QueryMol(
                smarts='[*:4][N:5]([H])[*:6]',
                checker_class=RGroupChecker,
                checker_kwargs={'checkers': {alkyl_checker, h_checker}}
            )
        ],
        product=QueryMol('O=C([N:2]([*:1])[H:3])[N:5]([*:4])[*:6]'),
        reaction_id=3,
        real_ids={2430},
        counting_fn=count_two_same_reagents
    ),
    Reaction(
        reagents=[
            QueryMol('[*:1][N:2]([H])[H:3]'),
            QueryMol('[*:4][N:5]([H])[H:6]')
        ],
        product=QueryMol('O=C([N:2]([*:1])[H:3])[N:5]([*:4])[H:6]'),
        reaction_id=4,
        real_ids={2708},
        counting_fn=count_two_same_reagents
    ),
    Reaction(
        reagents=[
            QueryMol(
                smarts='[*:1][N:2]([H])[*:3]',
                checker_class=RGroupChecker,
                checker_kwargs={'checkers': {alkyl_checker, aryl_checker, cycle_checker}}
            ),
            QueryMol(
                smarts='[F,Cl,Br,I][*:4]',
                checker_class=RGroupChecker,
                checker_kwargs={'checkers': {alkyl_checker}}
            )
        ],
        product=QueryMol('[*:1][N:2]([*:3])[*:4]'),
        reaction_id=5,
        real_ids={2230},
        counting_fn=count_two_different_reagents
    ),
    Reaction(
        reagents=[
            QueryMol(
                smarts='[*:1][N:2]([H])[H:3]',
                checker_class=RGroupChecker,
                checker_kwargs={'checkers': {aryl_checker}}
            ),
            QueryMol(
                smarts='[*:4][N:5]([H])[H:6]',
                checker_class=RGroupChecker,
                checker_kwargs={'checkers': {alkyl_checker}}
            )
        ],
        product=QueryMol('O=C(C(=O)[N:2]([*:1])[H:3])[N:5]([*:4])[H:6]'),
        reaction_id=6,
        real_ids={2718},
        counting_fn=count_two_same_reagents
    ),
    Reaction(
        reagents=[
            QueryMol(
                smarts='[*:1][N:2]([H])[*:3]',
                checker_class=RGroupChecker,
                checker_kwargs={'checkers': {alkyl_checker, h_checker}}
            ),
            QueryMol('[O:4]=[S:5](=[O:6])([F,Cl,Br,I])[*:7]')
        ],
        product=QueryMol('[O:4]=[S:5](=[O:6])([*:7])[N:2]([*:1])[*:3]'),
        reaction_id=7,
        real_ids={40},
        counting_fn=count_two_different_reagents
    ),
    Reaction(
        reagents=[
            QueryMol('[OH1:1][C:2]([*:3])=[O:4]'),
            QueryMol(
                smarts='[F,Cl,Br,I][*:5]',
                checker_class=RGroupChecker,
                checker_kwargs={'checkers': {alkyl_checker}}
            )
        ],
        product=QueryMol('[O:4]=[C:2]([*:3])[O:1][*:5]'),
        reaction_id=8,
        real_ids={1458},
        counting_fn=count_two_different_reagents
    ),
    Reaction(
        reagents=[
            QueryMol(
                smarts='[*:1][N:2]([H])[*:3]',
                checker_class=RGroupChecker,
                checker_kwargs={'checkers': {alkyl_checker}}
            ),
            QueryMol('[*:4][N:5]([H])[H:6]')
        ],
        product=QueryMol('O=C(C(=O)[N:2]([*:1])[*:3])[N:5]([*:4])[H:6]'),
        reaction_id=9,
        real_ids={2718},
        counting_fn=count_two_same_reagents
    ),
    Reaction(
        reagents=[
            QueryMol(
                smarts='[*:1][N:2]([H])[*:3]',
                checker_class=RGroupChecker,
                checker_kwargs={'checkers': {alkyl_checker}}
            ),
            QueryMol(
                smarts='[F,Cl,Br,I][*:4]',
                checker_class=RGroupChecker,
                checker_kwargs={'checkers': {aryl_checker}}
            )
        ],
        product=QueryMol('[*:1][N:2]([*:3])[*:4]'),
        reaction_id=10,
        real_ids={27},
        counting_fn=count_two_different_reagents
    )
]

# Check that reaction IDs are unique and count from 1 to the number of reactions
reaction_ids = {reaction.id for reaction in REAL_REACTIONS}

if reaction_ids != set(range(1, len(REAL_REACTIONS) + 1)):
    raise ValueError('REAL reaction IDs are not unique and/or do not count from 1 to the number of reactions.')
