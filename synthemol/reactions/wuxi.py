"""SMARTS representations of the WuXi GalaXi reactions."""
from synthemol.reactions.query_mol import QueryMol
from synthemol.reactions.reaction import Reaction

WUXI_REACTIONS_PHASE_1 = (
    # AF/AK/AG/AN + SA/SB/SC/SD/SE
    Reaction(
        reactants=[
            QueryMol('[N&!H0;!$(N[C,S]=[O,S,N]):4]'),  # AF/AK/AG/AN
            QueryMol('[O:2]=[C:1]([OH,O-])[*:3]')  # SA/SB/SC/SD/SE
        ],
        product=QueryMol('[O:2]=[C:1]([N:4])[*:3]'),
        reaction_id=1,
        chemical_space='WuXi'
    ),
    # AF/AK/AG/AN + HU
    Reaction(
        reactants=[
            QueryMol('[N&!H0;!$(N[C,S]=[O,S,N]):5]'),  # AF/AK/AG/AN
            QueryMol('Cl[S:1]([*:3])(=[O:2])=[O:4]'),  # HU
        ],
        product=QueryMol('[N:5][S:1]([*:3])(=[O:2])=[O:4]'),
        reaction_id=2,
        chemical_space='WuXi',
    ),
    # AF/AK/AG/AN + YV
    Reaction(
        reactants=[
            QueryMol('[#6:1]-[N&!H0;!$(N[C,S]=[O,S,N]):2]'),  # AF/AK/AG/AN
            QueryMol('[#6:3]-[#7]=C=O'),  # YV
        ],
        product=QueryMol('[#6:3]-[#7]-[#6](=O)-[#7:2]-[#6:1]'),
        reaction_id=3,
        chemical_space='WuXi',
    ),
    # LR + BT
    Reaction(
        reactants=[
            QueryMol('[c:2]-[#35]'),  # LR
            QueryMol('[#6:1]-[#5](-[#8])-[#8]'),  # BT
        ],
        product=QueryMol('[#6:1]-[#6:2]'),
        reaction_id=4,
        chemical_space='WuXi',
    ),
)

WUXI_REACTIONS_PHASE_2 = (
    # NH_N-Boc + SA/SB/SC/SD/SE + SA/SB/SC/SD/SE
    Reaction(
        reactants=[
            QueryMol('[N&!H0;!$(N[C,S]=[O,S,N]):1].CC(C)(C)OC(=O)[N:2]'),  # NH_N-Boc
            QueryMol('[O:4]=[C:3]([OH,O-])[*:5]'),  # SA/SB/SC/SD/SE
            QueryMol('[O:7]=[C:6]([OH,O-])[*:8]')  # SA/SB/SC/SD/SE
        ],
        product=QueryMol('[O:4]=[C:3]([N:1])[*:5].[O:7]=[C:6]([N:2])[*:8]'),
        reaction_id=5,
        chemical_space='WuXi'
    ),
    # COOH_N-Boc + AF/AG/AK/AN + SA/SB/SC/SD/SE
    Reaction(
        reactants=[
            QueryMol('[O:2]=[C:1]([OH,O-])[*:3].CC(C)(C)OC(=O)[N:4]'),  # COOH_N-Boc
            QueryMol('[N&!H0;!$(N[C,S]=[O,S,N]):5]'),  # AF/AG/AK/AN
            QueryMol('[O:7]=[C:6]([OH,O-])[*:8]')  # SA/SB/SC/SD/SE
        ],
        product=QueryMol('[O:2]=[C:1]([N:5])[*:3].[O:7]=[C:6]([N:4])[*:8]'),
        reaction_id=6,
        chemical_space='WuXi'
    ),
    # COOH_COOMe/COOH_COOEt + AF/AG/AK/AN + AF/AG/AK/AN
    Reaction(
        reactants=[
            QueryMol('[O:2]=[C:1]([OH,O-])[*:3].[O:5]=[C:4](O[C,CC])[*:6]'),  # COOH_COOMe/COOH_COOEt
            QueryMol('[N&!H0;!$(N[C,S]=[O,S,N]):7]'),  # AF/AG/AK/AN
            QueryMol('[N&!H0;!$(N[C,S]=[O,S,N]):8]')  # AF/AG/AK/AN
        ],
        product=QueryMol('[O:2]=[C:1]([N:7])[*:3].[O:5]=[C:4]([N:8])[*:6]'),
        reaction_id=7,
        chemical_space='WuXi',
    ),
)

WUXI_REACTIONS_PHASE_3 = (
    Reaction(
        reactants=[
            QueryMol(''),
            QueryMol(''),
            QueryMol('')
        ],
        product=QueryMol(''),
        reaction_id=8,
        chemical_space='WuXi'
    ),
)

WUXI_REACTIONS = WUXI_REACTIONS_PHASE_1 + WUXI_REACTIONS_PHASE_2 + WUXI_REACTIONS_PHASE_3
