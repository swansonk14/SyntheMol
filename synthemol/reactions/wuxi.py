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
    # BR_COOH/COOMe/COOEt + BT + AF/AG/AK/AN
    Reaction(
        reactants=[
            QueryMol('[c:1]-[#35].[O:3]=[C:2]([OH,O-])[*:4]'),  # BR_COOH/COOMe/COOEt
            QueryMol('[#6:5]-[#5](-[#8])-[#8]'),  # BT
            QueryMol('[N&!H0;!$(N[C,S]=[O,S,N]):6]')  # AF/AG/AK/AN
        ],
        product=QueryMol('[#6:5]-[#6:1].[O:3]=[C:2]([N:6])[*:4]'),
        reaction_id=8,
        chemical_space='WuXi'
    ),
    # BR_N-Boc + BT + SA/SB/SC/SD/SE
    Reaction(
        reactants=[
            QueryMol('[c:1]-[#35].CC(C)(C)OC(=O)[N:2]'),  # BR_N-Boc
            QueryMol('[#6:3]-[#5](-[#8])-[#8]'),  # BT
            QueryMol('[O:5]=[C:4]([OH,O-])[*:6]')  # SA/SB/SC/SD/SE
        ],
        product=QueryMol('[#6:3]-[#6:1].[O:5]=[C:4]([N:2])[*:6]'),
        reaction_id=9,
        chemical_space='WuXi'
    ),
    # BR_N-Boc + BT + YV
        Reaction(
        reactants=[
            QueryMol('[c:1]-[#35].CC(C)(C)OC(=O)[N:2]-[C:3]'),  # BR_N-Boc
            QueryMol('[#6:4]-[#5](-[#8])-[#8]'),  # BT
            QueryMol('[#6:5]-[#7]=C=O')  # YV
        ],
        product=QueryMol('[#6:4]-[#6:1].[#6:5]-[#7]-[#6](=O)-[#7:2]-[#6:3]'),
        reaction_id=10,
        chemical_space='WuXi'
    ),
    # BR_N-Boc + BT + LP
        Reaction(
        reactants=[
            QueryMol('[c:1]-[#35].CC(C)(C)OC(=O)[N:2]-[C:3]'),  # BR_N-Boc
            QueryMol('[#6:4]-[#5](-[#8])-[#8]'),  # BT
            QueryMol('[#6:5]-[#17]')  # LP
        ],
        product=QueryMol('[#6:4]-[#6:1].[#6:3]-[#7:2]-[#6:5]'),
        reaction_id=11,
        chemical_space='WuXi'
    ),
    # BR_N-Boc + BT + HU
        Reaction(
        reactants=[
            QueryMol('[c:1]-[#35].CC(C)(C)OC(=O)[N:2]'),  # BR_N-Boc
            QueryMol('[#6:3]-[#5](-[#8])-[#8]'),  # BT
            QueryMol('Cl[S:4]([*:6])(=[O:5])=[O:7]')  # HU
        ],
        product=QueryMol('[#6:3]-[#6:1].[N:2][S:4]([*:6])(=[O:5])=[O:7]'),
        reaction_id=12,
        chemical_space='WuXi'
    ),
    # BR_N-Boc + BT + QS
        Reaction(
        reactants=[
            QueryMol('[c:1]-[#35].CC(C)(C)OC(=O)[N:2]-[C:3]'),  # BR_N-Boc
            QueryMol('[#6:4]-[#5](-[#8])-[#8]'),  # BT
            QueryMol('[#6:5]-[#6;D2:6]=O')  # QS
        ],
        product=QueryMol('[#6:4]-[#6:1].[#6:5]-[#6:6]-[#7:2]-[#6:3]'),
        reaction_id=13,
        chemical_space='WuXi'
    ),
    # COOH_N-Boc + AF/AG/AK/AN + YV
        Reaction(
        reactants=[
            QueryMol('[O:2]=[C:1]([OH,O-])[*:3].CC(C)(C)OC(=O)[N:4]-[C:5]'),  # COOH_N-Boc
            QueryMol('[N&!H0;!$(N[C,S]=[O,S,N]):6]'),  # AF/AG/AK/AN
            QueryMol('[#6:7]-[#7]=C=O')  # YV
        ],
        product=QueryMol('[O:2]=[C:1]([N:6])[*:3].[#6:7]-[#7]-[#6](=O)-[#7:4]-[#6:5]'),
        reaction_id=14,
        chemical_space='WuXi'
    ),
    # COOH_N-Boc + AF/AG/AK/AN + LP
        Reaction(
        reactants=[
            QueryMol('[O:2]=[C:1]([OH,O-])[*:3].CC(C)(C)OC(=O)[N:4]-[C:5]'),  # COOH_N-Boc
            QueryMol('[N&!H0;!$(N[C,S]=[O,S,N]):6]'),  # AF/AG/AK/AN
            QueryMol('[#6:7]-[#17]')  # LP
        ],
        product=QueryMol('[O:2]=[C:1]([N:6])[*:3].[#6:5]-[#7:4]-[#6:7]'),
        reaction_id=15,
        chemical_space='WuXi'
    ),
    # COOH_N-Boc + AF/AG/AK/AN + HU
        Reaction(
        reactants=[
            QueryMol('[O:2]=[C:1]([OH,O-])[*:3].CC(C)(C)OC(=O)[N:4]'),  # COOH_N-Boc
            QueryMol('[N&!H0;!$(N[C,S]=[O,S,N]):5]'),  # AF/AG/AK/AN
            QueryMol('Cl[S:6]([*:8])(=[O:7])=[O:9]')  # HU
        ],
        product=QueryMol('[O:2]=[C:1]([N:5])[*:3].[N:4][S:6]([*:8])(=[O:7])=[O:9]'),
        reaction_id=16,
        chemical_space='WuXi'
    ),
    # COOH_N-Boc + AF/AG/AK/AN + QS
        Reaction(
        reactants=[
            QueryMol('[O:2]=[C:1]([OH,O-])[*:3].CC(C)(C)OC(=O)[N:4]-[C:5]'),  # COOH_N-Boc
            QueryMol('[N&!H0;!$(N[C,S]=[O,S,N]):6]'),  # AF/AG/AK/AN
            QueryMol('[#6:7]-[#6;D2:8]=O')  # QS
        ],
        product=QueryMol('[O:2]=[C:1]([N:6])[*:3].[#6:7]-[#6:8]-[#7:4]-[#6:5]'),
        reaction_id=17,
        chemical_space='WuXi'
    ),
    # NH_N-Boc + YV + SA/SB/SC/SD/SE
        Reaction(
        reactants=[
            QueryMol(''),  # NH_N-Boc
            QueryMol(''),  # YV
            QueryMol('')  # SA/SB/SC/SD/SE
        ],
        product=QueryMol(''),
        reaction_id=18,
        chemical_space='WuXi'
    ),
    # NH_N-Boc + YV + YV
        Reaction(
        reactants=[
            QueryMol(''),  # NH_N-Boc
            QueryMol(''),  # YV
            QueryMol('')  # YV
        ],
        product=QueryMol(''),
        reaction_id=19,
        chemical_space='WuXi'
    ),
    # NH_N-Boc + YV + LP
        Reaction(
        reactants=[
            QueryMol(''),  # NH_N-Boc
            QueryMol(''),  # YV
            QueryMol('')  # LP
        ],
        product=QueryMol(''),
        reaction_id=20,
        chemical_space='WuXi'
    ),
    # NH_N-Boc + YV + HU
        Reaction(
        reactants=[
            QueryMol(''),  # NH_N-Boc
            QueryMol(''),  # YV
            QueryMol('')  # HU
        ],
        product=QueryMol(''),
        reaction_id=21,
        chemical_space='WuXi'
    ),
    # NH_N-Boc + YV + QS
        Reaction(
        reactants=[
            QueryMol(''),  # NH_N-Boc
            QueryMol(''),  # YV
            QueryMol('')  # QS
        ],
        product=QueryMol(''),
        reaction_id=22,
        chemical_space='WuXi'
    ),
    # NH_N-Boc + LP + SA/SB/SC/SD/SE
        Reaction(
        reactants=[
            QueryMol(''),  # NH_N-Boc
            QueryMol(''),  # LP
            QueryMol('')  # SA/SB/SC/SD/SE
        ],
        product=QueryMol(''),
        reaction_id=23,
        chemical_space='WuXi'
    ),
    # NH_N-Boc + LP + YV
        Reaction(
        reactants=[
            QueryMol(''),  # NH_N-Boc
            QueryMol(''),  # LP
            QueryMol('')  # YV
        ],
        product=QueryMol(''),
        reaction_id=24,
        chemical_space='WuXi'
    ),
    # NH_N-Boc + LP + LP
        Reaction(
        reactants=[
            QueryMol(''),  # NH_N-Boc
            QueryMol(''),  # LP
            QueryMol('')  # LP
        ],
        product=QueryMol(''),
        reaction_id=25,
        chemical_space='WuXi'
    ),
    # NH_N-Boc + LP + HU
        Reaction(
        reactants=[
            QueryMol(''),  # NH_N-Boc
            QueryMol(''),  # LP
            QueryMol('')  # HU
        ],
        product=QueryMol(''),
        reaction_id=26,
        chemical_space='WuXi'
    ),
    # NH_N-Boc + LP + QS
        Reaction(
        reactants=[
            QueryMol(''),  # NH_N-Boc
            QueryMol(''),  # LP
            QueryMol('')  # QS
        ],
        product=QueryMol(''),
        reaction_id=27,
        chemical_space='WuXi'
    ),
    # NH_N-Boc + HU + SA/SB/SC/SD/SE
        Reaction(
        reactants=[
            QueryMol(''),  # NH_N-Boc
            QueryMol(''),  # HU
            QueryMol('')  # SA/SB/SC/SD/SE
        ],
        product=QueryMol(''),
        reaction_id=28,
        chemical_space='WuXi'
    ),
    # NH_N-Boc + HU + YV
        Reaction(
        reactants=[
            QueryMol(''),  # NH_N-Boc
            QueryMol(''),  # HU
            QueryMol('')  # YV
        ],
        product=QueryMol(''),
        reaction_id=29,
        chemical_space='WuXi'
    ),
    # NH_N-Boc + HU + LP
        Reaction(
        reactants=[
            QueryMol(''),  # NH_N-Boc
            QueryMol(''),  # HU
            QueryMol('')  # LP
        ],
        product=QueryMol(''),
        reaction_id=30,
        chemical_space='WuXi'
    ),
    # NH_N-Boc + HU + HU
        Reaction(
        reactants=[
            QueryMol(''),  # NH_N-Boc
            QueryMol(''),  # HU
            QueryMol('')  # HU
        ],
        product=QueryMol(''),
        reaction_id=31,
        chemical_space='WuXi'
    ),
    # NH_N-Boc + HU + QS
        Reaction(
        reactants=[
            QueryMol(''),  # NH_N-Boc
            QueryMol(''),  # HU
            QueryMol('')  # QS
        ],
        product=QueryMol(''),
        reaction_id=32,
        chemical_space='WuXi'
    ),
    # NH_N-Boc + SA/SB/SC/SD/SE + YV
        Reaction(
        reactants=[
            QueryMol(''),  # NH_N-Boc
            QueryMol(''),  # SA/SB/SC/SD/SE
            QueryMol('')  # YV
        ],
        product=QueryMol(''),
        reaction_id=33,
        chemical_space='WuXi'
    ),
    # NH_N-Boc + SA/SB/SC/SD/SE + LP
        Reaction(
        reactants=[
            QueryMol(''),  # NH_N-Boc
            QueryMol(''),  # SA/SB/SC/SD/SE
            QueryMol('')  # LP
        ],
        product=QueryMol(''),
        reaction_id=34,
        chemical_space='WuXi'
    ),
    # NH_N-Boc + SA/SB/SC/SD/SE + HU
        Reaction(
        reactants=[
            QueryMol(''),  # NH_N-Boc
            QueryMol(''),  # SA/SB/SC/SD/SE
            QueryMol('')  # HU
        ],
        product=QueryMol(''),
        reaction_id=35,
        chemical_space='WuXi'
    ),
    # NH_N-Boc + SA/SB/SC/SD/SE + QS
        Reaction(
        reactants=[
            QueryMol(''),  # NH_N-Boc
            QueryMol(''),  # SA/SB/SC/SD/SE
            QueryMol('')  # QS
        ],
        product=QueryMol(''),
        reaction_id=36,
        chemical_space='WuXi'
    ),
)

WUXI_REACTIONS = WUXI_REACTIONS_PHASE_1 + WUXI_REACTIONS_PHASE_2 + WUXI_REACTIONS_PHASE_3
