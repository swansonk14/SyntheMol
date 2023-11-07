"""SMARTS representations of the WuXi GalaXi reactions."""
from synthemol.reactions.query_mol import QueryMol
from synthemol.reactions.reaction import Reaction

A_SUBSETS = ("AF", "AK", "AG", "AN")
S_SUBSETS = ("SA", "SB", "SC", "SD", "SE")
BT_SUBSETS = ("BT",)
HU_SUBSETS = ("HU",)
LP_SUBSETS = ("LP",)
LR_SUBSETS = ("LR",)
QS_SUBSETS = ("QS",)
YV_SUBSETS = ("YV",)
PHASE_2_NH_N_BOC_SUBSETS = ("phase_2_NH_N-Boc",)
PHASE_2_COOH_N_BOC_SUBSETS = ("phase_2_COOH_N-Boc",)
PHASE_2_COOH_COOMeEt_SUBSETS = ("phase_2_COOH_COOMe", "phase_2_COOH_COOEt")
PHASE_3_BR_COOHMeEt_SUBSETS = (
    "phase_3_Br_COOH",
    "phase_3_Br_COOHMe",
    "phase_3_Br_COOHEt",
)
PHASE_3_BR_N_BOC_SUBSETS = ("phase_3_Br_N-Boc",)
PHASE_3_COOH_N_BOC_SUBSETS = ("phase_3_COOH_N-Boc",)
PHASE_3_N_BOC_NH_SUBSETS = ("phase_3_N-Boc_NH",)


WUXI_REACTIONS_PHASE_1 = (
    # AF/AK/AG/AN + SA/SB/SC/SD/SE
    Reaction(
        reactants=[
            QueryMol("[N&!H0;!$(N[C,S]=[O,S,N]):11]", subsets=A_SUBSETS),  # AF/AK/AG/AN
            QueryMol(
                "[O:22]=[C:21]([OH,O-])[*:23]", subsets=S_SUBSETS
            ),  # SA/SB/SC/SD/SE
        ],
        product=QueryMol("[O:22]=[C:21]([N:11])[*:23]"),
        reaction_id=1,
        chemical_space="wuxi",
    ),
    # AF/AK/AG/AN + HU
    Reaction(
        reactants=[
            QueryMol("[N&!H0;!$(N[C,S]=[O,S,N]):11]", subsets=A_SUBSETS),  # AF/AK/AG/AN
            QueryMol("Cl[S:21]([*:23])(=[O:22])=[O:24]", subsets=HU_SUBSETS),  # HU
        ],
        product=QueryMol("[N:11][S:21]([*:23])(=[O:22])=[O:24]"),
        reaction_id=2,
        chemical_space="wuxi",
    ),
    # AF/AK/AG/AN + YV
    Reaction(
        reactants=[
            QueryMol(
                "[#6:11]-[N&!H0;!$(N[C,S]=[O,S,N]):12]", subsets=A_SUBSETS
            ),  # AF/AK/AG/AN
            QueryMol("[#6:21]-[#7]=C=O", subsets=YV_SUBSETS),  # YV
        ],
        product=QueryMol("[#6:21]-[#7]-[#6](=O)-[#7:12]-[#6:11]"),
        reaction_id=3,
        chemical_space="wuxi",
    ),
    # LR + BT
    Reaction(
        reactants=[
            QueryMol("[c:11]-[#35]", subsets=LR_SUBSETS),  # LR
            QueryMol("[#6:21]-[#5](-[#8])-[#8]", subsets=BT_SUBSETS),  # BT
        ],
        product=QueryMol("[#6:21]-[#6:11]"),
        reaction_id=4,
        chemical_space="wuxi",
    ),
)

WUXI_REACTIONS_PHASE_2 = (
    # NH_N-Boc + SA/SB/SC/SD/SE + SA/SB/SC/SD/SE
    Reaction(
        reactants=[
            QueryMol(
                "[N&!H0;!$(N[C,S]=[O,S,N]):11].CC(C)(C)OC(=O)[N:12]",
                subsets=PHASE_2_NH_N_BOC_SUBSETS,
            ),  # NH_N-Boc
            QueryMol(
                "[O:22]=[C:21]([OH,O-])[*:23]", subsets=S_SUBSETS
            ),  # SA/SB/SC/SD/SE
            QueryMol(
                "[O:32]=[C:31]([OH,O-])[*:33]", subsets=S_SUBSETS
            ),  # SA/SB/SC/SD/SE
        ],
        product=QueryMol("[O:22]=[C:21]([N:11])[*:23].[O:32]=[C:31]([N:12])[*:33]"),
        reaction_id=5,
        chemical_space="wuxi",
    ),
    # COOH_N-Boc + AF/AG/AK/AN + SA/SB/SC/SD/SE
    Reaction(
        reactants=[
            QueryMol(
                "[O:12]=[C:11]([OH,O-])[*:13].CC(C)(C)OC(=O)[N:14]",
                subsets=PHASE_2_COOH_N_BOC_SUBSETS,
            ),  # COOH_N-Boc
            QueryMol("[N&!H0;!$(N[C,S]=[O,S,N]):21]", subsets=A_SUBSETS),  # AF/AG/AK/AN
            QueryMol(
                "[O:32]=[C:31]([OH,O-])[*:33]", subsets=S_SUBSETS
            ),  # SA/SB/SC/SD/SE
        ],
        product=QueryMol("[O:12]=[C:11]([N:21])[*:13].[O:32]=[C:31]([N:14])[*:33]"),
        reaction_id=6,
        chemical_space="wuxi",
    ),
    # COOH_COOMe/COOH_COOEt + AF/AG/AK/AN + AF/AG/AK/AN
    Reaction(
        reactants=[
            QueryMol(
                "[O:12]=[C:11]([OH,O-])[*:13].[O:15]=[C:14](O[C,CC])[*:16]",
                subsets=PHASE_2_COOH_COOMeEt_SUBSETS,
            ),  # COOH_COOMe/COOH_COOEt
            QueryMol("[N&!H0;!$(N[C,S]=[O,S,N]):21]", subsets=A_SUBSETS),  # AF/AG/AK/AN
            QueryMol("[N&!H0;!$(N[C,S]=[O,S,N]):31]", subsets=A_SUBSETS),  # AF/AG/AK/AN
        ],
        product=QueryMol("[O:12]=[C:11]([N:21])[*:13].[O:15]=[C:14]([N:31])[*:16]"),
        reaction_id=7,
        chemical_space="wuxi",
    ),
)

WUXI_REACTIONS_PHASE_3 = (
    # BR_COOH/COOMe/COOEt + BT + AF/AG/AK/AN
    Reaction(
        reactants=[
            QueryMol(
                "[c:11]-[#35].[O:13]=[C:12]([$([OH,O-]),$(O[C,CC])])[*:14]",
                subsets=PHASE_3_BR_COOHMeEt_SUBSETS,
            ),  # BR_COOH/COOMe/COOEt
            QueryMol("[#6:21]-[#5](-[#8])-[#8]", subsets=BT_SUBSETS),  # BT
            QueryMol("[N&!H0;!$(N[C,S]=[O,S,N]):31]", subsets=A_SUBSETS),  # AF/AG/AK/AN
        ],
        product=QueryMol("[#6:21]-[#6:11].[O:13]=[C:12]([N:31])[*:14]"),
        reaction_id=8,
        chemical_space="wuxi",
    ),
    # BR_N-Boc + BT + SA/SB/SC/SD/SE
    Reaction(
        reactants=[
            QueryMol(
                "[c:11]-[#35].CC(C)(C)OC(=O)[N:12]", subsets=PHASE_3_BR_N_BOC_SUBSETS
            ),  # BR_N-Boc
            QueryMol("[#6:21]-[#5](-[#8])-[#8]", subsets=BT_SUBSETS),  # BT
            QueryMol(
                "[O:32]=[C:31]([OH,O-])[*:33]", subsets=S_SUBSETS
            ),  # SA/SB/SC/SD/SE
        ],
        product=QueryMol("[#6:21]-[#6:11].[O:32]=[C:31]([N:12])[*:33]"),
        reaction_id=9,
        chemical_space="wuxi",
    ),
    # BR_N-Boc + BT + YV
    Reaction(
        reactants=[
            QueryMol(
                "[c:11]-[#35].CC(C)(C)OC(=O)[N:12]-[C:13]",
                subsets=PHASE_3_BR_N_BOC_SUBSETS,
            ),  # BR_N-Boc
            QueryMol("[#6:21]-[#5](-[#8])-[#8]", subsets=BT_SUBSETS),  # BT
            QueryMol("[#6:31]-[#7]=C=O", subsets=YV_SUBSETS),  # YV
        ],
        product=QueryMol("[#6:21]-[#6:11].[#6:31]-[#7]-[#6](=O)-[#7:12]-[#6:13]"),
        reaction_id=10,
        chemical_space="wuxi",
    ),
    # BR_N-Boc + BT + LP
    Reaction(
        reactants=[
            QueryMol(
                "[c:11]-[#35].CC(C)(C)OC(=O)[N:12]-[C:13]",
                subsets=PHASE_3_BR_N_BOC_SUBSETS,
            ),  # BR_N-Boc
            QueryMol("[#6:21]-[#5](-[#8])-[#8]", subsets=BT_SUBSETS),  # BT
            QueryMol("[#6:31]-[#17]", subsets=LP_SUBSETS),  # LP
        ],
        product=QueryMol("[#6:21]-[#6:11].[#6:13]-[#7:12]-[#6:31]"),
        reaction_id=11,
        chemical_space="wuxi",
    ),
    # BR_N-Boc + BT + HU
    Reaction(
        reactants=[
            QueryMol(
                "[c:11]-[#35].CC(C)(C)OC(=O)[N:12]", subsets=PHASE_3_BR_N_BOC_SUBSETS
            ),  # BR_N-Boc
            QueryMol("[#6:21]-[#5](-[#8])-[#8]", subsets=BT_SUBSETS),  # BT
            QueryMol("Cl[S:31]([*:33])(=[O:32])=[O:34]", subsets=HU_SUBSETS),  # HU
        ],
        product=QueryMol("[#6:21]-[#6:11].[N:12][S:31]([*:33])(=[O:32])=[O:34]"),
        reaction_id=12,
        chemical_space="wuxi",
    ),
    # BR_N-Boc + BT + QS
    Reaction(
        reactants=[
            QueryMol(
                "[c:11]-[#35].CC(C)(C)OC(=O)[N:12]-[C:13]",
                subsets=PHASE_3_BR_N_BOC_SUBSETS,
            ),  # BR_N-Boc
            QueryMol("[#6:21]-[#5](-[#8])-[#8]", subsets=BT_SUBSETS),  # BT
            QueryMol("[#6:31]-[#6;D3:32]=O", subsets=QS_SUBSETS),  # QS
        ],
        product=QueryMol("[#6:21]-[#6:11].[#6:31]-[#6:32]-[#7:12]-[#6:13]"),
        reaction_id=13,
        chemical_space="wuxi",
    ),
    # COOH_N-Boc + AF/AG/AK/AN + YV
    Reaction(
        reactants=[
            QueryMol(
                "[O:12]=[C:11]([OH,O-])[*:13].CC(C)(C)OC(=O)[N:14]-[C:15]",
                subsets=PHASE_3_COOH_N_BOC_SUBSETS,
            ),  # COOH_N-Boc
            QueryMol("[N&!H0;!$(N[C,S]=[O,S,N]):21]", subsets=A_SUBSETS),  # AF/AG/AK/AN
            QueryMol("[#6:31]-[#7]=C=O", subsets=YV_SUBSETS),  # YV
        ],
        product=QueryMol(
            "[O:12]=[C:11]([N:21])[*:13].[#6:31]-[#7]-[#6](=O)-[#7:14]-[#6:15]"
        ),
        reaction_id=14,
        chemical_space="wuxi",
    ),
    # COOH_N-Boc + AF/AG/AK/AN + LP
    Reaction(
        reactants=[
            QueryMol(
                "[O:12]=[C:11]([OH,O-])[*:13].CC(C)(C)OC(=O)[N:14]-[C:15]",
                subsets=PHASE_3_COOH_N_BOC_SUBSETS,
            ),  # COOH_N-Boc
            QueryMol("[N&!H0;!$(N[C,S]=[O,S,N]):21]", subsets=A_SUBSETS),  # AF/AG/AK/AN
            QueryMol("[#6:31]-[#17]", subsets=LP_SUBSETS),  # LP
        ],
        product=QueryMol("[O:12]=[C:11]([N:21])[*:13].[#6:15]-[#7:14]-[#6:31]"),
        reaction_id=15,
        chemical_space="wuxi",
    ),
    # COOH_N-Boc + AF/AG/AK/AN + HU
    Reaction(
        reactants=[
            QueryMol(
                "[O:12]=[C:11]([OH,O-])[*:13].CC(C)(C)OC(=O)[N:14]",
                subsets=PHASE_3_COOH_N_BOC_SUBSETS,
            ),  # COOH_N-Boc
            QueryMol("[N&!H0;!$(N[C,S]=[O,S,N]):21]", subsets=A_SUBSETS),  # AF/AG/AK/AN
            QueryMol("Cl[S:31]([*:33])(=[O:32])=[O:34]", subsets=HU_SUBSETS),  # HU
        ],
        product=QueryMol(
            "[O:12]=[C:11]([N:21])[*:13].[N:14][S:31]([*:33])(=[O:32])=[O:34]"
        ),
        reaction_id=16,
        chemical_space="wuxi",
    ),
    # COOH_N-Boc + AF/AG/AK/AN + QS
    Reaction(
        reactants=[
            QueryMol(
                "[O:12]=[C:11]([OH,O-])[*:13].CC(C)(C)OC(=O)[N:14]-[C:15]",
                subsets=PHASE_3_COOH_N_BOC_SUBSETS,
            ),  # COOH_N-Boc
            QueryMol("[N&!H0;!$(N[C,S]=[O,S,N]):21]", subsets=A_SUBSETS),  # AF/AG/AK/AN
            QueryMol("[#6:31]-[#6;D3:32]=O", subsets=QS_SUBSETS),  # QS
        ],
        product=QueryMol("[O:12]=[C:11]([N:21])[*:13].[#6:31]-[#6:32]-[#7:14]-[#6:15]"),
        reaction_id=17,
        chemical_space="wuxi",
    ),
    # NH_N-Boc + YV + SA/SB/SC/SD/SE
    Reaction(
        reactants=[
            QueryMol(
                "[#6:11]-[N&!H0;!$(N[C,S]=[O,S,N]):12].CC(C)(C)OC(=O)[N:13]",
                subsets=PHASE_3_N_BOC_NH_SUBSETS,
            ),  # NH_N-Boc
            QueryMol("[#6:21]-[#7]=C=O", subsets=YV_SUBSETS),  # YV
            QueryMol(
                "[O:32]=[C:31]([OH,O-])[*:33]", subsets=S_SUBSETS
            ),  # SA/SB/SC/SD/SE
        ],
        product=QueryMol(
            "[#6:21]-[#7]-[#6](=O)-[#7:12]-[#6:11].[O:32]=[C:31]([N:13])[*:33]"
        ),
        reaction_id=18,
        chemical_space="wuxi",
    ),
    # NH_N-Boc + YV + YV
    Reaction(
        reactants=[
            QueryMol(
                "[#6:11]-[N&!H0;!$(N[C,S]=[O,S,N]):12].CC(C)(C)OC(=O)[N:13]-[C:14]",
                subsets=PHASE_3_N_BOC_NH_SUBSETS,
            ),  # NH_N-Boc
            QueryMol("[#6:21]-[#7]=C=O", subsets=YV_SUBSETS),  # YV
            QueryMol("[#6:31]-[#7]=C=O", subsets=YV_SUBSETS),  # YV
        ],
        product=QueryMol(
            "[#6:21]-[#7]-[#6](=O)-[#7:12]-[#6:11].[#6:31]-[#7]-[#6](=O)-[#7:13]-[#6:14]"
        ),
        reaction_id=19,
        chemical_space="wuxi",
    ),
    # NH_N-Boc + YV + LP
    Reaction(
        reactants=[
            QueryMol(
                "[#6:11]-[N&!H0;!$(N[C,S]=[O,S,N]):12].CC(C)(C)OC(=O)[N:13]-[C:14]",
                subsets=PHASE_3_N_BOC_NH_SUBSETS,
            ),  # NH_N-Boc
            QueryMol("[#6:21]-[#7]=C=O", subsets=YV_SUBSETS),  # YV
            QueryMol("[#6:31]-[#17]", subsets=LP_SUBSETS),  # LP
        ],
        product=QueryMol(
            "[#6:21]-[#7]-[#6](=O)-[#7:12]-[#6:11].[#6:14]-[#7:13]-[#6:31]"
        ),
        reaction_id=20,
        chemical_space="wuxi",
    ),
    # NH_N-Boc + YV + HU
    Reaction(
        reactants=[
            QueryMol(
                "[#6:11]-[N&!H0;!$(N[C,S]=[O,S,N]):12].CC(C)(C)OC(=O)[N:13]",
                subsets=PHASE_3_N_BOC_NH_SUBSETS,
            ),  # NH_N-Boc
            QueryMol("[#6:21]-[#7]=C=O", subsets=YV_SUBSETS),  # YV
            QueryMol("Cl[S:31]([*:33])(=[O:32])=[O:34]", subsets=HU_SUBSETS),  # HU
        ],
        product=QueryMol(
            "[#6:21]-[#7]-[#6](=O)-[#7:12]-[#6:11].[N:13][S:31]([*:33])(=[O:32])=[O:34]"
        ),
        reaction_id=21,
        chemical_space="wuxi",
    ),
    # NH_N-Boc + YV + QS
    Reaction(
        reactants=[
            QueryMol(
                "[#6:11]-[N&!H0;!$(N[C,S]=[O,S,N]):12].CC(C)(C)OC(=O)[N:13]-[C:14]",
                subsets=PHASE_3_N_BOC_NH_SUBSETS,
            ),  # NH_N-Boc
            QueryMol("[#6:21]-[#7]=C=O", subsets=YV_SUBSETS),  # YV
            QueryMol("[#6:31]-[#6;D3:32]=O", subsets=QS_SUBSETS),  # QS
        ],
        product=QueryMol(
            "[#6:21]-[#7]-[#6](=O)-[#7:12]-[#6:11].[#6:31]-[#6:32]-[#7:13]-[#6:14]"
        ),
        reaction_id=22,
        chemical_space="wuxi",
    ),
    # NH_N-Boc + LP + SA/SB/SC/SD/SE
    Reaction(
        reactants=[
            QueryMol(
                "[#6:11]-[N&!H0;!$(N[C,S]=[O,S,N]):12].CC(C)(C)OC(=O)[N:13]",
                subsets=PHASE_3_N_BOC_NH_SUBSETS,
            ),  # NH_N-Boc
            QueryMol("[#6:21]-[#17]", subsets=LP_SUBSETS),  # LP
            QueryMol(
                "[O:32]=[C:31]([OH,O-])[*:33]", subsets=S_SUBSETS
            ),  # SA/SB/SC/SD/SE
        ],
        product=QueryMol("[#6:11]-[#7:12]-[#6:21].[O:32]=[C:31]([N:13])[*:33]"),
        reaction_id=23,
        chemical_space="wuxi",
    ),
    # NH_N-Boc + LP + YV
    Reaction(
        reactants=[
            QueryMol(
                "[#6:11]-[N&!H0;!$(N[C,S]=[O,S,N]):12].CC(C)(C)OC(=O)[N:13]-[C:14]",
                subsets=PHASE_3_N_BOC_NH_SUBSETS,
            ),  # NH_N-Boc
            QueryMol("[#6:21]-[#17]", subsets=LP_SUBSETS),  # LP
            QueryMol("[#6:31]-[#7]=C=O", subsets=YV_SUBSETS),  # YV
        ],
        product=QueryMol(
            "[#6:11]-[#7:12]-[#6:21].[#6:31]-[#7]-[#6](=O)-[#7:13]-[#6:14]"
        ),
        reaction_id=24,
        chemical_space="wuxi",
    ),
    # NH_N-Boc + LP + LP
    Reaction(
        reactants=[
            QueryMol(
                "[#6:11]-[N&!H0;!$(N[C,S]=[O,S,N]):12].CC(C)(C)OC(=O)[N:13]-[C:14]",
                subsets=PHASE_3_N_BOC_NH_SUBSETS,
            ),  # NH_N-Boc
            QueryMol("[#6:21]-[#17]", subsets=LP_SUBSETS),  # LP
            QueryMol("[#6:31]-[#17]", subsets=LP_SUBSETS),  # LP
        ],
        product=QueryMol("[#6:11]-[#7:12]-[#6:21].[#6:14]-[#7:13]-[#6:31]"),
        reaction_id=25,
        chemical_space="wuxi",
    ),
    # NH_N-Boc + LP + HU
    Reaction(
        reactants=[
            QueryMol(
                "[#6:11]-[N&!H0;!$(N[C,S]=[O,S,N]):12].CC(C)(C)OC(=O)[N:13]",
                subsets=PHASE_3_N_BOC_NH_SUBSETS,
            ),  # NH_N-Boc
            QueryMol("[#6:21]-[#17]", subsets=LP_SUBSETS),  # LP
            QueryMol("Cl[S:31]([*:33])(=[O:32])=[O:34]", subsets=HU_SUBSETS),  # HU
        ],
        product=QueryMol(
            "[#6:11]-[#7:12]-[#6:21].[N:13][S:31]([*:33])(=[O:32])=[O:34]"
        ),
        reaction_id=26,
        chemical_space="wuxi",
    ),
    # NH_N-Boc + LP + QS
    Reaction(
        reactants=[
            QueryMol(
                "[#6:11]-[N&!H0;!$(N[C,S]=[O,S,N]):12].CC(C)(C)OC(=O)[N:13]-[C:14]",
                subsets=PHASE_3_N_BOC_NH_SUBSETS,
            ),  # NH_N-Boc
            QueryMol("[#6:21]-[#17]", subsets=LP_SUBSETS),  # LP
            QueryMol("[#6:31]-[#6;D3:32]=O", subsets=QS_SUBSETS),  # QS
        ],
        product=QueryMol("[#6:11]-[#7:12]-[#6:21].[#6:31]-[#6:32]-[#7:13]-[#6:14]"),
        reaction_id=27,
        chemical_space="wuxi",
    ),
    # NH_N-Boc + HU + SA/SB/SC/SD/SE
    Reaction(
        reactants=[
            QueryMol(
                "[N&!H0;!$(N[C,S]=[O,S,N]):11].CC(C)(C)OC(=O)[N:12]",
                subsets=PHASE_3_N_BOC_NH_SUBSETS,
            ),  # NH_N-Boc
            QueryMol("Cl[S:21]([*:23])(=[O:22])=[O:24]", subsets=HU_SUBSETS),  # HU
            QueryMol(
                "[O:32]=[C:31]([OH,O-])[*:33]", subsets=S_SUBSETS
            ),  # SA/SB/SC/SD/SE
        ],
        product=QueryMol(
            "[N:11][S:21]([*:23])(=[O:22])=[O:24].[O:32]=[C:31]([N:12])[*:33]"
        ),
        reaction_id=28,
        chemical_space="wuxi",
    ),
    # NH_N-Boc + HU + YV
    Reaction(
        reactants=[
            QueryMol(
                "[N&!H0;!$(N[C,S]=[O,S,N]):11].CC(C)(C)OC(=O)[N:12]-[C:13]",
                subsets=PHASE_3_N_BOC_NH_SUBSETS,
            ),  # NH_N-Boc
            QueryMol("Cl[S:21]([*:23])(=[O:22])=[O:24]", subsets=HU_SUBSETS),  # HU
            QueryMol("[#6:31]-[#7]=C=O", subsets=YV_SUBSETS),  # YV
        ],
        product=QueryMol(
            "[N:11][S:21]([*:23])(=[O:22])=[O:24].[#6:31]-[#7]-[#6](=O)-[#7:12]-[#6:13]"
        ),
        reaction_id=29,
        chemical_space="wuxi",
    ),
    # NH_N-Boc + HU + LP
    Reaction(
        reactants=[
            QueryMol(
                "[N&!H0;!$(N[C,S]=[O,S,N]):11].CC(C)(C)OC(=O)[N:12]-[C:13]",
                subsets=PHASE_3_N_BOC_NH_SUBSETS,
            ),  # NH_N-Boc
            QueryMol("Cl[S:21]([*:23])(=[O:22])=[O:24]", subsets=HU_SUBSETS),  # HU
            QueryMol("[#6:31]-[#17]", subsets=LP_SUBSETS),  # LP
        ],
        product=QueryMol(
            "[N:11][S:21]([*:23])(=[O:22])=[O:24].[#6:13]-[#7:12]-[#6:31]"
        ),
        reaction_id=30,
        chemical_space="wuxi",
    ),
    # NH_N-Boc + HU + HU
    Reaction(
        reactants=[
            QueryMol(
                "[N&!H0;!$(N[C,S]=[O,S,N]):11].CC(C)(C)OC(=O)[N:12]",
                subsets=PHASE_3_N_BOC_NH_SUBSETS,
            ),  # NH_N-Boc
            QueryMol("Cl[S:21]([*:23])(=[O:22])=[O:24]", subsets=HU_SUBSETS),  # HU
            QueryMol("Cl[S:31]([*:33])(=[O:32])=[O:34]", subsets=HU_SUBSETS),  # HU
        ],
        product=QueryMol(
            "[N:11][S:21]([*:23])(=[O:22])=[O:24].[N:12][S:31]([*:33])(=[O:32])=[O:34]"
        ),
        reaction_id=31,
        chemical_space="wuxi",
    ),
    # NH_N-Boc + HU + QS
    Reaction(
        reactants=[
            QueryMol(
                "[N&!H0;!$(N[C,S]=[O,S,N]):11].CC(C)(C)OC(=O)[N:12]-[C:13]",
                subsets=PHASE_3_N_BOC_NH_SUBSETS,
            ),  # NH_N-Boc
            QueryMol("Cl[S:21]([*:23])(=[O:22])=[O:24]", subsets=HU_SUBSETS),  # HU
            QueryMol("[#6:31]-[#6;D3:32]=O", subsets=QS_SUBSETS),  # QS
        ],
        product=QueryMol(
            "[N:11][S:21]([*:23])(=[O:22])=[O:24].[#6:31]-[#6:32]-[#7:12]-[#6:13]"
        ),
        reaction_id=32,
        chemical_space="wuxi",
    ),
    # NH_N-Boc + SA/SB/SC/SD/SE + YV
    Reaction(
        reactants=[
            QueryMol(
                "[N&!H0;!$(N[C,S]=[O,S,N]):11].CC(C)(C)OC(=O)[N:12]-[C:13]",
                subsets=PHASE_3_N_BOC_NH_SUBSETS,
            ),  # NH_N-Boc
            QueryMol(
                "[O:22]=[C:21]([OH,O-])[*:23]", subsets=S_SUBSETS
            ),  # SA/SB/SC/SD/SE
            QueryMol("[#6:31]-[#7]=C=O", subsets=YV_SUBSETS),  # YV
        ],
        product=QueryMol(
            "[O:22]=[C:21]([N:11])[*:23].[#6:31]-[#7]-[#6](=O)-[#7:12]-[#6:13]"
        ),
        reaction_id=33,
        chemical_space="wuxi",
    ),
    # NH_N-Boc + SA/SB/SC/SD/SE + LP
    Reaction(
        reactants=[
            QueryMol(
                "[N&!H0;!$(N[C,S]=[O,S,N]):11].CC(C)(C)OC(=O)[N:12]-[C:13]",
                subsets=PHASE_3_N_BOC_NH_SUBSETS,
            ),  # NH_N-Boc
            QueryMol(
                "[O:22]=[C:21]([OH,O-])[*:23]", subsets=S_SUBSETS
            ),  # SA/SB/SC/SD/SE
            QueryMol("[#6:31]-[#17]", subsets=LP_SUBSETS),  # LP
        ],
        product=QueryMol("[O:22]=[C:21]([N:11])[*:23].[#6:13]-[#7:12]-[#6:31]"),
        reaction_id=34,
        chemical_space="wuxi",
    ),
    # NH_N-Boc + SA/SB/SC/SD/SE + HU
    Reaction(
        reactants=[
            QueryMol(
                "[N&!H0;!$(N[C,S]=[O,S,N]):11].CC(C)(C)OC(=O)[N:12]",
                subsets=PHASE_3_N_BOC_NH_SUBSETS,
            ),  # NH_N-Boc
            QueryMol(
                "[O:22]=[C:21]([OH,O-])[*:23]", subsets=S_SUBSETS
            ),  # SA/SB/SC/SD/SE
            QueryMol("Cl[S:31]([*:33])(=[O:32])=[O:34]", subsets=HU_SUBSETS),  # HU
        ],
        product=QueryMol(
            "[O:22]=[C:21]([N:11])[*:23].[N:12][S:31]([*:33])(=[O:32])=[O:34]"
        ),
        reaction_id=35,
        chemical_space="wuxi",
    ),
    # NH_N-Boc + SA/SB/SC/SD/SE + QS
    Reaction(
        reactants=[
            QueryMol(
                "[N&!H0;!$(N[C,S]=[O,S,N]):11].CC(C)(C)OC(=O)[N:12]-[C:13]",
                subsets=PHASE_3_N_BOC_NH_SUBSETS,
            ),  # NH_N-Boc
            QueryMol(
                "[O:22]=[C:21]([OH,O-])[*:23]", subsets=S_SUBSETS
            ),  # SA/SB/SC/SD/SE
            QueryMol("[#6:31]-[#6;D3:32]=O", subsets=QS_SUBSETS),  # QS
        ],
        product=QueryMol("[O:22]=[C:21]([N:11])[*:23].[#6:31]-[#6:32]-[#7:12]-[#6:13]"),
        reaction_id=36,
        chemical_space="wuxi",
    ),
)

WUXI_REACTIONS = (
    WUXI_REACTIONS_PHASE_1 + WUXI_REACTIONS_PHASE_2 + WUXI_REACTIONS_PHASE_3
)
