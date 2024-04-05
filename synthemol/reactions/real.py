"""SMARTS representations of some of the Enamine REAL reactions."""
from synthemol.reactions.query_mol import QueryMol
from synthemol.reactions.reaction import Reaction

HALIDE = "F,Cl,Br,I"
BOC = QueryMol("[H]C([H])([H])C([O][C](=[O])[*:1])(C([H])([H])[H])C([H])([H])[H]")
CO2TBU = QueryMol(
    "[H]C([H])([H])C([O:2][C:3](=[O:4])[*:1])(C([H])([H])[H])C([H])([H])[H]"
)
CO2ME = QueryMol("[H]C([H])([H])[O:2][C:3]([*:1])=[O:4]")
CO2ET = QueryMol("[H]C([H])([H])C([H])([H])[O:2][C:3]([*:1])=[O:4]")
CO2R = QueryMol("*[O:2][C:3]([*:1])=[O:4]")
POST_BOC_CLEAVAGE = QueryMol("[*:1]")
POST_ESTER_HYDROLYSIS = QueryMol("[H][O:2][C:3]([*:1])=[O:4]")

# Define post-reactions
BOC_CLEAVAGE = Reaction(reactants=[BOC], product=POST_BOC_CLEAVAGE,)
ESTER_HYDROLYSIS_CO2ME = Reaction(reactants=[CO2ME], product=POST_ESTER_HYDROLYSIS,)
ESTER_HYDROLYSIS_CO2ET = Reaction(reactants=[CO2ET], product=POST_ESTER_HYDROLYSIS,)
ESTER_HYDROLYSIS_CO2TBU = Reaction(reactants=[CO2TBU], product=POST_ESTER_HYDROLYSIS,)
ESTER_HYDROLYSIS_CO2R = Reaction(reactants=[CO2R], product=POST_ESTER_HYDROLYSIS,)


# Define REAL reactions
REAL_REACTIONS = (
    Reaction(
        reactants=[
            QueryMol(
                "[H]C([H])([H])C(OC(=O)[N:1]([*:2])[*:3])(C([H])([H])[H])C([H])([H])[H].[*:4][N:5]([H])[*:6]"
            ),
            QueryMol("[OH1][C:7]([*:8])=[O:9]"),
            QueryMol("[OH1][C:10]([*:11])=[O:12]"),
        ],
        product=QueryMol(
            "[*:4][N:5]([*:6])[C:7](=[O:9])[*:8].[*:3][N:1]([*:2])[C:10](=[O:12])[*:11]"
        ),
        chemical_space="real",
        reaction_id=275592,
    ),
    Reaction(
        reactants=[
            QueryMol("[*:1][N:2]([H])[*:3]"),
            QueryMol("[OH1][C:4]([*:5])=[O:6]"),
        ],
        product=QueryMol("[*:5][C:4](=[O:6])[N:2]([*:1])[*:3]"),
        chemical_space="real",
        reaction_id=22,
    ),
    Reaction(
        reactants=[
            QueryMol("[*:1][N:2]([H])[*:3]"),
            QueryMol("[OH1][C:4]([*:5])=[O:6]"),
        ],
        product=QueryMol("[*:5][C:4](=[O:6])[N:2]([*:1])[*:3]"),
        chemical_space="real",
        reaction_id=11,
    ),
    Reaction(
        reactants=[
            QueryMol("[*:1][N:2]([H])[*:3]"),
            QueryMol("[OH1][C:4]([*:5])=[O:6]"),
        ],
        product=QueryMol("[*:5][C:4](=[O:6])[N:2]([*:1])[*:3]"),
        chemical_space="real",
        reaction_id=527,
    ),
    Reaction(
        reactants=[
            QueryMol("[*:1][N:2]([H])[*:3]"),
            QueryMol("[OH1][C:4]([*:5])=[O:6]"),
        ],
        product=QueryMol("[*:5][C:4](=[O:6])[N:2]([*:1])[*:3]"),
        chemical_space="real",
        reaction_id=240690,
        post_reaction=BOC_CLEAVAGE,
    ),
    Reaction(
        reactants=[QueryMol("[*:1][N:2]([H])[H:3]"), QueryMol("[*:4][N:5]([H])[*:6]")],
        product=QueryMol("O=C([N:2]([*:1])[H:3])[N:5]([*:4])[*:6]"),
        chemical_space="real",
        reaction_id=2430,
    ),
    Reaction(
        reactants=[QueryMol("[*:1][N:2]([H])[H:3]"), QueryMol("[*:4][N:5]([H])[H:6]")],
        product=QueryMol("O=C([N:2]([*:1])[H:3])[N:5]([*:4])[H:6]"),
        chemical_space="real",
        reaction_id=2708,
    ),
    Reaction(
        reactants=[QueryMol("[*:1][N:2]([H])[*:3]"), QueryMol(f"[{HALIDE}][*:4]")],
        product=QueryMol("[*:1][N:2]([*:3])[*:4]"),
        chemical_space="real",
        reaction_id=2230,
    ),
    Reaction(
        reactants=[QueryMol("[*:1][N:2]([H])[H:3]"), QueryMol("[*:4][N:5]([H])[H:6]")],
        product=QueryMol("O=C(C(=O)[N:2]([*:1])[H:3])[N:5]([*:4])[H:6]"),
        chemical_space="real",
        reaction_id=2718,
    ),
    Reaction(
        reactants=[
            QueryMol("[*:1][N:2]([H])[*:3]"),
            QueryMol(f"[O:4]=[S:5](=[O:6])([{HALIDE}])[*:7]"),
        ],
        product=QueryMol("[O:4]=[S:5](=[O:6])([*:7])[N:2]([*:1])[*:3]"),
        chemical_space="real",
        reaction_id=40,
    ),
    Reaction(
        reactants=[QueryMol("[*:1][N:2]([H])[*:3]"), QueryMol(f"[{HALIDE}][*:4]")],
        product=QueryMol("[*:1][N:2]([*:3])[*:4]"),
        chemical_space="real",
        reaction_id=27,
    ),
    Reaction(
        reactants=[QueryMol("[*:4][N:5]([H])[H:6]"), QueryMol("[*:1][N:2]([H])[*:3]")],
        product=QueryMol("O=C(C(=O)[N:2]([*:1])[*:3])[N:5]([*:4])[H:6]"),
        chemical_space="real",
        reaction_id=271948,
        post_reaction=BOC_CLEAVAGE,
    ),
    Reaction(
        reactants=[QueryMol("[OH1:1][C:2]([*:3])=[O:4]"), QueryMol(f"[{HALIDE}][*:5]")],
        product=QueryMol("[O:4]=[C:2]([*:3])[O:1][*:5]"),
        chemical_space="real",
        reaction_id=1458,
    ),
    Reaction(
        reactants=[QueryMol("[*:1][NH,S,O:2]([H])"), QueryMol(f"[{HALIDE}][*:3]")],
        product=QueryMol("[*:1][NH,S,O:2]([*:3])"),
        chemical_space="real",
        reaction_id=7,
    ),
    Reaction(
        reactants=[
            QueryMol("[H]O[C:3]([C:4]([H])([H])[N:5]([H])[C:6]([*:7])=[O:8])=[O:9]"),
            QueryMol("[*:1][C:2]([H])=O"),
        ],
        product=QueryMol("[H]/[C:2](=[C:4]1/[N:5]=[C:6]([*:7])[O:8][C:3]1=[O:9])[*:1]"),
        chemical_space="real",
        reaction_id=10,
    ),
    Reaction(
        reactants=[
            QueryMol("[NH2:5][*:6]"),
            QueryMol("[H]C([H])([H])C([H])([H])O[CH:1]=[C:2]([*:3])[*:4]"),
        ],
        product=QueryMol("[*:6][NH:5][CH:1]=[C:2]([*:3])[*:4]"),
        chemical_space="real",
        reaction_id=17,
    ),
    Reaction(
        reactants=[QueryMol("[*:1][NH,S,O:2]([H])"), QueryMol(f"[{HALIDE}][*:3]")],
        product=QueryMol("[*:1][NH,S,O:2]([*:3])"),
        chemical_space="real",
        reaction_id=34,
    ),
    Reaction(
        reactants=[QueryMol("[*:1][N:2]([H])[*:3]"), QueryMol(f"[{HALIDE}][*:4]")],
        product=QueryMol("[*:1][N:2]([*:3])[*:4]"),
        chemical_space="real",
        reaction_id=38,
    ),
    Reaction(
        reactants=[
            QueryMol("[H][N:1]([*:2])[C:3](=[S:4])[N:5]([H])[N:6]([H])[H]"),
            QueryMol(f"[H][C:7]([{HALIDE}])([*:8])[C:9](=O)[*:10]"),
        ],
        product=QueryMol(
            "[H][N:1]([*:2])[C:3]1=[N:5][N:6]=[C:9]([*:10])[C:7]([H])([*:8])[S:4]1"
        ),
        chemical_space="real",
        reaction_id=41,
    ),
    Reaction(
        reactants=[
            QueryMol("[H][N:5]([H])[c:6][c:7][C:8]([#6:9])=O"),
            QueryMol("[H][C:1]([H])([*:2])[C:3](=O)[*:4]"),
        ],
        product=QueryMol("[#6:9][c:8]1[c:7][c:6][n:5][c:3]([*:4])[c:1]1[*:2]"),
        chemical_space="real",
        reaction_id=42,
    ),
    Reaction(
        reactants=[
            QueryMol(f"[H][C:1]([{HALIDE}])([*:2])[C:3](=O)[*:4]")
        ],  # TODO: how to make this appear as reactant 2 instead of 1?
        product=QueryMol("[H]N([H])c1n[c:3]([*:4])[c:1]([*:2])s1"),
        chemical_space="real",
        reaction_id=43,
    ),
    Reaction(
        reactants=[QueryMol("[*:1][N:2]([H])[*:3]"), QueryMol(f"[{HALIDE}][*:4]")],
        product=QueryMol("[*:1][N:2]([*:3])[*:4]"),
        chemical_space="real",
        reaction_id=44,
    ),
    Reaction(
        reactants=[
            QueryMol("[H][N:1]([H])[N:2]([H])[*:3]"),
            QueryMol("[H][C:4]([H])([C:5](=O)[*:6])[C:7](=O)[*:8]"),
        ],
        product=QueryMol(
            "[H][c:4]1[c:7]([*:8]([H])([H])[H])[n:1][n:2]([*:3])[c:5]1[*:6]"
        ),
        chemical_space="real",
        reaction_id=49,
        sub_reaction_id=1,
    ),
    Reaction(
        reactants=[
            QueryMol("[H][N:1]([H])[N:2]([H])[*:3]"),
            QueryMol("[H][C:4]([H])([C:5](=O)[*:6])[C:7]([H])([H])[C:8](=O)[*:9]"),
        ],
        product=QueryMol(
            "[H][C:7]1([H])[C:8]([*:9])=[N:1][N:2]([*:3])[C:5]([H])([*:6])[C:4]1([H])[H]"
        ),
        chemical_space="real",
        reaction_id=49,
        sub_reaction_id=2,
    ),
    Reaction(
        reactants=[QueryMol("[H][N:1]([H])[N:2]([H])[*:3]")],
        product=QueryMol(
            "[H]Oc1c([H])c([H])c([H])c([H])c1C(=O)c1c([H])[n:1][n:2]([*:3])c1[H]"
        ),
        chemical_space="real",
        reaction_id=52,
    ),
    Reaction(
        reactants=[
            QueryMol("[H][O:1][c:2][c:3][C:4]([H])=O"),
            QueryMol(f"[H][C:5]([H])([*:6])[{HALIDE}]"),
        ],
        product=QueryMol("[H][C:4]1=[C:5]([*:6])[O:1][c:2][c:3]1"),
        chemical_space="real",
        reaction_id=55,
    ),
    Reaction(
        reactants=[
            QueryMol("[H][N:3]([H])[c:4][c:5][N:6]"),
            QueryMol("[H][C:1](=O)[*:2]"),
        ],
        product=QueryMol("[H][N:3]1[C:1]([*:2])=[N:6][c:5][c:4]1"),
        chemical_space="real",
        reaction_id=60,
    ),
    Reaction(
        reactants=[QueryMol("[*:1][N:2]([H])[*:3]"), QueryMol(f"[{HALIDE}][*:4]")],
        product=QueryMol("[*:1][N:2]([*:3])[*:4]"),
        chemical_space="real",
        reaction_id=61,
    ),
    Reaction(
        reactants=[QueryMol(f"[{HALIDE}][*:1]"), QueryMol(f"[{HALIDE}][*:2]")],
        product=QueryMol("[*:1]S[*:2]"),
        chemical_space="real",
        reaction_id=62,
    ),
    Reaction(
        reactants=[QueryMol("[*:1][Cl]"), QueryMol("[*:2][Cl]")],
        product=QueryMol("O=S(=O)([*:1])[*:2]"),
        chemical_space="real",
        reaction_id=1500,
        sub_reaction_id=1,
    ),
    Reaction(
        reactants=[QueryMol("[*:1][Cl]"), QueryMol("[H]OS(=O)(=O)[*:2]")],
        product=QueryMol("O=S(=O)([*:1])[*:2]"),
        chemical_space="real",
        reaction_id=1500,
        sub_reaction_id=2,
    ),
    Reaction(
        reactants=[QueryMol("[H]OS(=O)(=O)[*:1]"), QueryMol("[*:2][Cl]")],
        product=QueryMol("O=S(=O)([*:1])[*:2]"),
        chemical_space="real",
        reaction_id=1500,
        sub_reaction_id=3,
    ),
    Reaction(
        reactants=[QueryMol("[H]OS(=O)(=O)[*:1]"), QueryMol("[H]OS(=O)(=O)[*:2]")],
        product=QueryMol("O=S(=O)([*:1])[*:2]"),
        chemical_space="real",
        reaction_id=1500,
        sub_reaction_id=4,
    ),
    Reaction(
        reactants=[
            QueryMol("[H][N:6]([H])[*:7]"),
            QueryMol("[N:1][c:2][c:3][C:4]#[N:5]"),
        ],
        product=QueryMol("[H]c1[n:1][c:2][c:3][c:4]([N:6]([H])[*:7])[n:5]1"),
        chemical_space="real",
        reaction_id=1982,
    ),
    Reaction(
        reactants=[QueryMol("[*:1][SH]"), QueryMol("[*:2][Cl]")],
        product=QueryMol("O=S(=O)([*:1])[*:2]"),
        chemical_space="real",
        reaction_id=2714,
        sub_reaction_id=1,
    ),
    Reaction(
        reactants=[QueryMol("[*:1][SH]"), QueryMol("[H]OS(=O)(=O)[*:2]")],
        product=QueryMol("O=S(=O)([*:1])[*:2]"),
        chemical_space="real",
        reaction_id=2714,
        sub_reaction_id=2,
    ),
    Reaction(
        reactants=[
            QueryMol(r"[H][O:1]/[N:2]=[C:3](\[*:4])[N:5]([H])[H]"),
            QueryMol("[H]O[C:6](=O)[*:7]"),
        ],
        product=QueryMol("[*:4][c:3]1[n:2][o:1][c:6]([*:7])[n:5]1"),
        chemical_space="real",
        reaction_id=265764,
    ),
    Reaction(
        reactants=[QueryMol("[*:1][N:2]([H])[*:3]"), QueryMol(f"[{HALIDE}][*:4]")],
        product=QueryMol("[*:1][N:2]([*:3])[*:4]"),
        chemical_space="real",
        reaction_id=269956,
        post_reaction=ESTER_HYDROLYSIS_CO2TBU,
    ),
    Reaction(
        reactants=[QueryMol("[*:1][N:2]([H])[*:3]"), QueryMol(f"[{HALIDE}][*:4]")],
        product=QueryMol("[*:1][N:2]([*:3])[*:4]"),
        chemical_space="real",
        reaction_id=269982,
        post_reaction=BOC_CLEAVAGE,
    ),
    Reaction(
        reactants=[
            QueryMol(r"[H][O:1]/[N:2]=[C:3](\[*:4])[N:5]([H])[H]"),
            QueryMol("[H]O[C:6](=O)[*:7]"),
        ],
        product=QueryMol("[*:4][c:3]1[n:2][o:1][c:6]([*:7])[n:5]1"),
        chemical_space="real",
        reaction_id=270036,
    ),
    Reaction(
        reactants=[
            QueryMol("[*:1][N:2]([H])[*:3]"),
            QueryMol("O=C(O[C:4](=[O:5])[*:6])[*]"),
        ],
        product=QueryMol("[O:5]=[C:4]([*:6])[N:2]([*:1])[*:3]"),
        chemical_space="real",
        reaction_id=270062,
    ),
    Reaction(
        reactants=[QueryMol("[H][N:1]([H])[*:2]"), QueryMol("[H]O[C:3](=[O:4])[*:5]"),],
        product=QueryMol("[H][N:1]([*:2])[C:3](=[O:4])[*:5]"),
        chemical_space="real",
        reaction_id=270112,
        post_reaction=ESTER_HYDROLYSIS_CO2ME,
    ),
    Reaction(
        reactants=[QueryMol("[*:1][N:2]([H])[*:3]"), QueryMol(f"[{HALIDE}][*:4]")],
        product=QueryMol("[*:1][N:2]([*:3])[*:4]"),
        chemical_space="real",
        reaction_id=270122,
        post_reaction=ESTER_HYDROLYSIS_CO2R,
    ),
    Reaction(
        reactants=[QueryMol("[*:1][N:2]([H])[*:3]"), QueryMol(f"[{HALIDE}][*:4]")],
        product=QueryMol("[*:1][N:2]([*:3])[*:4]"),
        chemical_space="real",
        reaction_id=270166,
        post_reaction=BOC_CLEAVAGE,
    ),
    Reaction(
        reactants=[QueryMol("[N:1]#[C:2][*:3]"), QueryMol("[H]O[C:4](=O)[*:5]")],
        product=QueryMol("[*:3][c:2]1no[c:4]([*:5])[n:1]1"),
        chemical_space="real",
        reaction_id=270196,
    ),
    Reaction(
        reactants=[QueryMol("[*:1][N:2]([H])[*:3]"), QueryMol(f"[{HALIDE}][*:4]")],
        product=QueryMol("[*:1][N:2]([*:3])[*:4]"),
        chemical_space="real",
        reaction_id=270344,
        sub_reaction_id=1,
        post_reaction=ESTER_HYDROLYSIS_CO2ET,
    ),
    Reaction(
        reactants=[QueryMol("[*:1][SH:2]"), QueryMol(f"[{HALIDE}][*:3]")],
        product=QueryMol("[*:1]S[*:3]"),
        chemical_space="real",
        reaction_id=270344,
        sub_reaction_id=2,
        post_reaction=ESTER_HYDROLYSIS_CO2ET,
    ),
    Reaction(
        reactants=[QueryMol(r"[H][O:1]/[N:2]=[C:3](\[*:4])[N:5]([H])[H]")],
        product=QueryMol("[H][n:5]1[c:3]([*:4])[n:2][o:1]c1=O"),
        chemical_space="real",
        reaction_id=270690,
    ),
    Reaction(
        reactants=[
            QueryMol("S=[C:1]=[N:2][*:3]"),
            QueryMol("[H][N:4]([*:5])[*:6]"),
            QueryMol("[H][N:7]([H])[N:8]([H])[C:9](=O)[*:10]"),
        ],
        product=QueryMol("[*:5][N:4]([*:6])[c:1]1[n:7][n:8][c:9]([*:10])[n:2]1[*:3]"),
        chemical_space="real",
        reaction_id=270942,
    ),
    Reaction(
        reactants=[QueryMol("[*:4][N:5]([H])[*:6]"), QueryMol("[O:1]=[C:2]([*:3])O*"),],
        product=QueryMol("[*:3][C:2]([N:5]([*:4])[*:6])=[O:1]"),
        chemical_space="real",
        reaction_id=271144,
        post_reaction=ESTER_HYDROLYSIS_CO2TBU,
    ),
    Reaction(
        reactants=[QueryMol("[N:1]#[C:2][*:3]"), QueryMol("[H]O[C:4](=O)[*:5]")],
        product=QueryMol("[*:3][c:2]1no[c:4]([*:5])[n:3]1"),
        chemical_space="real",
        reaction_id=271362,
        post_reaction=BOC_CLEAVAGE,
    ),
    Reaction(
        reactants=[QueryMol("[H]OB([*:1])O[H]"), QueryMol(f"[{HALIDE}][*:2]")],
        product=QueryMol("[*:1][*:2]"),
        chemical_space="real",
        reaction_id=271570,
    ),
    Reaction(
        reactants=[
            QueryMol("[H][N:1]([H])[N:2]([H])[C:3](=[O:4])[*:5]"),
            QueryMol("[H]O[C:6](=O)[*:7]"),
        ],
        product=QueryMol("[*:5][c:3]1[n:2][n:1][c:6]([*:7])[o:4]1"),
        chemical_space="real",
        reaction_id=271722,
    ),
    Reaction(
        reactants=[
            QueryMol(f"[H][C:1]([{HALIDE}])([*:2])[C:3](=O)[*:4]"),
        ],  # TODO: how to make this appear as reactant 2 instead of 1?
        product=QueryMol("[H]N([H])c1n[c:3]([*:4])[c:1]([*:2])s1"),
        chemical_space="real",
        reaction_id=272150,
        post_reaction=BOC_CLEAVAGE,
    ),
    Reaction(
        reactants=[QueryMol("[H][N:1]([H])[*:2]"), QueryMol("[H]O[C:3](=[O:4])[*:5]"),],
        product=QueryMol("[H][N:1]([*:2])[C:3](=[O:4])[*:5]"),
        chemical_space="real",
        reaction_id=272270,
        post_reaction=ESTER_HYDROLYSIS_CO2TBU,
    ),
    Reaction(
        reactants=[QueryMol("[H][N:1]([H])[*:2]"), QueryMol("[H]O[C:3](=[O:4])[*:5]")],
        product=QueryMol("[H][N:1]([*:2])[C:3](=[O:4])[*:5]"),
        chemical_space="real",
        reaction_id=272430,
        post_reaction=BOC_CLEAVAGE,
    ),
    Reaction(
        reactants=[QueryMol("[*:1][N:2]([H])[*:3]"), QueryMol(f"[{HALIDE}][*:4]")],
        product=QueryMol("[*:1][N:2]([*:3])[*:4]"),
        chemical_space="real",
        reaction_id=272692,
        sub_reaction_id=1,
        post_reaction=BOC_CLEAVAGE,
    ),
    Reaction(
        reactants=[QueryMol("[*:1][OH:2]"), QueryMol(f"[{HALIDE}][*:3]")],
        product=QueryMol("[*:1]O[*:3]"),
        chemical_space="real",
        reaction_id=272692,
        sub_reaction_id=2,
        post_reaction=BOC_CLEAVAGE,
    ),
    Reaction(
        reactants=[QueryMol("[*:1][N:2]([H])[*:3]"), QueryMol(f"[{HALIDE}][*:4]")],
        product=QueryMol("[*:1][N:2]([*:3])[*:4]"),
        chemical_space="real",
        reaction_id=272710,
        sub_reaction_id=1,
        post_reaction=ESTER_HYDROLYSIS_CO2TBU,
    ),
    Reaction(
        reactants=[QueryMol("[*:1][OH,SH:2]"), QueryMol(f"[{HALIDE}][*:3]")],
        product=QueryMol("[*:1][O,S][*:3]"),
        chemical_space="real",
        reaction_id=272710,
        sub_reaction_id=2,
        post_reaction=ESTER_HYDROLYSIS_CO2TBU,
    ),
    Reaction(
        reactants=[QueryMol("[H][N:1]([H])[*:2]"), QueryMol("[H][N:3]([H])[*:4]")],
        product=QueryMol("[H][N:1]([*:2])C(=O)[N:3]([H])[*:4]"),
        chemical_space="real",
        reaction_id=273390,
        post_reaction=ESTER_HYDROLYSIS_CO2TBU,
    ),
    Reaction(
        reactants=[QueryMol("[H][N:1]([H])[*:2]"), QueryMol("[H][N:3]([*:4])[*:5]")],
        product=QueryMol("[H][N:1](C(=O)[N:3]([*:4])[*:5])[*:2]"),
        chemical_space="real",
        reaction_id=273392,
        post_reaction=ESTER_HYDROLYSIS_CO2TBU,
    ),
    Reaction(
        reactants=[QueryMol("[H][N:1]([H])[*:2]"), QueryMol("[H][N:3]([*:4])[*:5]")],
        product=QueryMol("[H][N:1](C(=O)[N:3]([*:4])[*:5])[*:2]"),
        chemical_space="real",
        reaction_id=273452,
        post_reaction=ESTER_HYDROLYSIS_CO2TBU,
    ),
    Reaction(
        reactants=[
            QueryMol("[H][N:1]([H])[*:2]"),
            QueryMol("[H]C([H])(OC(=O)[N:3]([*:4])[*:5])C(F)(F)F"),
        ],
        product=QueryMol("[H][N:1](C(=O)[N:3]([*:4])[*:5])[*:2]"),
        chemical_space="real",
        reaction_id=273454,
        post_reaction=BOC_CLEAVAGE,
    ),
    Reaction(
        reactants=[QueryMol("[H][N:4]([H])[*:5]"), QueryMol("[H][N:1]([*:2])[*:3]")],
        product=QueryMol("[H][N:4](C(=O)[N:1]([*:2])[*:3])[*:5]"),
        chemical_space="real",
        reaction_id=273458,
        post_reaction=BOC_CLEAVAGE,
    ),
    Reaction(
        reactants=[QueryMol("[H][N:1]([H])[*:2]"), QueryMol("[H][N:3]([*:4])[*:5]")],
        product=QueryMol("[H][N:1](C(=O)[N:3]([*:4])[*:5])[*:2]"),
        chemical_space="real",
        reaction_id=273464,
        post_reaction=BOC_CLEAVAGE,
    ),
    Reaction(
        reactants=[QueryMol("[H][N:1]([*:2])[*:3]"), QueryMol("[H][N:4]([H])[*:5]")],
        product=QueryMol("[H][N:4]([*:5])C(=O)C(=O)[N:1]([*:2])[*:3]"),
        chemical_space="real",
        reaction_id=273494,
        post_reaction=ESTER_HYDROLYSIS_CO2TBU,
    ),
    Reaction(
        reactants=[QueryMol("[H][N:4]([H])[*:5]"), QueryMol("[H][N:1]([*:2])[*:3]")],
        product=QueryMol("[H][N:4]([*:5])C(=O)C(=O)[N:1]([*:2])[*:3]"),
        chemical_space="real",
        reaction_id=273496,
        post_reaction=BOC_CLEAVAGE,
    ),
    Reaction(
        reactants=[
            QueryMol("[H][N:1]([H])[*:2]"),
            QueryMol("[H]C([H])(OC(=O)[N:3]([*:4])[*:5])C(F)(F)F"),
        ],
        product=QueryMol("[H][N:1](C(=O)[N:3]([*:4])[*:5])[*:2]"),
        chemical_space="real",
        reaction_id=273574,
        post_reaction=BOC_CLEAVAGE,
    ),
    Reaction(
        reactants=[QueryMol("[H]OB([*:1])O[H]"), QueryMol(f"[{HALIDE}][*:2]")],
        product=QueryMol("[*:1][*:2]"),
        chemical_space="real",
        reaction_id=273584,
        post_reaction=BOC_CLEAVAGE,
    ),
    Reaction(
        reactants=[QueryMol("[H][N:1]([H])[*:2]"), QueryMol("[H]O[C:3](=[O:4])[*:5]")],
        product=QueryMol("[H][N:1]([*:2])[C:3](=[O:4])[*:5]"),
        chemical_space="real",
        reaction_id=273652,
        post_reaction=ESTER_HYDROLYSIS_CO2TBU,
    ),
    Reaction(
        reactants=[QueryMol("[H]OB([*:1])O[H]"), QueryMol(f"[{HALIDE}][*:2]")],
        product=QueryMol("[*:1][*:2]"),
        chemical_space="real",
        reaction_id=273712,
    ),
    Reaction(
        reactants=[
            QueryMol("[N-:1]=[N+:2]=[N:3][*:4]"),
            QueryMol("[H][C:5]#[C:6][*:7]"),
        ],
        product=QueryMol("[H][c:5]1[c:6]([*:7])[n:1][n:2][n:3]1[*:4]"),
        chemical_space="real",
        reaction_id=273910,
    ),
    Reaction(
        reactants=[
            QueryMol(
                "[H]C([H])([H])O[C:1]1=[N:2][C:3]([H])([H])[O,CH2:4][O,CH2:5][C:6]1([H])[H]"
            ),
            QueryMol("[H][N:7]([H])[N:8]([H])[C:9](=O)[*:10]"),
        ],
        product=QueryMol(
            "[H][C:6]1([H])[O,CH2:5][O,CH2:4][C:3]([H])([H])[n:2]2[c:9]([*:10])[n:8][n:7][c:1]21"
        ),
        chemical_space="real",
        reaction_id=274052,
    ),
    Reaction(
        reactants=[
            QueryMol("[H][N:1]([H])[C:2](=[S:3])[*:4]"),
            QueryMol(
                f"[H][C:5]([*:6])([{HALIDE}])[C:7](=O)[C:8]#[C:9][Si](C([H])([H])[H])(C([H])([H])[H])C([H])([H])[H]"
            ),
            QueryMol("[N-:10]=[N+:11]=[N:12][*:13]"),
        ],
        product=QueryMol(
            "[H][c:9]1[c:8](-[c:7]2[n:1][c:2]([*:4])[s:3][c:5]2[*:6])[n:10][n:11][n:12]1[*:13]"
        ),
        chemical_space="real",
        reaction_id=274972,
    ),
    Reaction(
        reactants=[
            QueryMol("[N-:6]=[N+:7]=[N:8][*:9]"),
            QueryMol("[H]O[C:10](=[O:11])[*:12]"),
            QueryMol(
                "[H][C:1]#[C:2][*:3].[H]C([H])([H])C(OC(=O)[N:4]([H])[*:5])(C([H])([H])[H])C([H])([H])[H]"
            ),
        ],
        product=QueryMol(
            "[H][N:4]([*:5])[C:10](=[O:11])[*:12].[H][c:1]1[c:2]([*:3])[n:6][n:7][n:8]1[*:9]"
        ),
        chemical_space="real",
        reaction_id=276010,
        sub_reaction_id=1,
    ),
    Reaction(
        reactants=[
            QueryMol(
                "[N-:1]=[N+:2]=[N:3][*:4].[H]C([H])([H])C(OC(=O)[N:5]([H])[*:6])(C([H])([H])[H])C([H])([H])[H]"
            ),
            QueryMol("[H]O[C:10](=[O:11])[*:12]"),
            QueryMol("[H][C:7]#[C:8][*:9]"),
        ],
        product=QueryMol(
            "[H][N:5]([*:6])[C:10](=[O:11])[*:12].[H][c:7]1[c:8]([*:9])[n:1][n:2][n:3]1[*:4]"
        ),
        chemical_space="real",
        reaction_id=276010,
        sub_reaction_id=2,
    ),
    Reaction(
        reactants=[
            QueryMol("[H][N:1]([H])[C:2](=S)[N:3]([*:4])[*:5]"),
            QueryMol("[H][N:6]([H])[N:7]([H])[C:8](=O)[*:9]"),
        ],
        product=QueryMol("[H][n:1]1[c:8]([*:9])[n:7][n:6][c:2]1[N:3]([*:4])[*:5]"),
        chemical_space="real",
        reaction_id=276090,
    ),
    Reaction(
        reactants=[
            QueryMol("[H]N=[C:1]([*:2])[N:3]([H])[H]"),
            QueryMol("[H]O[C:4](=O)[*:5]"),
            QueryMol("[H][N:6]([H])[N:7]([H])[*:8]"),
        ],
        product=QueryMol("[*:8][n:7]1[n:6][c:1]([*:2])[n:3][c:4]1[*:5]"),
        chemical_space="real",
        reaction_id=276630,
    ),
    Reaction(
        reactants=[
            QueryMol("[H][N:1]([H])[N:2]([H])[C:3](=O)[*:4]"),
            QueryMol("[H][N:5]=[C:6]([*:7])OC([H])([H])[H]"),
        ],
        product=QueryMol("[H][n:5]1[c:3]([*:4])[n:2][n:1][c:6]1[*:7]"),
        chemical_space="real",
        reaction_id=276670,
    ),
    Reaction(
        reactants=[
            QueryMol(
                "[H][c:1]1[c:2]([H])[c:3]([H])[c:4]([N:5]([H])[H])[c:6]([C:7](=O)[*:8])[c:9]1[H]"
            ),
            QueryMol(
                "[H]O[C:10](=[O:11])[C:12]([H])([*:13])[N:14]([H])C(=O)OC(C([H])([H])[H])(C([H])([H])[H])C([H])([H])[H]"
            ),
        ],
        product=QueryMol(
            "[H][c:1]1[c:2]([H])[c:3]([H])[c:4]2[c:6]([c:9]1[H])[C:7]([*:8])=[N:14][C:12]([H])([*:13])[C:10](=[O:11])[N:5]2[H]"
        ),
        chemical_space="real",
        reaction_id=276770,
    ),
    Reaction(
        reactants=[
            QueryMol("[H][N:1]([H])[N:2]([H])[*:3]"),
            QueryMol("[H][C:4](=[O:5])[c:6][c:7][C:8](=O)OC([H])([H])[H]"),
        ],
        product=QueryMol("[H][c:8]1[n:1][n:2]([*:3])[c:4](=[O:5])[c:6][c:7]1"),
        chemical_space="real",
        reaction_id=279312,  # TODO: check aromaticity
        post_reaction=BOC_CLEAVAGE,
    ),
    Reaction(
        reactants=[
            QueryMol(
                "[H]C([H])([H])O[C:1]1=[N:2][C:3]([H])([H])[O,CH2:4][O,CH2:5][C:6]1([H])[H]"
            ),
            QueryMol("[H][N:7]([H])[N:8]([H])[C:9](=O)[*:10]"),
        ],
        product=QueryMol(
            "[H][C:6]1([H])[O,CH2:5][O,CH2:4][C:3]([H])([H])[n:2]2[c:9]([*:10])[n:8][n:7][c:1]21"
        ),
        chemical_space="real",
        reaction_id=279370,
        post_reaction=BOC_CLEAVAGE,
    ),
)


# The 13 most common REAL reaction IDs (used in the SyntheMol MCTS and RL papers)
MOST_COMMON_REAL_IDS = {
    275592,
    22,
    11,
    527,
    240690,
    2430,
    2708,
    2230,
    2718,
    40,
    27,
    271948,
    1458,
}

# The 13 most common REAL reactions
MOST_COMMON_REAL_REACTIONS = tuple(
    reaction
    for reaction in REAL_REACTIONS
    if reaction.reaction_id in MOST_COMMON_REAL_IDS
)
