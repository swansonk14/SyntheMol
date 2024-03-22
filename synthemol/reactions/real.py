"""SMARTS representations of some of the Enamine REAL reactions."""
from synthemol.reactions.query_mol import QueryMol
from synthemol.reactions.reaction import Reaction

HALIDE = "F,Cl,Br,I"

REAL_REACTIONS = (
    Reaction(
        reactants=[
            QueryMol("CC(C)(C)OC(=O)[N:1]([*:2])[*:3].[*:4][N:5]([H])[*:6]"),
            QueryMol("[OH1][C:7]([*:8])=[O:9]"),
            QueryMol("[OH1][C:10]([*:11])=[O:12]"),
        ],
        product=QueryMol(
            "[*:4][N:5]([*:6])[C:7](=[O:9])[*:8].[*:3][N:1]([*:2])[C:10](=[O:12])[*:11]"
        ),
        reaction_id=275592,
        chemical_space="real",
    ),
    Reaction(
        reactants=[
            QueryMol("[*:1][N:2]([H])[*:3]"),
            QueryMol("[OH1][C:4]([*:5])=[O:6]"),
        ],
        product=QueryMol("[*:5][C:4](=[O:6])[N:2]([*:1])[*:3]"),
        reaction_id=22,
        chemical_space="real",
    ),
    Reaction(
        reactants=[
            QueryMol("[*:1][N:2]([H])[*:3]"),
            QueryMol("[OH1][C:4]([*:5])=[O:6]"),
        ],
        product=QueryMol("[*:5][C:4](=[O:6])[N:2]([*:1])[*:3]"),
        reaction_id=11,
        chemical_space="real",
    ),
    Reaction(
        reactants=[
            QueryMol("[*:1][N:2]([H])[*:3]"),
            QueryMol("[OH1][C:4]([*:5])=[O:6]"),
        ],
        product=QueryMol("[*:5][C:4](=[O:6])[N:2]([*:1])[*:3]"),
        reaction_id=527,
        chemical_space="real",
    ),
    Reaction(
        reactants=[
            QueryMol("[*:1][N:2]([H])[*:3]"),
            QueryMol("[OH1][C:4]([*:5])=[O:6]"),
        ],
        product=QueryMol("[*:5][C:4](=[O:6])[N:2]([*:1])[*:3]"),
        reaction_id=240690,
        chemical_space="real",
    ),  # TODO: Boc cleavage
    Reaction(
        reactants=[QueryMol("[*:1][N:2]([H])[H:3]"), QueryMol("[*:4][N:5]([H])[*:6]")],
        product=QueryMol("O=C([N:2]([*:1])[H:3])[N:5]([*:4])[*:6]"),
        reaction_id=2430,
        chemical_space="real",
    ),
    Reaction(
        reactants=[QueryMol("[*:1][N:2]([H])[H:3]"), QueryMol("[*:4][N:5]([H])[H:6]")],
        product=QueryMol("O=C([N:2]([*:1])[H:3])[N:5]([*:4])[H:6]"),
        reaction_id=2708,
        chemical_space="real",
    ),
    Reaction(
        reactants=[QueryMol("[*:1][N:2]([H])[*:3]"), QueryMol(f"[{HALIDE}][*:4]")],
        product=QueryMol("[*:1][N:2]([*:3])[*:4]"),
        reaction_id=2230,
        chemical_space="real",
    ),
    Reaction(
        reactants=[QueryMol("[*:1][N:2]([H])[H:3]"), QueryMol("[*:4][N:5]([H])[H:6]")],
        product=QueryMol("O=C(C(=O)[N:2]([*:1])[H:3])[N:5]([*:4])[H:6]"),
        reaction_id=2718,
        chemical_space="real",
    ),
    Reaction(
        reactants=[
            QueryMol("[*:1][N:2]([H])[*:3]"),
            QueryMol(f"[O:4]=[S:5](=[O:6])([{HALIDE}])[*:7]"),
        ],
        product=QueryMol("[O:4]=[S:5](=[O:6])([*:7])[N:2]([*:1])[*:3]"),
        reaction_id=40,
        chemical_space="real",
    ),
    Reaction(
        reactants=[QueryMol("[*:1][N:2]([H])[*:3]"), QueryMol(f"[{HALIDE}][*:4]")],
        product=QueryMol("[*:1][N:2]([*:3])[*:4]"),
        reaction_id=27,
        chemical_space="real",
    ),
    Reaction(
        reactants=[QueryMol("[*:1][N:2]([H])[*:3]"), QueryMol("[*:4][N:5]([H])[H:6]")],
        product=QueryMol("O=C(C(=O)[N:2]([*:1])[*:3])[N:5]([*:4])[H:6]"),
        reaction_id=271948,
        chemical_space="real",
    ),  # TODO: Boc cleavage
    Reaction(
        reactants=[QueryMol("[OH1:1][C:2]([*:3])=[O:4]"), QueryMol(f"[{HALIDE}][*:5]")],
        product=QueryMol("[O:4]=[C:2]([*:3])[O:1][*:5]"),
        reaction_id=1458,
        chemical_space="real",
    ),
    Reaction(
        reactants=[QueryMol("[*:1][NH,S,O:2]([H])"), QueryMol(f"[{HALIDE}][*:3]")],
        product=QueryMol("[*:1][NH,S,O:2]([*:3])"),
        reaction_id=7,
        chemical_space="real",
    ),
    Reaction(
        reactants=[
            QueryMol("[*:1][C:2]([H])=O"),
            QueryMol("[H]O[C:3]([C:4]([H])([H])[N:5]([H])[C:6]([*:7])=[O:8])=[O:9]"),
        ],
        product=QueryMol(
            "[H]/[C:2](=[C:4]1/[N:5]=[C:6]([*:7])[O:8][C:3]1=[O:9])[*:1]"
        ),  # TODO: should this be an isomeric SMILES? Should it be aromatic?
        reaction_id=10,
        chemical_space="real",
    ),
    Reaction(
        reactants=[
            QueryMol("[H]C([H])([H])C([H])([H])O[CH:1]=[C:2]([*:3])[*:4]"),
            QueryMol("[NH2:5][*:6]"),
        ],
        product=QueryMol("[*:6][NH:5][CH:1]=[C:2]([*:3])[*:4]"),
        reaction_id=17,
        chemical_space="real",
    ),
    Reaction(
        reactants=[QueryMol("[*:1][NH,S,O:2]([H])"), QueryMol(f"[{HALIDE}][*:3]")],
        product=QueryMol("[*:1][NH,S,O:2]([*:3])"),
        reaction_id=34,
        chemical_space="real",
    ),
    Reaction(
        reactants=[QueryMol("[*:1][N:2]([H])[*:3]"), QueryMol(f"[{HALIDE}][*:4]")],
        product=QueryMol("[*:1][N:2]([*:3])[*:4]"),
        reaction_id=38,
        chemical_space="real",
    ),
    Reaction(
        reactants=[
            QueryMol("[H][N:1]([*:2])[C:3](=[S:4])[N:5]([H])[N:6]([H])[H]"),
            QueryMol(f"[H][C:7]([{HALIDE}])([*:8])[C:9](=O)[*:10]"),
        ],
        product=QueryMol(
            "[H][N:1]([*:2])[C:3]1=[N:5][N:6]=[C:9]([*:10])[C:7]([H])([*:8])[S:4]1"
        ),
        reaction_id=41,
        chemical_space="real",
    ),
    Reaction(
        reactants=[
            QueryMol("[H][C:1]([H])([*:2])[C:3](=O)[*:4]"),
            QueryMol(""),
        ],  # TODO: figure out R group ring
        product=QueryMol(""),  # TODO: figure out R group ring
        reaction_id=42,
        chemical_space="real",
    ),
    Reaction(
        reactants=[QueryMol(f"[H][C:1]([{HALIDE}])([*:2])[C:3](=O)[*:4]")],
        product=QueryMol(
            "[H]N([H])c1n[c:3]([*:4])[c:1]([*:2])s1"
        ),  # TODO: check that this is aromatic
        reaction_id=43,
        chemical_space="real",
    ),
    Reaction(
        reactants=[QueryMol("[*:1][N:2]([H])[*:3]"), QueryMol(f"[{HALIDE}][*:4]")],
        product=QueryMol("[*:1][N:2]([*:3])[*:4]"),
        reaction_id=44,
        chemical_space="real",
    ),
    Reaction(
        reactants=[
            QueryMol("[H][N:1]([H])[N:2]([H])[*:3]"),
            QueryMol("[H][C:4]([H])([C:5](=O)[*:6])[C:7](=O)[*:8]"),
        ],
        product=QueryMol(
            "[H][c:4]1[c:7]([*:8]([H])([H])[H])[n:1][n:2]([*:3])[c:5]1[*:6]"
        ),  # TODO: check if aromatic
        reaction_id=49,
        chemical_space="real",
    ),
    Reaction(
        reactants=[
            QueryMol("[H][N:1]([H])[N:2]([H])[*:3]"),
            QueryMol("[H][C:4]([H])([C:5](=O)[*:6])[C:7]([H])([H])[C:8](=O)[*:9]"),
        ],
        product=QueryMol(
            "[H][C:7]1([H])[C:8]([*:9])=[N:1][N:2]([*:3])[C:5]([H])([*:6])[C:4]1([H])[H]"
        ),  # TODO: check if aromatic
        reaction_id=49,  # TODO: check if it's okay to reuse ID 49 for the two versions of it
        chemical_space="real",
    ),
    Reaction(
        reactants=[QueryMol("[H][N:1]([H])[N:2]([H])[*:3]")],
        product=QueryMol(
            "[H]Oc1c([H])c([H])c([H])c([H])c1C(=O)c1c([H])[n:1][n:2]([*:3])c1[H]"
        ),  # TODO: check if aromatic
        reaction_id=52,
        chemical_space="real",
    ),
    Reaction(
        reactants=[QueryMol(""), QueryMol("")],
        product=QueryMol(""),
        reaction_id=55,  # TODO: figure out R group ring
        chemical_space="real",
    ),
    Reaction(
        reactants=[QueryMol(""), QueryMol("")],
        product=QueryMol(""),
        reaction_id=60,  # TODO: figure out R group ring
        chemical_space="real",
    ),
    Reaction(
        reactants=[QueryMol("[*:1][N:2]([H])[*:3]"), QueryMol(f"[{HALIDE}][*:4]")],
        product=QueryMol("[*:1][N:2]([*:3])[*:4]"),
        reaction_id=61,
        chemical_space="real",
    ),
    Reaction(
        reactants=[QueryMol(f"[{HALIDE}][*:1]"), QueryMol(f"[{HALIDE}][*:2]")],
        product=QueryMol("[*:1]S[*:2]"),
        reaction_id=62,
        chemical_space="real",
    ),
    Reaction(
        reactants=[QueryMol("[*:1][Cl]"), QueryMol("[*:2][Cl]")],
        product=QueryMol("O=S(=O)([*:1])[*:2]"),
        reaction_id=1500,  # TODO: check if okay to reuse ID
        chemical_space="real",
    ),
    Reaction(
        reactants=[QueryMol("[*:1][Cl]"), QueryMol("[H]OS(=O)(=O)[*:2]")],
        product=QueryMol("O=S(=O)([*:1])[*:2]"),
        reaction_id=1500,  # TODO: check if okay to reuse ID
        chemical_space="real",
    ),
    Reaction(
        reactants=[QueryMol("[H]OS(=O)(=O)[*:1]"), QueryMol("[*:2][Cl]")],
        product=QueryMol("O=S(=O)([*:1])[*:2]"),
        reaction_id=1500,  # TODO: check if okay to reuse ID
        chemical_space="real",
    ),
    Reaction(
        reactants=[QueryMol("[H]OS(=O)(=O)[*:1]"), QueryMol("[H]OS(=O)(=O)[*:2]")],
        product=QueryMol("O=S(=O)([*:1])[*:2]"),
        reaction_id=1500,  # TODO: check if okay to reuse ID
        chemical_space="real",
    ),
    Reaction(
        reactants=[QueryMol("")],
        product=QueryMol(""),
        reaction_id=1982,  # TODO: figure out R group ring
        chemical_space="real",
    ),
    Reaction(
        reactants=[QueryMol("[*:1][SH]"), QueryMol("[*:2][Cl]")],
        product=QueryMol("O=S(=O)([*:1])[*:2]"),
        reaction_id=2714,  # TODO: check if okay to reuse ID
        chemical_space="real",
    ),
    Reaction(
        reactants=[QueryMol("[*:1][SH]"), QueryMol("[H]OS(=O)(=O)[*:2]")],
        product=QueryMol("O=S(=O)([*:1])[*:2]"),
        reaction_id=2714,  # TODO: check if okay to reuse ID
        chemical_space="real",
    ),
    Reaction(
        reactants=[
            QueryMol(
                r"[H][O:1]/[N:2]=[C:3](\[*:4])[N:5]([H])[H]"
            ),  # TODO: should this be isomeric?
            QueryMol("[H]O[C:6](=O)[*:7]"),
        ],
        product=QueryMol(
            "[*:4][c:3]1[n:2][o:1][c:6]([*:7])[n:5]1"
        ),  # TODO: check if aromatic
        reaction_id=265764,
        chemical_space="real",
    ),
    Reaction(
        reactants=[QueryMol("[*:1][N:2]([H])[*:3]"), QueryMol(f"[{HALIDE}][*:4]")],
        product=QueryMol("[*:1][N:2]([*:3])[*:4]"),
        reaction_id=269956,  # TODO: ester hydrolysis
        chemical_space="real",
    ),
    Reaction(
        reactants=[QueryMol("[*:1][N:2]([H])[*:3]"), QueryMol(f"[{HALIDE}][*:4]")],
        product=QueryMol("[*:1][N:2]([*:3])[*:4]"),
        reaction_id=269982,  # TODO: BOC cleavage
        chemical_space="real",
    ),
    Reaction(
        reactants=[
            QueryMol(r"[H][O:1]/[N:2]=[C:3](\[*:4])[N:5]([H])[H]"),
            QueryMol("[H]O[C:6](=O)[*:7]"),
        ],  # TODO: should this be isomeric?
        product=QueryMol(
            "[*:4][c:3]1[n:2][o:1][c:6]([*:7])[n:5]1"
        ),  # TODO: check if aromatic; BOC cleavage
        reaction_id=270036,
        chemical_space="real",
    ),
    Reaction(
        reactants=[
            QueryMol("[*:1][N:2]([H])[*:3]"),
            QueryMol("O=C(O[C:4](=[O:5])[*:6])[*]"),
        ],
        product=QueryMol("[O:5]=[C:4]([*:6])[N:2]([*:1])[*:3]"),
        reaction_id=270062,
        chemical_space="real",
    ),
    Reaction(
        reactants=[QueryMol("[H][N:1]([H])[*:2]"), QueryMol("[H]O[C:3](=[O:4])[*:5]"),],
        product=QueryMol("[H][N:1]([*:2])[C:3](=[O:4])[*:5]"),
        reaction_id=270112,  # TODO: ester hydrolysis
        chemical_space="real",
    ),
    Reaction(
        reactants=[QueryMol("[*:1][N:2]([H])[*:3]"), QueryMol(f"[{HALIDE}][*:4]")],
        product=QueryMol("[*:1][N:2]([*:3])[*:4]"),
        reaction_id=270122,  # TODO: ester hydrolysis
        chemical_space="real",
    ),
    Reaction(
        reactants=[QueryMol("[*:1][N:2]([H])[*:3]"), QueryMol(f"[{HALIDE}][*:4]")],
        product=QueryMol("[*:1][N:2]([*:3])[*:4]"),
        reaction_id=270166,  # TODO: BOC cleavage
        chemical_space="real",
    ),
    Reaction(
        reactants=[QueryMol("[N:1]#[C:2][*:3]"), QueryMol("[H]O[C:4](=O)[*:5]")],
        product=QueryMol("[*:3][c:2]1no[c:4]([*:5])[n:3]1"),  # TODO: check if aromatic
        reaction_id=270196,
        chemical_space="real",
    ),
    Reaction(
        reactants=[QueryMol("[*:1][N:2]([H])[*:3]"), QueryMol(f"[{HALIDE}][*:4]")],
        product=QueryMol("[*:1][N:2]([*:3])[*:4]"),
        reaction_id=270344,  # TODO: check if okay to reuse ID; ester hydrolysis
        chemical_space="real",
    ),
    Reaction(
        reactants=[QueryMol("[*:1][SH:2]"), QueryMol(f"[{HALIDE}][*:3]")],
        product=QueryMol("[*:1]S[*:3]"),
        reaction_id=270344,  # TODO: check if okay to reuse ID; ester hydrolysis
        chemical_space="real",
    ),
    Reaction(
        reactants=[
            QueryMol(r"[H][O:1]/[N:2]=[C:3](\[*:4])[N:5]([H])[H]")
        ],  # TODO: should this be isomeric?
        product=QueryMol("[H][n:5]1[c:3]([*:4])[n:2][o:1]c1=O"),
        reaction_id=270690,
        chemical_space="real",
    ),
    Reaction(
        reactants=[
            QueryMol("S=[C:1]=[N:2][*:3]"),
            QueryMol("[H][N:4]([*:5])[*:6]"),
            QueryMol("[H][N:7]([H])[N:8]([H])[C:9](=O)[*:10]"),
        ],
        product=QueryMol(
            "[*:5][N:4]([*:6])[c:1]1[n:7][n:8][c:9]([*:10])[n:2]1[*:3]"
        ),  # TODO: should this be aromatic?
        reaction_id=270942,
        chemical_space="real",
    ),
    Reaction(
        reactants=[QueryMol("[O:1]=[C:2]([*:3])O*"), QueryMol("[*:4][N:5]([H])[*:6]"),],
        product=QueryMol("[*:3][C:2]([N:5]([*:4])[*:6])=[O:1]"),
        reaction_id=271144,  # TODO: ester hydrolysis
        chemical_space="real",
    ),
    Reaction(
        reactants=[QueryMol("[N:1]#[C:2][*:3]"), QueryMol("[H]O[C:4](=O)[*:5]")],
        product=QueryMol("[*:3][c:2]1no[c:4]([*:5])[n:3]1"),  # TODO: check if aromatic
        reaction_id=271362,  # TODO: BOC cleavage
        chemical_space="real",
    ),
    Reaction(
        reactants=[QueryMol("[H]OB([*:1])O[H]"), QueryMol(f"[{HALIDE}][*:2]")],
        product=QueryMol("[*:1][*:2]"),
        reaction_id=271570,
        chemical_space="real",
    ),
    Reaction(
        reactants=[
            QueryMol("[H][N:1]([H])[N:2]([H])[C:3](=[O:4])[*:5]"),
            QueryMol("[H]O[C:6](=O)[*:7]"),
        ],
        product=QueryMol(
            "[*:5][c:3]1[n:2][n:1][c:6]([*:7])[o:4]1"
        ),  # TODO: check if aromatic
        reaction_id=271722,
        chemical_space="real",
    ),
    Reaction(
        reactants=[QueryMol(f"[H][C:1]([{HALIDE}])([*:2])[C:3](=O)[*:4]"),],
        product=QueryMol("[H]N([H])c1n[c:3]([*:4])[c:1]([*:2])s1"),
        reaction_id=272150,  # TODO: BOC cleavage
        chemical_space="real",
    ),
    Reaction(
        reactants=[QueryMol("[H][N:1]([H])[*:2]"), QueryMol("[H]O[C:3](=[O:4])[*:5]"),],
        product=QueryMol("[H][N:1]([*:2])[C:3](=[O:4])[*:5]"),
        reaction_id=272270,  # TODO: ester hydrolysis
        chemical_space="real",
    ),
    Reaction(
        reactants=[QueryMol("[H][N:1]([H])[*:2]"), QueryMol("[H]O[C:3](=[O:4])[*:5]")],
        product=QueryMol("[H][N:1]([*:2])[C:3](=[O:4])[*:5]"),
        reaction_id=272430,  # TODO: BOC cleavage
        chemical_space="real",
    ),
    Reaction(
        reactants=[QueryMol("[*:1][N:2]([H])[*:3]"), QueryMol(f"[{HALIDE}][*:4]")],
        product=QueryMol("[*:1][N:2]([*:3])[*:4]"),
        reaction_id=272692,  # TODO: check if okay to reuse ID; BOC cleavage
        chemical_space="real",
    ),
    Reaction(
        reactants=[QueryMol("[*:1][OH:2]"), QueryMol(f"[{HALIDE}][*:3]")],
        product=QueryMol("[*:1]O[*:3]"),
        reaction_id=272692,  # TODO: check if okay to reuse ID; BOC cleavage
        chemical_space="real",
    ),
    Reaction(
        reactants=[QueryMol("[*:1][N:2]([H])[*:3]"), QueryMol(f"[{HALIDE}][*:4]")],
        product=QueryMol("[*:1][N:2]([*:3])[*:4]"),
        reaction_id=272710,  # TODO: check if okay to reuse ID; ester hydrolysis
        chemical_space="real",
    ),
    Reaction(
        reactants=[QueryMol("[*:1][OH,SH:2]"), QueryMol(f"[{HALIDE}][*:3]")],
        product=QueryMol("[*:1][O,S][*:3]"),
        reaction_id=272710,  # TODO: check if okay to reuse ID; ester hydrolysis
        chemical_space="real",
    ),
    Reaction(
        reactants=[QueryMol("[H]OB([*:1])O[H]"), QueryMol(f"[{HALIDE}][*:2]")],
        product=QueryMol("[*:1][*:2]"),
        reaction_id=273584,  # TODO: BOC cleavage
        chemical_space="real",
    ),
    Reaction(
        reactants=[QueryMol("[H][N:1]([H])[*:2]"), QueryMol("[H]O[C:3](=[O:4])[*:5]")],
        product=QueryMol("[H][N:1]([*:2])[C:3](=[O:4])[*:5]"),
        reaction_id=273652,  # TODO: ester hydrolysis
        chemical_space="real",
    ),
    Reaction(
        reactants=[QueryMol("[H]OB([*:1])O[H]"), QueryMol(f"[{HALIDE}][*:2]")],
        product=QueryMol("[*:1][*:2]"),
        reaction_id=273712,
        chemical_space="real",
    ),
    Reaction(
        reactants=[
            QueryMol("[N-:1]=[N+:2]=[N:3][*:4]"),  # TODO: is the azide correct?
            QueryMol("[H][C:5]#[C:6][*:7]"),
        ],
        product=QueryMol(
            "[H][c:5]1[c:6]([*:7])[n:1][n:2][n:3]1[*:4]"
        ),  # TODO: check if aromatic
        reaction_id=273910,
        chemical_space="real",
    ),
    Reaction(
        reactants=[QueryMol(""), QueryMol("")],
        product=QueryMol(""),
        reaction_id=274052,  # TODO: X group with BOC
        chemical_space="real",
    ),
    Reaction(
        reactants=[
            QueryMol("[H][N:1]([H])[C:2](=[S:3])[*:4]"),
            QueryMol(
                f"[H][C:5]([*:6])([{HALIDE}])[C:7](=O)[C:8]#[C:9][Si](C([H])([H])[H])(C([H])([H])[H])C([H])([H])[H]"
            ),
            QueryMol("[N-:10]=[N+:11]=[N:12][*:13]"),  # TODO: is the azide correct?
        ],
        product=QueryMol("[H][c:9]1[c:8](-[c:7]2[n:1][c:2]([*:4])[s:3][c:5]2[*:6])[n:10][n:11][n:12]1[*:13]"),  # TODO: check if aromatic
        reaction_id=274972,
        chemical_space="real",
    ),
    Reaction(
        reactants=[QueryMol(""), QueryMol("")],
        product=QueryMol(""),
        reaction_id=0,
        chemical_space="real",
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
    reaction for reaction in REAL_REACTIONS if reaction.id in MOST_COMMON_REAL_IDS
)
