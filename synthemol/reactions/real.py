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
        reactants=[QueryMol("[*:1][C:2]([H])=O"), QueryMol("[H]O[C:3]([C:4]([H])([H])[N:5]([H])[C:6]([*:7])=[O:8])=[O:9]")],
        product=QueryMol("[H]/[C:2](=[C:4]1/[N:5]=[C:6]([*:7])[O:8][C:3]1=[O:9])[*:1]"),  # TODO: should this be an isomeric SMILES? Should it be aromatic?
        reaction_id=10,
        chemical_space="real",
    ),
    Reaction(
        reactants=[QueryMol("[H]C([H])([H])C([H])([H])O[CH:1]=[C:2]([*:3])[*:4]"), QueryMol("[NH2:5][*:6]")],
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
        reactants=[QueryMol("[H][N:1]([*:2])[C:3](=[S:4])[N:5]([H])[N:6]([H])[H]"), QueryMol(f"[H][C:7]([{HALIDE}])([*:8])[C:9](=O)[*:10]")],
        product=QueryMol("[H][N:1]([*:2])[C:3]1=[N:5][N:6]=[C:9]([*:10])[C:7]([H])([*:8])[S:4]1"),
        reaction_id=41,
        chemical_space="real",
    ),
    Reaction(
        reactants=[QueryMol("[H][C:1]([H])([*:2])[C:3](=O)[*:4]"), QueryMol("")],  # TODO: figure out R group ring
        product=QueryMol(""),  # TODO: figure out R group ring
        reaction_id=42,
        chemical_space="real",
    ),
    Reaction(
        reactants=[QueryMol(f"[H][C:1]([{HALIDE}])([*:2])[C:3](=O)[*:4]")],
        product=QueryMol("[H]N([H])c1n[c:3]([*:4])[c:1]([*:2])s1"),  # TODO: check that this is aromatic
        reaction_id=43,
        chemical_space="real",
    ),
    Reaction(
        reactants=[QueryMol("[*:1][N:2]([H])[*:3]"), QueryMol(f"[{HALIDE}][*:4]")],
        product=QueryMol("[*:1][N:2]([*:3])[*:4]"),
        reaction_id=44,
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
