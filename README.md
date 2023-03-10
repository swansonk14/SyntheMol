# SyntheMol: Generative AI for Drug Discovery

SyntheMol is a generative AI method for designing structurally novel and diverse drug candidates with predicted bioactivity that are easy to synthesize.

SyntheMol consists of a Monte Carlo tree search (MCTS) that explores a combinatorial chemical space consisting of molecular building blocks and chemical reactions. The MCTS is guided by a molecular property prediction AI model, such as a graph neural network or a random forest. Currently, SyntheMol is designed to use 139,444 building blocks and 13 chemical reactions from the [Enamine REAL Space](https://enamine.net/compound-collections/real-compounds/real-space-navigator), which can produce over 30 billion molecules. However, SyntheMol can be easily adapted to use any set of building blocks and reactions.

SyntheMol is described in the paper [TODO](TODO), where we applied SyntheMol to design novel antibiotics for the Gram-negative bacterium _Acinetobacter baumannii_. Full detail for reproducing the results in the paper are provided in the [docs](https://github.com/swansonk14/SyntheMol/docs) directory.

If you use SyntheMol in your research, please cite:

```
TODO
```

## Table of contents

TODO: create table of contents


TODO: do we want to provide data directly?

## Installation

TODO: put SyntheMol (and chem_utils) on pip

TODO: remove environment.yml?

Install SyntheMol via pip.
```bash
pip install SyntheMol
```

Install SyntheMol from GitHub.
```bash
git clone git@github.com:swansonk14/SyntheMol.git
cd SyntheMol
pip install -e .
```

Download the necessary data files.
```bash
gdown "https://drive.google.com/drive/folders/1LLLwxe_nQAnsRSQpIRq_ngyCm1txS-Sq" -O SyntheMol/files --folder
```


## Combinatorial chemical space

SyntheMol is currently designed to use 139,444 unique building blocks and 13 chemical reactions from the [Enamine REAL Space](https://enamine.net/compound-collections/real-compounds/real-space-navigator), which can produce over 30 billion molecules.

TODO: check numbers with new building blocks

However, an alternate combinatorial chemical space can be used by replacing the building blocks and chemical reactions as follows.
- **Building blocks:** Replace `SyntheMol/files/building_blocks.csv.gz` with a custom file containing the building blocks. The file should be a CSV file with a header row and two columns: `smiles` and `ID`. The `smiles` column should contain the SMILES string for each building block, and the `ID` column should contain a unique ID for each building block.
- **Chemical reactions:** In `SyntheMol/reactions/custom.py`, replace set `CUSTOM_REACTIONS` to a list of `Reaction` objects similar to the `REAL_REACTIONS` list in `SyntheMol/reactions/real.py`. If `CUSTOM_REACTIONS` is not None, then it will automatically be used instead of the `REAL_REACTIONS`.



## Property prediction model

TODO


### Train property prediction model

TODO


### Map building blocks to model scores

TODO


## Generate molecules

TODO

## Filter generated molecules

TODO
