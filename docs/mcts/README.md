# SyntheMol-MCTS

This directory contains instructions for reproducing the results from the paper [TODO](TODO), which uses the Monte Carlo tree search (MCTS) version of SyntheMol.


## Data

The relevant data should be downloaded from [this Zenodo record](https://zenodo.org/doi/10.5281/zenodo.10257838), unzipped, and moved to `SyntheMol/data`.

Note that here, we use the 2021 q3-4 version of the building blocks and the 2022 q1-2 version of the enumerated REAL Space molecules. However, by default, SyntheMol now downloads a newer version of the building blocks during installation.


## Instructions

[real.md](real.md): Instructions for processing Enamine REAL building blocks, reactions, and molecules.

[clogp.md](clogp.md): Instructions for performing an _in silico_ study of SyntheMol using a computational molecular property, cLogP, which is the computed octanol-water partition coefficient.

[antibiotics.md](antibiotics.md): Instructions for generating antibiotic candidates for _Acinetobacter baumannii_ using a Monte Carlo tree search (MCTS). Includes instructions for processing antibiotics data, training antibacterial activity prediction models, generating molecules with SyntheMol-MCTS, and selecting candidates.

[plots.md](plots.md) Instructions for producing plots analyzing the data and results.