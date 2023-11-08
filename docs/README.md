# Documentation

This directory contains instructions for reproducing the results in our two antibiotic generation papers: [TODO](TODO) and [TODO](TODO).

## SyntheMol-MCTS

Below are instruments for reproducing the results from the paper [TODO](TODO), which uses the Monte Carlo tree search (MCTS) version of SyntheMol.

The relevant data should be downloaded from [this Google Drive file](https://drive.google.com/uc?id=1ziZJqk1PhDtc-kRFxNA3CDM728zjUfHH) to `SyntheMol/data`. This can be done as follows:

```bash
gdown "https://drive.google.com/uc?id=1ziZJqk1PhDtc-kRFxNA3CDM728zjUfHH" -O $(python -c "import synthemol; from pathlib import Path; print(Path(synthemol.__path__[0]).parent)")/data.zip
unzip data.zip && rm data.zip
```

Note that here, we use the 2021 q3-4 version of the building blocks and the 2022 q1-2 version of the enumerated REAL Space molecules. However, by default, SyntheMol now downloads a newer version of the building blocks during installation.

[real_2022.md](real_2022.md): Instructions for processing Enamine REAL building blocks, reactions, and molecules.

[clogp.md](clogp.md): Instructions for performing an _in silico_ study of SyntheMol using a computational molecular property, cLogP, which is the computed octanol-water partition coefficient.

[antibiotics_mcts.md](antibiotics_mcts.md): Instructions for generating antibiotic candidates for _Acinetobacter baumannii_ using a Monte Carlo tree search (MCTS). Includes instructions for processing antibiotics data, training antibacterial activity prediction models, generating molecules with SyntheMol-MCTS, and selecting candidates.

[plots.md](plots.md) Instructions for producing plots analyzing the data and results.


## SyntheMol-RL

Below are instruments for reproducing the results from the paper [TODO](TODO), which uses the reinforcement learning (RL) version of SyntheMol.

The relevant data should be downloaded from this Google Drive file (TODO: link) to `SyntheMol/data`. This can be done as follows:

```bash
gdown "https://drive.google.com/drive/folders/1ssJhy0ZZoh4P-ELG0xiG5iO53_Jts9du?usp=sharing" --folder -O $(python -c "import synthemol; from pathlib import Path; print(Path(synthemol.__path__[0]).parent)")/rl
```

Note that here, we use the 2022 q1-2 version of the building blocks and the September 25, 2023 release of the enumerated REAL Space molecules.

[real_2023.md](real_2023.md): Instructions for processing Enamine REAL building blocks, reactions, and molecules.

[wuxi.md](wuxi.md): Instructions for processing the WuXi building blocks, reactions, and molecules (May 2023 release).

[antibiotics_rl.md](antibiotics_rl.md): Instructions for generating antibiotic candidates for _Acinetobacter baumannii_ using reinforcement learning (RL). Includes instructions for processing antibiotics data, training antibacterial activity prediction models, generating molecules with SyntheMol-RL, and selecting candidates.
