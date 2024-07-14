# SyntheMol-RL

Below are instruments for reproducing the results from the paper [TODO](TODO), which uses the reinforcement learning (
RL) version of SyntheMol.

## Data

The relevant data should be downloaded from this Google Drive file (TODO: link) to `SyntheMol/data`. This can be done as
follows:

```bash
gdown "https://drive.google.com/drive/folders/1ssJhy0ZZoh4P-ELG0xiG5iO53_Jts9du?usp=sharing" --folder -O $(python -c "import synthemol; from pathlib import Path; print(Path(synthemol.__path__[0]).parent)")/rl
```

Note that here, we use the 2022 q1-2 version of the building blocks and the 2022 q1-2 version of the enumerated REAL
Space molecules.

## Instructions

[real.md](real.md): Instructions for processing Enamine REAL building blocks, reactions, and molecules.

[wuxi.md](wuxi.md): Instructions for processing the WuXi building blocks, reactions, and molecules (May 2023 release).

[antibiotics.md](antibiotics.md): Instructions for generating antibiotic candidates for _Staphyloccocus aureus_ using
reinforcement learning (RL). Includes instructions for processing antibiotics data, training antibacterial activity
prediction models, generating molecules with SyntheMol-RL, and selecting candidates.

[ablations.md](ablations.md): Instructions for reproducing the ablation experiments in the paper to determine the
importance of various components of SyntheMol-RL.

[gflownet.md](gflownet.md): Instructions for reproducing the GFlowNet experiments in the paper.