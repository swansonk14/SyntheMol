# Documentation

TODO: put data on Zenodo and clean up google drive

This directory contains instructions for reproducing the results in our antibiotic generation paper: [TODO](TODO).

The relevant data should be downloaded from [this Google Drive file](https://drive.google.com/uc?id=1ziZJqk1PhDtc-kRFxNA3CDM728zjUfHH) to `SyntheMol/data`. This can be done as follows:

```bash
gdown "https://drive.google.com/uc?id=1ziZJqk1PhDtc-kRFxNA3CDM728zjUfHH" -O $(python -c "import synthemol; from pathlib import Path; print(Path(synthemol.__path__[0]).parent)")/data.zip
unzip data.zip && rm data.zip
```

Note that here, we use the 2021 q3-4 version of the building blocks and the 2022 q1-2 version of the enumerated REAL Space molecules. However, by default, SyntheMol now downloads a newer version of the building blocks during installation.

[real.md](real.md): Instructions for processing Enamine REAL building blocks, reactions, and molecules.

[clogp.md](clogp.md): Instructions for performing an _in silico_ study of SyntheMol using a computational molecular property, cLogP, which is the computed octanol-water partition coefficient.

[antibiotics.md](antibiotics.md): Instructions for generating antibiotic candidates for _Acinetobacter baumannii_. Includes instructions for processing antibiotics data, training antibacterial activity prediction models, generating molecules with SyntheMol, and selecting candidates.

[plots.md](plots.md) Instructions for producing plots analyzing the data and results.