# Combinatorial Antibiotic Generation

## Installation

Install conda environment.
```
conda env create -f environment.yml
```

Activate conda environment.
```
conda activate combinatorial_antibiotics
```

## Download Data

Download the REAL Enamine building blocks SDF file from here: https://enamine.net/compound-collections/real-compounds/real-database#

The `2021q3-4_Enamine_REAL_reagents_SDF.sdf` contains 138,085 molecules.

## Process Data

Convert the building blocks from SDF to (unique) SMILES.
```
python sdf_to_smiles.py \
    --data_path ../data/2021q3-4_Enamine_REAL_reagents_SDF.sdf \
    --save_path ../data/2021q3-4_Enamine_REAL_reagents_SMILES.csv
```

All molecules were successfully converted from SDF to SMILES, and among those 134,609 are unique.
