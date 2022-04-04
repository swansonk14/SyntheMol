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

Note: This seems to be because the SMILES are not capturing any stereochemistry information even though it is annotated with the `CFG` tag in the SDF file (although 3D coordinates are missing).


## Train chemprop on screening data

The following commands are from [Chemprop](https://github.com/chemprop/chemprop).


### Generate RDKit features

Generate RDKit features for 7,500 AB set.
```
python scripts/save_features.py \
    --data_path ../combinatorial_antibiotics/data/Screening_data/AB_training/2500+5000_training_set_2sd.csv \
    --save_path ../combinatorial_antibiotics/features/Screening_data/AB_training/2500+5000_training_set_2sd.npz
```

Generate RDKit features for 2,500 AB+EC+PA+SA set.
```
python scripts/save_features.py \
    --data_path ../combinatorial_antibiotics/data/Screening_data/AB+EC+PA+SA/isomeric_master.csv \
    --save_path ../combinatorial_antibiotics/features/Screening_data/AB+EC+PA+SA/isomeric_master.npz
```


### Train Chemprop models

Train Chemprop model on 7,500 AB set.
```
python train.py \
    --data_path ../combinatorial_antibiotics/data/Screening_data/AB_training/2500+5000_training_set_2sd.csv \
    --dataset_type classification \
    --smiles_column SMILES \
    --target_columns Activity \
    --features_path ../combinatorial_antibiotics/features/Screening_data/AB_training/2500+5000_training_set_2sd.npz \
    --no_features_scaling \
    --split_type cv \
    --num_folds 10 \
    --save_dir ../combinatorial_antibiotics/ckpt/AB_training \
    --quiet
```

Train Chemprop model on 2,500 AB+EC+PA+SA set (multitask).
```
python train.py \
    --data_path ../combinatorial_antibiotics/data/Screening_data/AB+EC+PA+SA/isomeric_master.csv \
    --dataset_type classification \
    --smiles_column smiles \
    --target_columns AB_Activity EC_Activity PA_Activity RN4220Staph_Activity \
    --features_path ../combinatorial_antibiotics/features/Screening_data/AB+EC+PA+SA/isomeric_master.npz \
    --no_features_scaling \
    --split_type cv \
    --num_folds 10 \
    --save_dir ../combinatorial_antibiotics/ckpt/AB+EC+PA+SA \
    --quiet
```

Train Chemprop model on 2,500 AB set.
```
python train.py \
    --data_path ../combinatorial_antibiotics/data/Screening_data/AB+EC+PA+SA/isomeric_master.csv \
    --dataset_type classification \
    --smiles_column smiles \
    --target_columns AB_Activity \
    --features_path ../combinatorial_antibiotics/features/Screening_data/AB+EC+PA+SA/isomeric_master.npz \
    --no_features_scaling \
    --split_type cv \
    --num_folds 10 \
    --save_dir ../combinatorial_antibiotics/ckpt/AB \
    --quiet
```

Train Chemprop model on 2,500 EC set.
```
python train.py \
    --data_path ../combinatorial_antibiotics/data/Screening_data/AB+EC+PA+SA/isomeric_master.csv \
    --dataset_type classification \
    --smiles_column smiles \
    --target_columns EC_Activity \
    --features_path ../combinatorial_antibiotics/features/Screening_data/AB+EC+PA+SA/isomeric_master.npz \
    --no_features_scaling \
    --split_type cv \
    --num_folds 10 \
    --save_dir ../combinatorial_antibiotics/ckpt/EC \
    --quiet
```

Train Chemprop model on 2,500 PA set.
```
python train.py \
    --data_path ../combinatorial_antibiotics/data/Screening_data/AB+EC+PA+SA/isomeric_master.csv \
    --dataset_type classification \
    --smiles_column smiles \
    --target_columns PA_Activity \
    --features_path ../combinatorial_antibiotics/features/Screening_data/AB+EC+PA+SA/isomeric_master.npz \
    --no_features_scaling \
    --split_type cv \
    --num_folds 10 \
    --save_dir ../combinatorial_antibiotics/ckpt/PA \
    --quiet
```

Train Chemprop model on 2,500 SA set.
```
python train.py \
    --data_path ../combinatorial_antibiotics/data/Screening_data/AB+EC+PA+SA/isomeric_master.csv \
    --dataset_type classification \
    --smiles_column smiles \
    --target_columns RN4220Staph_Activity \
    --features_path ../combinatorial_antibiotics/features/Screening_data/AB+EC+PA+SA/isomeric_master.npz \
    --no_features_scaling \
    --split_type cv \
    --num_folds 10 \
    --save_dir ../combinatorial_antibiotics/ckpt/SA \
    --quiet
```
