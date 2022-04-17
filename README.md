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

Download the REAL Enamine building blocks SDF file and CXSMILES database files from here: https://enamine.net/compound-collections/real-compounds/real-database#

The `2021q3-4_Enamine_REAL_reagents_SDF.sdf` contains 138,085 molecules.

The 11 CXSMILES files contain 4,492,114,676 molecules and 170 chemical reactions.

## Process Data

Convert the building blocks from SDF to (unique) SMILES using the `sdf_to_smiles.py` script from [chem_utils](https://github.com/swansonk14/chem_utils).
```
python sdf_to_smiles.py \
    --data_path ../data/2021q3-4_Enamine_REAL_reagents_SDF.sdf \
    --save_path ../data/2021q3-4_Enamine_REAL_reagents_SMILES.csv \
    --properties Reagent_ID Catalog_ID
```

All molecules were successfully converted from SDF to SMILES, and among those 134,609 are unique.

Note: This seems to be because the SMILES are not capturing any stereochemistry information even though it is annotated with the `CFG` tag in the SDF file (although 3D coordinates are missing).


## Plot regression values

```
python plot_regression_values.py \
    --data_path ../../../chemistry/antibiotic_moa/data/broad_organism.csv \
    --rep1_column Acinetobacter_baumannii_50uM_R1 \
    --rep2_column Acinetobacter_baumannii_50uM_R2 \
    --remove_outliers \
    --save_dir ../plots/Acinetobacter_baumannii
```

Repeat for `Acinetobacter_baumannii`, `Escherichia_coli`, `Pseudomonas_aeruginosa`, and `Staphylococcus_aureus`.


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
    --metric prc-auc \
    --extra_metrics auc \
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
    --metric prc-auc \
    --extra_metrics auc \
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
    --metric prc-auc \
    --extra_metrics auc \
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
    --metric prc-auc \
    --extra_metrics auc \
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
    --metric prc-auc \
    --extra_metrics auc \
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
    --metric prc-auc \
    --extra_metrics auc \
    --save_dir ../combinatorial_antibiotics/ckpt/SA \
    --quiet
```

## Count REAL reactions and reagents

Count all REAL reactions and reagents in the 4.5 billion REAL database.
```
python count_real_database.py \
    --data_dir ../data/Enamine_REAL_SMILES \
    --save_dir ../data/Enamine_REAL_counts
```

Count the top 10 unique reactions and reagents in the 4.5 billion REAL database. (Note: Reactions 22, 11, 527, and 240690 represent the same reactants + products, so we consider them one reaction.)
```
python count_real_database.py \
    --data_dir ../data/Enamine_REAL_SMILES \
    --reactions 275592 22 11 527 240690 2430 2708 2230 2718 40 1458 271948 27 \
    --save_dir ../data/Enamine_REAL_counts_top_10
```

## Sample REAL database molecules

Sample REAL database molecules for quick analysis, testing, and visualization.
```
python sample_real_database.py \
    --data_dir ../data/Enamine_REAL_SMILES \
    --num_molecules 10000 \
    --save_dir ../data/Enamine_REAL_SMILES_sampled
```
