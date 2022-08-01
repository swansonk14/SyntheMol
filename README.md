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

Install [chem_utils](https://github.com/swansonk14/chem_utils).
```
cd ..
git clone git@github.com:swansonk14/chem_utils.git
cd chem_utils
pip install -e .
cd ../combinatorial_antibiotics
```

Install [chemprop](https://github.com/chemprop/chemprop).
```
cd ..
git clone git@github.com:chemprop/chemprop.git
cd chemprop
pip install -e .
cd ../combinatorial_antibiotics
```


## Download Data

Download the REAL Enamine building blocks SDF file and CXSMILES database files from here: https://enamine.net/compound-collections/real-compounds/real-database#

The `2021q3-4_Enamine_REAL_reagents_SDF.sdf` contains 138,085 molecules.

The 11 CXSMILES files contain 4,492,114,676 molecules and 170 chemical reactions.

## Process Data

### SDF to SMILES

Convert the building blocks from SDF to (unique) SMILES using the `sdf_to_smiles.py` script from [chem_utils](https://github.com/swansonk14/chem_utils).
```
python sdf_to_smiles.py \
    --data_path ../../combinatorial_antibiotics/data/2021q3-4_Enamine_REAL_reagents_SDF.sdf \
    --save_path ../../combinatorial_antibiotics/data/2021q3-4_Enamine_REAL_reagents_SMILES.csv \
    --properties Reagent_ID Catalog_ID
```

All molecules were successfully converted from SDF to SMILES, and among those 134,609 are unique.

Note: This seems to be because the SMILES are not capturing any stereochemistry information even though it is annotated with the `CFG` tag in the SDF file (although 3D coordinates are missing).

### Remove Salts

Remove the salts from the building blocks using the `canonicalize_smiles.py` script from [chem_utils](https://github.com/swansonk14/chem_utils). This will also canonicalize the SMILES using RDKit's canonicalization method.
```
python canonicalize_smiles.py \
    --data_path ../../combinatorial_antibiotics/data/2021q3-4_Enamine_REAL_reagents_SMILES.csv \
    --save_path ../../combinatorial_antibiotics/data/2021q3-4_Enamine_REAL_reagents_SMILES_no_salts.csv \
    --remove_salts \
    --delete_disconnected_mols
```

NOTE: This step is crucial to prevent errors in running reactions. Salts can cause reactions to create products that are the same as the reactants, leading to undesired infinite loops.


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

## Count fragments in training molecules

Count the number of building blocks that appear in each of the training set molecules.
```
python find_fragments_in_molecules.py \
    --fragment_path ../data/2021q3-4_Enamine_REAL_reagents_SMILES_no_salts.csv \
    --molecule_path ../data/Screening_data/AB_training/2500+5000_training_set_2sd.csv \
    --molecule_smiles_column SMILES \
    --counts_save_path ../data/Screening_data/AB_training/2500+5000_training_set_2sd_fragment_counts.csv \
    --plot_save_path ../plots/train_7500_fragment_counts.pdf \
```

## Count number of feasible molecules

Count the number of molecules that could be feasibly produced by the REAL reactions using the building blocks.
```
python count_feasible_molecules.py \
    --fragment_path ../data/2021q3-4_Enamine_REAL_reagents_SMILES_no_salts.csv \
    --reagent_to_fragments_path ../data/reagents_to_fragments.json
```

## Map reagents to fragments

Map reagents (reactants) to REAL fragments (building blocks). This pre-computation step saves time when generating molecules.
```
python map_reagents_to_fragments.py \
    --fragment_path ../data/2021q3-4_Enamine_REAL_reagents_SMILES_no_salts.csv \
    --save_path ../data/reagents_to_fragments.json
```

## Generate random molecules

Generate random molecules using combinatorial molecule construction.
```
python generate_random_molecules.py \
    --fragment_path ../data/2021q3-4_Enamine_REAL_reagents_SMILES_no_salts.csv \
    --reagent_to_fragments_path ../data/reagents_to_fragments.json \
    --num_molecules 10 \
    --max_num_reactions 3 \
    --save_path ../generations/random.csv
```


## Selecting representative Enamine molecules

Download and extract the Enamine screening collection.
```
wge thttps://enamine.net/files/Stock_Screening_Collections/screening_collection.zip
unzip screening_collection.zip
rm screening_collection.zip
```

Convert SDF to SMILES. (Command from [chem_utils](https://github.com/swansonk14/chem_utils).)
```
python sdf_to_smiles.py \
    --data_path ../../combinatorial_antibiotics/data/screening_collection.sdf \
    --save_path ../../combinatorial_antibiotics/data/screening_collection.csv \
    --properties idnumber
```

2,968,443 of 2,968,445 compounds were converted successfully.

Cluster molecules. (Command from [chem_utils](https://github.com/swansonk14/chem_utils).)
```
python cluster_molecules.py \
    --data_path ../../combinatorial_antibiotics/data/screening_collection.csv \
    --save_path ../../combinatorial_antibiotics/data/screening_collection_clustered.csv \
    --num_clusters 5000 \
    --mini_batch
```

Randomly select a molecule from each cluster. (Command from [chem_utils](https://github.com/swansonk14/chem_utils).)
```
python sample_molecules.py \
    --data_path ../../combinatorial_antibiotics/data/screening_collection_clustered.csv \
    --save_path ../../combinatorial_antibiotics/data/screening_collection_clustered_sampled.csv \
    --num_molecules 1 \
    --cluster_column cluster_label
```


## Process AB training data

```
python process_data.py \
    --data_paths ../data/Screening_data/AB_original/AB_2560_normalized.csv ../data/Screening_data/AB_original/AB_Mar27_normalized.csv ../data/Screening_data/AB_original/For_gen_AB_DRH.csv \
    --save_path ../data/Screening_data/AB_combined.csv \
    --save_hits_path ../data/Screening_data/AB_combined_hits.csv
```


## Train random forest model on fingerprints

```
python train_model.py \
    --data_path ../data/Screening_data/AB_combined.csv \
    --save_path ../ckpt/AB_combined_RF_rdkit.pkl \
    --model_type rf \
    --fingerprint_type rdkit
```


## Map fragments to model scores and train similarities

```
python map_fragments_to_model_scores.py \
    --fragment_path ../data/2021q3-4_Enamine_REAL_reagents_SMILES_no_salts.csv \
    --model_path ../ckpt/AB_combined_RF_rdkit.pkl \
    --save_path ../ckpt/AB_combined_RF_rdkit_fragments_to_model_scores.json \
    --model_type rf \
    --fingerprint_type rdkit
```

```
python map_fragments_to_train_similarities.py \
    --fragment_path ../data/2021q3-4_Enamine_REAL_reagents_SMILES_no_salts.csv \
    --train_hits_path ../data/Screening_data/AB_combined_hits.csv \
    --save_path ../data/AB_combined_fragments_to_train_similarities.json
```


## Generate molecules using MCTS

```
python tree_search.py \
    --model_path ../ckpt/AB_combined_RF_rdkit.pkl \
    --fragment_path ../data/2021q3-4_Enamine_REAL_reagents_SMILES_no_salts.csv \
    --fragment_to_model_score_path ../ckpt/AB_combined_RF_rdkit_fragments_to_model_scores.json \
    --fragment_to_train_similarity_path ../data/AB_combined_fragments_to_train_similarities.json \
    --reagent_to_fragments_path ../data/reagents_to_fragments.json \
    --train_hits_path ../data/Screening_data/AB_combined_hits.csv \
    --save_path ../generations/tree_search/mcts.csv \
    --search_type mcts \
    --model_type rf \
    --fingerprint_type rdkit
```


## Assess generated molecules

```
python assess_generated_molecules.py \
    --data_path ../generations/tree_search/mcts.csv \
    --save_dir ../generations/tree_search/mcts \
    --train_hits_path ../data/Screening_data/AB_combined_hits.csv
```

## Select generated molecules

Compute nearest neighbor Tversky distance to train set. (Command from [chem_utils](https://github.com/swansonk14/chem_utils).)

```
python nearest_neighbor.py \
    --data_path ../../combinatorial_antibiotics/generations/tree_search/mcts/molecules.csv \
    --reference_data_path ../../combinatorial_antibiotics/data/Screening_data/AB_combined_hits.csv \
    --reference_name train_hits \
    --metrics tversky
```

Filter to only keep molecules with similarity <= 0.4 and model score >= 0.2 (Command from [chem_utils](https://github.com/swansonk14/chem_utils).)

```
python filter_molecules.py \
    --data_path ../../combinatorial_antibiotics/generations/tree_search/mcts/molecules.csv \
    --save_path ../../combinatorial_antibiotics/generations/tree_search/mcts/molecules_filtered.csv \
    --filter_column train_hits_tversky_nearest_neighbor_similarity \
    --max_value 0.4
```

```
python filter_molecules.py \
    --data_path ../../combinatorial_antibiotics/generations/tree_search/mcts/molecules_filtered.csv \
    --save_path ../../combinatorial_antibiotics/generations/tree_search/mcts/molecules_filtered.csv \
    --filter_column model_score \
    --min_value 0.2
```

Cluster molecules. (Command from [chem_utils](https://github.com/swansonk14/chem_utils).)

```
python cluster_molecules.py \
    --data_path ../../combinatorial_antibiotics/generations/tree_search/mcts/molecules_filtered.csv \
    --num_clusters 100
```

Select top molecule from each cluster. (Command from [chem_utils](https://github.com/swansonk14/chem_utils).)

```
python select_from_clusters.py \
    --data_path ../../combinatorial_antibiotics/generations/tree_search/mcts/molecules_filtered.csv \
    --save_path ../../combinatorial_antibiotics/generations/tree_search/mcts/molecules_selected_100.csv \
    --value_column model_score
```

Visualize selected molecules. (Command from [chem_utils](https://github.com/swansonk14/chem_utils).)

```
python visualize_molecules.py \
    --data_path ../../combinatorial_antibiotics/generations/tree_search/mcts/molecules_selected_100.csv \
    --save_dir ../../combinatorial_antibiotics/generations/tree_search/mcts/molecules_selected_100
```
