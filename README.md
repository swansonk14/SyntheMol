# Combinatorial Antibiotic Generation

## Installation

Download code.
```bash
git clone git@github.com:swansonk14/combinatorial_antibiotics.git
cd combinatorial_antibiotics
```

Install conda environment.
```bash
conda env create -f environment.yml
```

Note: If there are any conda installation or version issues, `environment-frozen.yml` lists the explicit versions of the full list of installed packages used to run this code.

Activate conda environment.
```bash
conda activate combinatorial_antibiotics
```

Install [chem_utils](https://github.com/swansonk14/chem_utils).
```bash
cd ..
git clone git@github.com:swansonk14/chem_utils.git
cd chem_utils
pip install -e .
cd ../combinatorial_antibiotics
```

(Note: Had some issues with pip so commented out the `install_requires` lines in `setup.py`.)

Install [chemprop](https://github.com/chemprop/chemprop).
```bash
cd ..
git clone git@github.com:chemprop/chemprop.git
cd chemprop
pip install -e .
cd ../combinatorial_antibiotics
```

(Note: Had some issues with pip so commented out the `install_requires` and `extras_require` lines in `setup.py`.)


## Process Data

### Download REAL space

Download the reagent and reaction IDs used for the full REAL space of 31B compounds (as of 8/30/22).

```bash
lftp -c "open -p 21 -u username,password ftp-rdb-fr.chem-space.com; mirror -c --parallel=16 . data/Enamine_REAL_space"
```

In the above command, replace `username` and `password` with the appropriate values.


### Map REAL Reactions to Reagents

Determine which reagents are valid in which REAL reactions.

```bash
python map_real_reactions_to_reagents.py \
    --data_dir data/Enamine_REAL_space \
    --save_path data/reaction_to_reagents_REAL_space.json \
    --parallel
```

Total number of molecules = 31,507,987,117


### Count REAL Reactions

Determine which reactions (and reagents) are most common in REAL space.

```bash
python count_real_space.py \
    --data_dir data/Enamine_REAL_space \
    --save_dir data/Enamine_REAL_space_counts
```


### Sample REAL Molecules

Randomly sample REAL space molecules for analysis of a representative sample of REAL space molecules.

```bash
python sample_real_space.py \
    --data_dir data/Enamine_REAL_space \
    --save_path data/Enamine_REAL_space_sampled_25k.csv \
    --num_molecules 25000 \
    --parallel
```


### Download Building Blocks

All data, raw and processed, is available in this Google Drive folder: https://drive.google.com/drive/folders/1sbl1gL1d3acVJ1RZVtJV90uLgW1j6ee9?usp=sharing. Any references to data paths are relative to this directory.

Download the REAL Enamine building blocks SDF file from https://enamine.net/compound-collections/real-compounds/real-database

Note: The building blocks SDF file appears to have been removed from their website but can be found in the Google Drive folder as `2021q3-4_Enamine_REAL_reagents_SDF.sdf`

The `2021q3-4_Enamine_REAL_reagents_SDF.sdf` file contains 138,085 molecules.


### SDF to SMILES

Convert the building blocks from SDF to (unique) SMILES using [chem_utils](https://github.com/swansonk14/chem_utils).
```bash
python sdf_to_smiles.py \
    --data_path ../../combinatorial_antibiotics/data/2021q3-4_Enamine_REAL_reagents_SDF.sdf \
    --save_path ../../combinatorial_antibiotics/data/2021q3-4_Enamine_REAL_reagents_SMILES.csv \
    --properties Reagent_ID Catalog_ID
```

All 138,085 molecules were successfully converted from SDF to SMILES, and among those 134,609 are unique molecules.

Note: This seems to be because the SMILES are not capturing any stereochemistry information even though it is annotated with the `CFG` tag in the SDF file (although 3D coordinates are missing).


### Remove Salts

Remove the salts from the building blocks using [chem_utils](https://github.com/swansonk14/chem_utils). This will also canonicalize the SMILES using RDKit's canonicalization method.
```bash
python canonicalize_smiles.py \
    --data_path ../../combinatorial_antibiotics/data/2021q3-4_Enamine_REAL_reagents_SMILES.csv \
    --save_path ../../combinatorial_antibiotics/data/2021q3-4_Enamine_REAL_reagents_SMILES_no_salts.csv \
    --remove_salts \
    --delete_disconnected_mols
```

This removes 25 molecules whose salts cannot be stripped, leaving 138,060 molecules, of which 132,479 are unique.

Note: This step is crucial to prevent errors in running reactions. Salts can cause reactions to create products that are the same as the reactants, leading to undesired infinite loops during molecule generation.


### Count REAL molecules limited to valid fragments

Count REAL space when limiting to the reactions and fragments that we have post-processing.

```bash
python count_real_space.py \
    --data_dir data/Enamine_REAL_space \
    --save_dir data/Enamine_REAL_space_counts_with_selected_fragments_and_reactions \
    --fragment_path data/2021q3-4_Enamine_REAL_reagents_SMILES_no_salts.csv \
    --only_selected_reactions
```

Number of files = 2,993
Number of fragments = 138,060 (132,479 unique molecules)
Number of selected reactions = 13
Total number of molecules = 31,507,987,117
Total number of molecules with selected fragments/reactions = 29,575,293,692


### Process AB Training Data

The data referred to here is in the `screening_data` subfolder from the Google Drive folder.

```bash
python process_data.py \
    --data_paths data/screening_data/AB_original/AB_2560_normalized.csv data/screening_data/AB_original/AB_Mar27_normalized.csv data/screening_data/AB_original/For_gen_AB_DRH.csv \
    --save_path data/screening_data/AB_combined.csv \
    --save_hits_path data/screening_data/AB_combined_hits.csv
```

Output:
```
AB_2560_normalized
Data size = 2,371
Mean activity = 0.9759337199493885
Std activity = 0.2673771821337539
Activity threshold of mean - 2 std = 0.4411793556818807
Number of hits = 130
Number of non-hits = 2,241

AB_Mar27_normalized
Data size = 5,376
Mean activity = 0.9980530249533108
Std activity = 0.13303517604898704
Activity threshold of mean - 2 std = 0.7319826728553367
Number of hits = 112
Number of non-hits = 5,264

For_gen_AB_DRH
Data size = 6,680
Mean activity = 0.9701813713618264
Std activity = 0.17623388953330232
Activity threshold of mean - 2 std = 0.6177135922952217
Number of hits = 294
Number of non-hits = 6,386

Full data size = 14,427
Data size after dropping non-conflicting duplicates = 13,594
Data size after dropping conflicting duplicates = 13,524

Final data size = 13,524
Number of hits = 470
Number of non-hits = 13,054
```


### Download ChEMBL Antibiotics

Download lists of known antibiotic-related compounds from ChEMBL using the following search terms. For each, click the CSV download button, unzip the downloaded file, and rename the CSV file appropriately.

- [antibacterial](https://www.ebi.ac.uk/chembl/g/#search_results/compounds/query=antibacterial)
- [antibiotic](https://www.ebi.ac.uk/chembl/g/#search_results/compounds/query=antibiotic)

The results of the queries on August 8, 2022, are available in the Google Drive folder in the `chembl` subfolder.

The files are:

- `chembl_antibacterial.csv` with 604 molecules (589 with SMILES)
- `chembl_antibiotic.csv` with 636 molecules (587 with SMILES)

These four files can be combined, processed, and deduplicated to form a single collection of antibiotic-related compounds:

```bash
python merge_chembl_downloads.py \
    --data_paths data/chembl/chembl_antibacterial.csv data/chembl/chembl_antibiotic.csv \
    --labels antibacterial antibiotic \
    --save_path data/chembl/chembl_antibacterial_antibiotic.csv
```

The file `chembl_antibacterial_antibiotic.csv` contains 1,005 molecules.


## Build Models

### Train Models

Train 10 random forest with RDKit features and 10 chemprop models (with or without RDKit features) using 10-fold cross-validation on the AB training data.

Note: For both random forest and chemprop models, we used CPU-only machines (no GPUs).

Random forest with RDKit
```bash
python train_model.py \
    --data_path data/screening_data/AB_combined.csv \
    --save_dir ckpt/AB_combined_RF_rdkit \
    --model_type rf \
    --dataset_type classification \
    --fingerprint_type rdkit \
    --property_column activity \
    --num_models 10
```

Chemprop
```bash
python train_model.py \
    --data_path data/screening_data/AB_combined.csv \
    --save_dir ckpt/AB_combined_chemprop \
    --dataset_type classification \
    --model_type chemprop \
    --property_column activity \
    --num_models 10
```

Chemprop with RDKit
```bash
python train_model.py \
    --data_path data/screening_data/AB_combined.csv \
    --save_dir ckpt/AB_combined_chemprop_rdkit \
    --model_type chemprop \
    --dataset_type classification \
    --fingerprint_type rdkit \
    --property_column activity \
    --num_models 10
```

| Model          | ROC-AUC         | PRC-AUC         |
|----------------|-----------------|-----------------|
| RF RDKit       | 0.835 +/- 0.035 | 0.401 +/- 0.099 |
| Chemprop       | 0.803 +/- 0.036 | 0.354 +/- 0.086 |
| Chemprop RDKit | 0.827 +/- 0.028 | 0.388 +/- 0.078 |


### Map Fragments to Model Scores

Random forest with RDKit
```bash
python map_fragments_to_model_scores.py \
    --fragment_path data/2021q3-4_Enamine_REAL_reagents_SMILES_no_salts.csv \
    --model_path ckpt/AB_combined_RF_rdkit \
    --save_path ckpt/AB_combined_RF_rdkit/fragments_to_model_scores.json \
    --model_type rf \
    --fingerprint_type rdkit
```

Chemprop
```bash
python map_fragments_to_model_scores.py \
    --fragment_path data/2021q3-4_Enamine_REAL_reagents_SMILES_no_salts.csv \
    --model_path ckpt/AB_combined_chemprop \
    --save_path ckpt/AB_combined_chemprop/fragments_to_model_scores.json \
    --model_type chemprop
```

Chemprop with RDKit
```bash
python map_fragments_to_model_scores.py \
    --fragment_path data/2021q3-4_Enamine_REAL_reagents_SMILES_no_salts.csv \
    --model_path ckpt/AB_combined_chemprop_rdkit \
    --save_path ckpt/AB_combined_chemprop_rdkit/fragments_to_model_scores.json \
    --model_type chemprop \
    --fingerprint_type rdkit
```


## Generate Molecules

Run MCTS with the random forest and chemprop models to generate molecules with one reaction.

Random forest with RDKit
```bash
python tree_search.py \
    --model_path ckpt/AB_combined_RF_rdkit \
    --fragment_path data/2021q3-4_Enamine_REAL_reagents_SMILES_no_salts.csv \
    --reaction_to_reagents_path data/reaction_to_reagents_REAL_space.json \
    --fragment_to_model_score_path ckpt/AB_combined_RF_rdkit/fragments_to_model_scores.json \
    --save_dir generations/mcts_AB_combined_RF_rdkit \
    --search_type mcts \
    --model_type rf \
    --fingerprint_type rdkit \
    --n_rollout 20000 \
    --fragment_diversity \
    --max_reactions 1
```

Chemprop
```bash
python tree_search.py \
    --model_path ckpt/AB_combined_chemprop \
    --fragment_path data/2021q3-4_Enamine_REAL_reagents_SMILES_no_salts.csv \
    --reaction_to_reagents_path data/reaction_to_reagents_REAL_space.json \
    --fragment_to_model_score_path ckpt/AB_combined_chemprop/fragments_to_model_scores.json \
    --save_dir generations/mcts_AB_combined_chemprop \
    --search_type mcts \
    --model_type chemprop \
    --n_rollout 20000 \
    --fragment_diversity \
    --max_reactions 1
```

Chemprop with RDKit
```bash
python tree_search.py \
    --model_path ckpt/AB_combined_chemprop_rdkit \
    --fragment_path data/2021q3-4_Enamine_REAL_reagents_SMILES_no_salts.csv \
    --reaction_to_reagents_path data/reaction_to_reagents_REAL_space.json \
    --fragment_to_model_score_path ckpt/AB_combined_chemprop_rdkit/fragments_to_model_scores.json \
    --save_dir generations/mcts_AB_combined_chemprop_rdkit \
    --search_type mcts \
    --model_type chemprop \
    --fingerprint_type rdkit \
    --n_rollout 20000 \
    --fragment_diversity \
    --max_reactions 1
```

Random
```bash
python tree_search.py \
    --fragment_path data/2021q3-4_Enamine_REAL_reagents_SMILES_no_salts.csv \
    --reaction_to_reagents_path data/reaction_to_reagents_REAL_space.json \
    --save_dir generations/random \
    --search_type random \
    --n_rollout 20000 \
    --max_reactions 1
```


## Make Predictions on Random Molecules

Use all three models in ensemble to make predictions on the randomly generated molecules. First, each model is used separately to make predictions. Then the average prediction is computed.

Random forest on random molecules
```bash
python predict_model.py \
    --data_path generations/random/molecules.csv \
    --model_path ckpt/AB_combined_RF_rdkit \
    --model_type rf \
    --fingerprint_type rdkit \
    --average_preds
```

Chemprop on random molecules
```bash
python predict_model.py \
    --data_path generations/random/molecules.csv \
    --model_path ckpt/AB_combined_chemprop \
    --model_type chemprop \
    --average_preds
```

Chemprop with RDKit on random molecules
```bash
python predict_model.py \
    --data_path generations/random/molecules.csv \
    --model_path ckpt/AB_combined_chemprop_rdkit \
    --model_type chemprop \
    --fingerprint_type rdkit \
    --average_preds
```

Average model predictions
```bash
python average_scores.py \
    --data_path generations/random/molecules.csv \
    --score_columns rf_rdkit_ensemble_preds chemprop_ensemble_preds chemprop_rdkit_ensemble_preds \
    --average_column score
```


## Assess Generated Molecules

Optionally, assess various quantities about the generated molecules.

```bash
#!/bin/bash

for NAME in mcts_AB_combined_RF_rdkit mcts_AB_combined_chemprop mcts_AB_combined_chemprop_rdkit random
do
python assess_generated_molecules.py \
    --data_path generations/${NAME}/molecules.csv \
    --save_dir generations/${NAME} \
    --reference_paths data/screening_data/AB_combined_hits.csv data/chembl/chembl_antibacterial_antibiotic.csv
done
```


## Select Generated Molecules

Run a series of filtering steps to select molecules for experimental testing


### Train Hit Nearest Neighbor

Compute the nearest neighbor Tversky similarity to train hits using [chem_utils](https://github.com/swansonk14/chem_utils).

```bash
#!/bin/bash

for NAME in mcts_AB_combined_RF_rdkit mcts_AB_combined_chemprop mcts_AB_combined_chemprop_rdkit random
do
python nearest_neighbor.py \
    --data_path ../../combinatorial_antibiotics/generations/${NAME}/molecules.csv \
    --reference_data_path ../../combinatorial_antibiotics/data/screening_data/AB_combined_hits.csv \
    --reference_name train_hits \
    --metrics tversky
done
```


### ChEMBL Antibiotic Nearest Neighbor

Compute the nearest neighbor Tversky similarity to ChEMBL antibiotics using [chem_utils](https://github.com/swansonk14/chem_utils).

```bash
#!/bin/bash

for NAME in mcts_AB_combined_RF_rdkit mcts_AB_combined_chemprop mcts_AB_combined_chemprop_rdkit random
do
python nearest_neighbor.py \
    --data_path ../../combinatorial_antibiotics/generations/${NAME}/molecules.csv \
    --reference_data_path ../../combinatorial_antibiotics/data/chembl/chembl_antibacterial_antibiotic.csv \
    --reference_name chembl_antibacterial_antibiotic \
    --metrics tversky
done
```


### Filter by Similarity to Train Hits

Filter to only keep molecules with nearest neighbor Tverksy similarity to the train hits <= 0.5 using [chem_utils](https://github.com/swansonk14/chem_utils).

Random forest with RDKit
```bash
#!/bin/bash

for NAME in mcts_AB_combined_RF_rdkit mcts_AB_combined_chemprop mcts_AB_combined_chemprop_rdkit random
do
python filter_molecules.py \
    --data_path ../../combinatorial_antibiotics/generations/${NAME}/molecules.csv \
    --save_path ../../combinatorial_antibiotics/generations/${NAME}/molecules_train_sim_below_0.5.csv \
    --filter_column train_hits_tversky_nearest_neighbor_similarity \
    --max_value 0.5
done
```


### Filter by Similarity to ChEMBL Antibiotics

Filter to only keep molecules with nearest neighbor Tverksy similarity to the ChEMBL antibiotics <= 0.5 using [chem_utils](https://github.com/swansonk14/chem_utils).

```bash
#!/bin/bash

for NAME in mcts_AB_combined_RF_rdkit mcts_AB_combined_chemprop mcts_AB_combined_chemprop_rdkit random
do
python filter_molecules.py \
    --data_path ../../combinatorial_antibiotics/generations/${NAME}/molecules_train_sim_below_0.5.csv \
    --save_path ../../combinatorial_antibiotics/generations/${NAME}/molecules_train_sim_below_0.5_chembl_sim_below_0.5.csv \
    --filter_column chembl_antibacterial_antibiotic_tversky_nearest_neighbor_similarity \
    --max_value 0.5
done
```


### Filter by Model Score

Filter to only keep molecules with the top 20% of model scores using [chem_utils](https://github.com/swansonk14/chem_utils).

```bash
#!/bin/bash

for NAME in mcts_AB_combined_RF_rdkit mcts_AB_combined_chemprop mcts_AB_combined_chemprop_rdkit random
do
python filter_molecules.py \
    --data_path ../../combinatorial_antibiotics/generations/${NAME}/molecules_train_sim_below_0.5_chembl_sim_below_0.5.csv \
    --save_path ../../combinatorial_antibiotics/generations/${NAME}/molecules_train_sim_below_0.5_chembl_sim_below_0.5_top_20_percent.csv \
    --filter_column score \
    --top_proportion 0.2
done
```


### Cluster Molecules

Cluster molecules based on their Morgan fingerprint using [chem_utils](https://github.com/swansonk14/chem_utils).

```bash
#!/bin/bash

for NAME in mcts_AB_combined_RF_rdkit mcts_AB_combined_chemprop mcts_AB_combined_chemprop_rdkit random
do
python cluster_molecules.py \
    --data_path ../../combinatorial_antibiotics/generations/${NAME}/molecules_train_sim_below_0.5_chembl_sim_below_0.5_top_20_percent.csv \
    --num_clusters 50
done
```


### Select Molecules from Clusters

Select the top scoring molecule from each cluster using [chem_utils](https://github.com/swansonk14/chem_utils).

```bash
#!/bin/bash

for NAME in mcts_AB_combined_RF_rdkit mcts_AB_combined_chemprop mcts_AB_combined_chemprop_rdkit random
do
python select_from_clusters.py \
    --data_path ../../combinatorial_antibiotics/generations/${NAME}/molecules_train_sim_below_0.5_chembl_sim_below_0.5_top_20_percent.csv \
    --save_path ../../combinatorial_antibiotics/generations/${NAME}/molecules_train_sim_below_0.5_chembl_sim_below_0.5_top_20_percent_selected_50.csv \
    --value_column score
done
```

We now have 50 molecules selected from each model's search that meet our desired train hit similarity and model score criteria.


### Visualize Molecules

Optionally, visualize the selecting molecules using [chem_utils](https://github.com/swansonk14/chem_utils).

```bash
#!/bin/bash

for NAME in mcts_AB_combined_RF_rdkit mcts_AB_combined_chemprop mcts_AB_combined_chemprop_rdkit random
do
python visualize_molecules.py \
    --data_path ../../combinatorial_antibiotics/generations/${NAME}/molecules_train_sim_below_0.5_chembl_sim_below_0.5_top_20_percent_selected_50.csv \
    --save_dir ../../combinatorial_antibiotics/generations/${NAME}
mv 1.png molecules_train_sim_below_0.5_chembl_sim_below_0.5_top_20_percent_selected_50.png
done
```


### Map Molecules to REAL IDs

Map generated molecules to REAL IDs in the format expected by Enamine.

```bash
#!/bin/bash

for NAME in mcts_AB_combined_RF_rdkit mcts_AB_combined_chemprop mcts_AB_combined_chemprop_rdkit random
do
python map_generated_molecules_to_real_ids.py \
    --data_path generations/${NAME}/molecules_train_sim_below_0.5_chembl_sim_below_0.5_top_20_percent_selected_50.csv \
    --save_dir generations/${NAME}/molecules_train_sim_below_0.5_chembl_sim_below_0.5_top_20_percent_selected_50_real_ids
done
```


### Assess Selected Molecules

Optionally, assess various quantities about the selected molecules.

```bash
#!/bin/bash

for NAME in mcts_AB_combined_RF_rdkit mcts_AB_combined_chemprop mcts_AB_combined_chemprop_rdkit random
do
python assess_generated_molecules.py \
    --data_path generations/${NAME}/molecules_train_sim_below_0.5_chembl_sim_below_0.5_top_20_percent_selected_50.csv \
    --save_dir generations/${NAME}/analysis_molecules_train_sim_below_0.5_chembl_sim_below_0.5_top_20_percent_selected_50 \
    --reference_paths data/screening_data/AB_combined_hits.csv data/chembl/chembl_antibacterial_antibiotic.csv
done
```

## Apply Toxicity Model

ClinTox is a clinical trial toxicity binary classification dataset with 1,478 molecules with 112 (7.58%) toxic.

Train toxicity model using the `CT_TOX` property of the ClinTox dataset. (Note that the ClinTox properties `CT_TOX` and `FDA_APPROVED` are nearly always identical, so we only use one.)
```bash
python train_model.py \
    --data_path data/clintox.csv \
    --save_dir ckpt/clintox \
    --model_type chemprop \
    --dataset_type classification \
    --fingerprint_type rdkit \
    --property_column CT_TOX \
    --num_models 10
```

```
Overall test ROC-AUC = 0.881 +/- 0.045
Overall test PRC-AUC = 0.514 +/- 0.141
```

Make predictions on synthesized molecules.
```bash
python predict_model.py \
    --data_path generations/enamine_synthesized_results_formatted.csv \
    --model_path ckpt/clintox \
    --model_type chemprop \
    --fingerprint_type rdkit \
    --average_preds
```
