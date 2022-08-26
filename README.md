# Combinatorial Antibiotic Generation

## Installation

Install conda environment.
```
conda env create -f environment.yml
```

Note: If there are any conda installation or version issues, `environment-frozen.yml` lists the explicit versions of the full list of installed packages used to run this code.

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

(Note: Had some issues with pip so commented out the `install_requires` lines in `setup.py`.)

Install [chemprop](https://github.com/chemprop/chemprop).
```
cd ..
git clone git@github.com:chemprop/chemprop.git
cd chemprop
pip install -e .
cd ../combinatorial_antibiotics
```

(Note: Had some issues with pip so commented out the `install_requires` and `extras_require` lines in `setup.py`.)


## Process Data

### Download Building Blocks

All data, raw and processed, is available in this Google Drive folder: https://drive.google.com/drive/folders/1sbl1gL1d3acVJ1RZVtJV90uLgW1j6ee9?usp=sharing. Any references to data paths are relative to this directory.

Download the REAL Enamine building blocks SDF file from https://enamine.net/compound-collections/real-compounds/real-database

Note: The building blocks SDF file appears to have been removed from their website but can be found in the Google Drive folder as `2021q3-4_Enamine_REAL_reagents_SDF.sdf`

The `2021q3-4_Enamine_REAL_reagents_SDF.sdf` file contains 138,085 molecules.


### SDF to SMILES

Convert the building blocks from SDF to (unique) SMILES using [chem_utils](https://github.com/swansonk14/chem_utils).
```
python sdf_to_smiles.py \
    --data_path ../../combinatorial_antibiotics/data/2021q3-4_Enamine_REAL_reagents_SDF.sdf \
    --save_path ../../combinatorial_antibiotics/data/2021q3-4_Enamine_REAL_reagents_SMILES.csv \
    --properties Reagent_ID Catalog_ID
```

All molecules were successfully converted from SDF to SMILES, and among those 134,609 are unique.

Note: This seems to be because the SMILES are not capturing any stereochemistry information even though it is annotated with the `CFG` tag in the SDF file (although 3D coordinates are missing).


### Remove Salts

Remove the salts from the building blocks using [chem_utils](https://github.com/swansonk14/chem_utils). This will also canonicalize the SMILES using RDKit's canonicalization method.
```
python canonicalize_smiles.py \
    --data_path ../../combinatorial_antibiotics/data/2021q3-4_Enamine_REAL_reagents_SMILES.csv \
    --save_path ../../combinatorial_antibiotics/data/2021q3-4_Enamine_REAL_reagents_SMILES_no_salts.csv \
    --remove_salts \
    --delete_disconnected_mols
```

This removes 25 molecules whose salts cannot be stripped, leaving 138,060 molecules.

Note: This step is crucial to prevent errors in running reactions. Salts can cause reactions to create products that are the same as the reactants, leading to undesired infinite loops during molecule generation.


### Map Reagents to Fragments

Map reagents (reactants) to REAL fragments (building blocks). This pre-computation step saves time when generating molecules.
```
python map_reagents_to_fragments.py \
    --fragment_path data/2021q3-4_Enamine_REAL_reagents_SMILES_no_salts.csv \
    --save_path data/reagents_to_fragments.json
```


### Process AB Training Data

The data referred to here is in the `screening_data` subfolder from the Google Drive folder.

```
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

The results of the queries on August 8, 2022 are available in the Google Drive folder in the `chembl` subfolder.

The files are:

- `antibacterial.csv` with 604 molecules
- `antibiotic.csv` with 636 molecules

These four files can be combined, processed, and deduplicated to form a single collection of antibiotic-related compounds:

```
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

Random forest
```
python train_model.py \
    --data_path data/screening_data/AB_combined.csv \
    --save_dir ckpt/AB_combined_RF_rdkit \
    --model_type rf \
    --fingerprint_type rdkit \
    --num_models 10
```

Chemprop
```
python train_model.py \
    --data_path data/screening_data/AB_combined.csv \
    --save_dir ckpt/AB_combined_chemprop \
    --model_type chemprop \
    --num_models 10
```

Chemprop with RDKit
```
python train_model.py \
    --data_path data/screening_data/AB_combined.csv \
    --save_dir ckpt/AB_combined_chemprop_rdkit \
    --model_type chemprop \
    --fingerprint_type rdkit \
    --num_models 10
```


### Map Fragments to Model Scores

TODO: left off with reproducibility here.

Random forest
```
python map_fragments_to_model_scores.py \
    --fragment_path data/2021q3-4_Enamine_REAL_reagents_SMILES_no_salts.csv \
    --model_path ckpt/AB_combined_RF_rdkit \
    --save_path ckpt/AB_combined_RF_rdkit/fragments_to_model_scores.json \
    --model_type rf \
    --fingerprint_type rdkit
```

Chemprop
```
python map_fragments_to_model_scores.py \
    --fragment_path data/2021q3-4_Enamine_REAL_reagents_SMILES_no_salts.csv \
    --model_path ckpt/AB_combined_chemprop \
    --save_path ckpt/AB_combined_chemprop/fragments_to_model_scores.json \
    --model_type chemprop
```

Chemprop with RDKit
```
python map_fragments_to_model_scores.py \
    --fragment_path data/2021q3-4_Enamine_REAL_reagents_SMILES_no_salts.csv \
    --model_path ckpt/AB_combined_chemprop_rdkit \
    --save_path ckpt/AB_combined_chemprop_rdkit/fragments_to_model_scores.json \
    --model_type chemprop \
    --fingerprint_type rdkit
```


## Generate Molecules

Run MCTS with the random forest and chemprop models to generate molecules with one reaction.

Random forest
```
python tree_search.py \
    --model_path ckpt/AB_combined_RF_rdkit \
    --fragment_path data/2021q3-4_Enamine_REAL_reagents_SMILES_no_salts.csv \
    --reagent_to_fragments_path data/reagents_to_fragments.json \
    --fragment_to_model_score_path ckpt/AB_combined_RF_rdkit/fragments_to_model_scores.json \
    --save_dir generations/mcts_AB_combined_RF_rdkit \
    --search_type mcts \
    --model_type rf \
    --fingerprint_type rdkit \
    --n_rollout 5000 \
    --fragment_diversity \
    --max_reactions 1
```

Chemprop
```
python tree_search.py \
    --model_path ckpt/AB_combined_chemprop \
    --fragment_path data/2021q3-4_Enamine_REAL_reagents_SMILES_no_salts.csv \
    --reagent_to_fragments_path data/reagents_to_fragments.json \
    --fragment_to_model_score_path ckpt/AB_combined_chemprop/fragments_to_model_scores.json \
    --save_dir generations/mcts_AB_combined_chemprop \
    --search_type mcts \
    --model_type chemprop \
    --n_rollout 5000 \
    --fragment_diversity \
    --max_reactions 1
```

Chemprop with RDKit
```
python tree_search.py \
    --model_path ckpt/AB_combined_chemprop_rdkit \
    --fragment_path data/2021q3-4_Enamine_REAL_reagents_SMILES_no_salts.csv \
    --reagent_to_fragments_path data/reagents_to_fragments.json \
    --fragment_to_model_score_path ckpt/AB_combined_chemprop_rdkit/fragments_to_model_scores.json \
    --save_dir generations/mcts_AB_combined_chemprop_rdkit \
    --search_type mcts \
    --model_type chemprop \
    --fingerprint_type rdkit \
    --n_rollout 5000 \
    --fragment_diversity \
    --max_reactions 1
```

TODO: random search


## Assess Generated Molecules

Optionally, assess various quantities about the generated molecules.

Random forest
```
python assess_generated_molecules.py \
    --data_path generations/mcts_AB_combined_RF_rdkit/molecules.csv \
    --save_dir generations/mcts_AB_combined_RF_rdkit \
    --train_path data/screening_data/AB_combined.csv
    --train_hits_path data/screening_data/AB_combined_hits.csv
```

Chemprop
```
python assess_generated_molecules.py \
    --data_path generations/mcts_AB_combined_chemprop/molecules.csv \
    --save_dir generations/mcts_AB_combined_chemprop \
    --train_path data/screening_data/AB_combined.csv
    --train_hits_path data/screening_data/AB_combined_hits.csv
```

Chemprop with RDKit
```
python assess_generated_molecules.py \
    --data_path generations/mcts_AB_combined_chemprop_rdkit/molecules.csv \
    --save_dir generations/mcts_AB_combined_chemprop_rdkit \
    --train_path data/screening_data/AB_combined.csv
    --train_hits_path data/screening_data/AB_combined_hits.csv
```


## Select Generated Molecules

Run a series of filtering steps to select molecules for experimental testing


### Train Hit Nearest Neighbor

Compute the nearest neighbor Tversky similarity to train hits using [chem_utils](https://github.com/swansonk14/chem_utils).

Random forest
```
python nearest_neighbor.py \
    --data_path ../../combinatorial_antibiotics/generations/mcts_AB_combined_RF_rdkit/molecules.csv \
    --reference_data_path ../../combinatorial_antibiotics/data/screening_data/AB_combined_hits.csv \
    --reference_name train_hits \
    --metrics tversky
```

Chemprop
```
python nearest_neighbor.py \
    --data_path ../../combinatorial_antibiotics/generations/mcts_AB_combined_chemprop/molecules.csv \
    --reference_data_path ../../combinatorial_antibiotics/data/screening_data/AB_combined_hits.csv \
    --reference_name train_hits \
    --metrics tversky
```

Chemprop with RDKit
```
python nearest_neighbor.py \
    --data_path ../../combinatorial_antibiotics/generations/mcts_AB_combined_chemprop_rdkit/molecules.csv \
    --reference_data_path ../../combinatorial_antibiotics/data/screening_data/AB_combined_hits.csv \
    --reference_name train_hits \
    --metrics tversky
```


### ChEMBL Antibiotic Nearest Neighbor

Compute the nearest neighbor Tversky similarity to ChEMBL antibiotics using [chem_utils](https://github.com/swansonk14/chem_utils).

Random forest
```
python nearest_neighbor.py \
    --data_path ../../combinatorial_antibiotics/generations/mcts_AB_combined_RF_rdkit/molecules.csv \
    --reference_data_path ../../combinatorial_antibiotics/data/chembl/chembl_antibacterial_antibiotic.csv \
    --reference_name chembl_antibacterial_antibiotic \
    --metrics tversky \
    --reference_smiles_column canonical_smiles
```

Chemprop
```
python nearest_neighbor.py \
    --data_path ../../combinatorial_antibiotics/generations/mcts_AB_combined_chemprop/molecules.csv \
    --reference_data_path ../../combinatorial_antibiotics/data/chembl/chembl_antibacterial_antibiotic.csv \
    --reference_name chembl_antibacterial_antibiotic \
    --metrics tversky \
    --reference_smiles_column canonical_smiles
```

Chemprop with RDKit
```
python nearest_neighbor.py \
    --data_path ../../combinatorial_antibiotics/generations/mcts_AB_combined_chemprop_rdkit/molecules.csv \
    --reference_data_path ../../combinatorial_antibiotics/data/chembl/chembl_antibacterial_antibiotic.csv \
    --reference_name chembl_antibacterial_antibiotic \
    --metrics tversky \
    --reference_smiles_column canonical_smiles
```


### Filter by Similarity to Train Hits

Filter to only keep molecules with nearest neighbor Tverksy similarity to the train hits <= 0.5 using [chem_utils](https://github.com/swansonk14/chem_utils).

TODO: check these similarities

Random forest
```
python filter_molecules.py \
    --data_path ../../combinatorial_antibiotics/generations/mcts_AB_combined_RF_rdkit/molecules.csv \
    --save_path ../../combinatorial_antibiotics/generations/mcts_AB_combined_RF_rdkit/molecules_train_sim_below_0.5.csv \
    --filter_column train_hits_tversky_nearest_neighbor_similarity \
    --max_value 0.5
```

Chemprop
```
python filter_molecules.py \
    --data_path ../../combinatorial_antibiotics/generations/mcts_AB_combined_chemprop/molecules.csv \
    --save_path ../../combinatorial_antibiotics/generations/mcts_AB_combined_chemprop/molecules_train_sim_below_0.5.csv \
    --filter_column train_hits_tversky_nearest_neighbor_similarity \
    --max_value 0.5
```

Chemprop with RDKit
```
python filter_molecules.py \
    --data_path ../../combinatorial_antibiotics/generations/mcts_AB_combined_chemprop_rdkit/molecules.csv \
    --save_path ../../combinatorial_antibiotics/generations/mcts_AB_combined_chemprop_rdkit/molecules_train_sim_below_0.5.csv \
    --filter_column train_hits_tversky_nearest_neighbor_similarity \
    --max_value 0.5
```


### Filter by Similarity to ChEMBL Antibiotics

Filter to only keep molecules with nearest neighbor Tverksy similarity to the ChEMBL antibiotics <= 0.5 using [chem_utils](https://github.com/swansonk14/chem_utils).

TODO: check if we want to do this

Random forest
```
python filter_molecules.py \
    --data_path ../../combinatorial_antibiotics/generations/mcts_AB_combined_RF_rdkit/molecules_train_sim_below_0.5.csv \
    --save_path ../../combinatorial_antibiotics/generations/mcts_AB_combined_RF_rdkit/molecules_train_sim_below_0.5_chembl_sim_below_0.5.csv \
    --filter_column chembl_antibacterial_antibiotic_tversky_nearest_neighbor_similarity \
    --max_value 0.5
```

Chemprop
```
python filter_molecules.py \
    --data_path ../../combinatorial_antibiotics/generations/mcts_AB_combined_chemprop/molecules_train_sim_below_0.5.csv \
    --save_path ../../combinatorial_antibiotics/generations/mcts_AB_combined_chemprop/molecules_train_sim_below_0.5_chembl_sim_below_0.5.csv \
    --filter_column chembl_antibacterial_antibiotic_tversky_nearest_neighbor_similarity \
    --max_value 0.5
```

Chemprop with RDKit
```
python filter_molecules.py \
    --data_path ../../combinatorial_antibiotics/generations/mcts_AB_combined_chemprop_rdkit/molecules_train_sim_below_0.5.csv \
    --save_path ../../combinatorial_antibiotics/generations/mcts_AB_combined_chemprop_rdkit/molecules_train_sim_below_0.5_chembl_sim_below_0.5.csv \
    --filter_column chembl_antibacterial_antibiotic_tversky_nearest_neighbor_similarity \
    --max_value 0.5
```


### Filter by Model Score

Filter to only keep molecules with the top 20% of model scores using [chem_utils](https://github.com/swansonk14/chem_utils).

Random forest
```
python filter_molecules.py \
    --data_path ../../combinatorial_antibiotics/generations/mcts_AB_combined_RF_rdkit/molecules_train_sim_below_0.5_chembl_sim_below_0.5.csv \
    --save_path ../../combinatorial_antibiotics/generations/mcts_AB_combined_RF_rdkit/molecules_train_sim_below_0.5_chembl_sim_below_0.5_top_20_percent.csv \
    --filter_column score \
    --top_proportion 0.2
```

Chemprop
```
python filter_molecules.py \
    --data_path ../../combinatorial_antibiotics/generations/mcts_AB_combined_chemprop/molecules_train_sim_below_0.5_chembl_sim_below_0.5.csv \
    --save_path ../../combinatorial_antibiotics/generations/mcts_AB_combined_chemprop/molecules_train_sim_below_0.5_chembl_sim_below_0.5_top_20_percent.csv \
    --filter_column score \
    --top_proportion 0.2
```

Chemprop with RDKit
```
python filter_molecules.py \
    --data_path ../../combinatorial_antibiotics/generations/mcts_AB_combined_chemprop_rdkit/molecules_train_sim_below_0.5_chembl_sim_below_0.5.csv \
    --save_path ../../combinatorial_antibiotics/generations/mcts_AB_combined_chemprop_rdkit/molecules_train_sim_below_0.5_chembl_sim_below_0.5_top_20_percent.csv \
    --filter_column score \
    --top_proportion 0.2
```


### Cluster Molecules

Cluster molecules based on their Morgan fingerprint using [chem_utils](https://github.com/swansonk14/chem_utils).

Random forest
```
python cluster_molecules.py \
    --data_path ../../combinatorial_antibiotics/generations/mcts_AB_combined_RF_rdkit/molecules_train_sim_below_0.5_chembl_sim_below_0.5_top_20_percent.csv \
    --num_clusters 300
```

Chemprop
```
python cluster_molecules.py \
    --data_path ../../combinatorial_antibiotics/generations/mcts_AB_combined_chemprop/molecules_train_sim_below_0.5_chembl_sim_below_0.5_top_20_percent.csv \
    --num_clusters 300
```

Chemprop with RDKit
```
python cluster_molecules.py \
    --data_path ../../combinatorial_antibiotics/generations/mcts_AB_combined_chemprop_rdkit/molecules_train_sim_below_0.5_chembl_sim_below_0.5_top_20_percent.csv \
    --num_clusters 300
```


### Select Molecules from Clusters

Select the top scoring molecule from each cluster using [chem_utils](https://github.com/swansonk14/chem_utils).

Random forest
```
python select_from_clusters.py \
    --data_path ../../combinatorial_antibiotics/generations/mcts_AB_combined_RF_rdkit/molecules_train_sim_below_0.5_chembl_sim_below_0.5_top_20_percent.csv \
    --save_path ../../combinatorial_antibiotics/generations/mcts_AB_combined_RF_rdkit/molecules_train_sim_below_0.5_chembl_sim_below_0.5_top_20_percent_selected_300.csv \
    --value_column score
```

Chemprop
```
python select_from_clusters.py \
    --data_path ../../combinatorial_antibiotics/generations/mcts_AB_combined_chemprop/molecules_train_sim_below_0.5_chembl_sim_below_0.5_top_20_percent.csv \
    --save_path ../../combinatorial_antibiotics/generations/mcts_AB_combined_chemprop/molecules_train_sim_below_0.5_chembl_sim_below_0.5_top_20_percent_selected_300.csv \
    --value_column score
```

Chemprop with RDKit
```
python select_from_clusters.py \
    --data_path ../../combinatorial_antibiotics/generations/mcts_AB_combined_chemprop_rdkit/molecules_train_sim_below_0.5_chembl_sim_below_0.5_top_20_percent.csv \
    --save_path ../../combinatorial_antibiotics/generations/mcts_AB_combined_chemprop_rdkit/molecules_train_sim_below_0.5_chembl_sim_below_0.5_top_20_percent_selected_300.csv \
    --value_column score
```

We now have 100 molecules selected from each model's search that meet our desired train hit similarity and model score criteria.


### Visualize Molecules

Optionally, visualize the selecting molecules using [chem_utils](https://github.com/swansonk14/chem_utils).

Random forest
```
python visualize_molecules.py \
    --data_path ../../combinatorial_antibiotics/generations/mcts_AB_combined_RF_rdkit/molecules_train_sim_below_0.5_chembl_sim_below_0.5_top_20_percent_selected_300.csv \
    --save_dir ../../combinatorial_antibiotics/generations/mcts_AB_combined_RF_rdkit/molecules_train_sim_below_0.5_chembl_sim_below_0.5_top_20_percent_selected_300
```

Chemprop
```
python visualize_molecules.py \
    --data_path ../../combinatorial_antibiotics/generations/mcts_AB_combined_chemprop/molecules_train_sim_below_0.5_chembl_sim_below_0.5_top_20_percent_selected_300.csv \
    --save_dir ../../combinatorial_antibiotics/generations/mcts_AB_combined_chemprop/molecules_train_sim_below_0.5_chembl_sim_below_0.5_top_20_percent_selected_300
```

Chemprop with RDKit
```
python visualize_molecules.py \
    --data_path ../../combinatorial_antibiotics/generations/mcts_AB_combined_chemprop_rdkit/molecules_train_sim_below_0.5_chembl_sim_below_0.5_top_20_percent_selected_300.csv \
    --save_dir ../../combinatorial_antibiotics/generations/mcts_AB_combined_chemprop_rdkit/molecules_train_sim_below_0.5_chembl_sim_below_0.5_top_20_percent_selected_300
```


### Map Molecules to REAL IDs

Map generated molecules to REAL IDs in the format expected by Enamine.

TODO: either remove the reaction counts limitation or add commands to perform the database counting

Random forest
```
python map_generated_molecules_to_real_ids.py \
    --data_path generations/mcts_AB_combined_RF_rdkit/molecules_train_sim_below_0.5_chembl_sim_below_0.5_top_20_percent_selected_300.csv \
    --smiles_save_path generations/mcts_AB_combined_RF_rdkit/molecules_train_sim_below_0.5_chembl_sim_below_0.5_top_20_percent_selected_300_real_ids.csv \
    --sdf_save_path generations/mcts_AB_combined_RF_rdkit/molecules_train_sim_below_0.5_chembl_sim_below_0.5_top_20_percent_selected_300_real_ids.sdf \
    --reaction_counts_path data/Enamine_REAL_counts/reaction_counts.csv \
    --reaction_percent_min 0.001
```

Chemprop
```
python map_generated_molecules_to_real_ids.py \
    --data_path generations/mcts_AB_combined_chemprop/molecules_train_sim_below_0.5_chembl_sim_below_0.5_top_20_percent_selected_300.csv \
    --smiles_save_path generations/mcts_AB_combined_chemprop/molecules_train_sim_below_0.5_chembl_sim_below_0.5_top_20_percent_selected_300_real_ids.csv \
    --sdf_save_path generations/mcts_AB_combined_chemprop/molecules_train_sim_below_0.5_chembl_sim_below_0.5_top_20_percent_selected_300_real_ids.sdf \
    --reaction_counts_path data/Enamine_REAL_counts/reaction_counts.csv \
    --reaction_percent_min 0.001
```

Chemprop with RDKit
```
python map_generated_molecules_to_real_ids.py \
    --data_path generations/mcts_AB_combined_chemprop_rdkit/molecules_train_sim_below_0.5_chembl_sim_below_0.5_top_20_percent_selected_300.csv \
    --smiles_save_path generations/mcts_AB_combined_chemprop_rdkit/molecules_train_sim_below_0.5_chembl_sim_below_0.5_top_20_percent_selected_300_real_ids.csv \
    --sdf_save_path generations/mcts_AB_combined_chemprop_rdkit/molecules_train_sim_below_0.5_chembl_sim_below_0.5_top_20_percent_selected_300_real_ids.sdf \
    --reaction_counts_path data/Enamine_REAL_counts/reaction_counts.csv \
    --reaction_percent_min 0.001
```
