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

Determine which reactions are most common in REAL space.

```bash
python count_real_space.py \
    --data_dir data/Enamine_REAL_space \
    --save_path data/reaction_counts_REAL_space.csv \
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

All 138,085 molecules were successfully converted from SDF to SMILES, and among those 134,609 are unique.

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

This removes 25 molecules whose salts cannot be stripped, leaving 138,060 molecules.

Note: This step is crucial to prevent errors in running reactions. Salts can cause reactions to create products that are the same as the reactants, leading to undesired infinite loops during molecule generation.


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
    --save_dir ckpt/AB_combined_${MODEL} \
    --model_type rf \
    --fingerprint_type rdkit \
    --num_models 10
```

Chemprop
```bash
python train_model.py \
    --data_path data/screening_data/AB_combined.csv \
    --save_dir ckpt/AB_combined_chemprop \
    --model_type chemprop \
    --num_models 10
```

Chemprop with RDKit
```bash
python train_model.py \
    --data_path data/screening_data/AB_combined.csv \
    --save_dir ckpt/AB_combined_chemprop_rdkit \
    --model_type chemprop \
    --fingerprint_type rdkit \
    --num_models 10
```

| Model          | ROC-AUC         | PRC-AUC         |
|----------------|-----------------|-----------------|
| RF RDKit       | 0.825 +/- 0.035 | 0.401 +/- 0.099 |
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
    --train_path data/screening_data/AB_combined.csv \
    --train_hits_path data/screening_data/AB_combined_hits.csv
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
    --metrics tversky \
    --reference_smiles_column canonical_smiles
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
    --save_dir ../../combinatorial_antibiotics/generations/${NAME}/molecules_train_sim_below_0.5_chembl_sim_below_0.5_top_20_percent_selected_50
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

## Train Toxicity Model

Perform these steps from the chemprop repo.

Generate RDKit features for toxicity data.
```bash
python scripts/save_features.py \
    --data_path data/clintox.csv \
    --save_path features/clintox.npz
```

Train chemprop + RDKit model on toxicity data.
```bash
python train.py \
    --data_path data/clintox.csv \
    --dataset_type classification \
    --metric prc-auc \
    --extra_metrics auc \
    --num_folds 10 \
    --features_path features/clintox.npz \
    --no_features_scaling \
    --save_dir ckpt/clintox \
    --quiet
```

```
Overall test prc-auc = 0.707140 +/- 0.062722
Overall test auc = 0.882426 +/- 0.040230
```

Generate RDKit features for synthesized molecules
```bash
python scripts/save_features.py \
    --data_path ../combinatorial_antibiotics/generations/enamine_synthesized_results_formatted.csv \
    --save_path ../combinatorial_antibiotics/generations/enamine_synthesized_results_formatted.npz \
    --smiles_column smiles
```

Make predictions on synthesized molecules.
```bash
python predict.py \
    --test_path ../combinatorial_antibiotics/generations/enamine_synthesized_results_formatted.csv \
    --preds_path ../combinatorial_antibiotics/generations/enamine_synthesized_results_formatted.csv \
    --checkpoint_dir ckpt/clintox \
    --features_path ../combinatorial_antibiotics/generations/enamine_synthesized_results_formatted.npz \
    --no_features_scaling \
    --smiles_column smiles
```

## Plots

### Data Plots

Plot data values for each training set.
```bash
python plot_regression_values.py \
    --data_path data/screening_data/AB_original/AB_2560_normalized.csv \
    --rep1_column Rep_1 \
    --rep2_column Rep_2 \
    --save_dir plots/paper/AB_2560_normalized
```

```bash
python plot_regression_values.py \
    --data_path data/screening_data/AB_original/AB_Mar27_normalized.csv \
    --rep1_column Rep_1 \
    --rep2_column Rep_2 \
    --save_dir plots/paper/AB_Mar27_normalized
```

```bash
python plot_regression_values.py \
    --data_path data/screening_data/AB_original/For_gen_AB_DRH.csv \
    --rep1_column Rep_1 \
    --rep2_column Rep_2 \
    --save_dir plots/paper/For_gen_AB_DRH
```

Plot t-SNE of training data and ChEMBL antibiotics using [chem_utils](https://github.com/swansonk14/chem_utils).
```bash
python dimensionality_reduction.py \
    --data_paths ../../combinatorial_antibiotics/data/screening_data/AB_original/AB_2560_normalized.csv \
    ../../combinatorial_antibiotics/data/screening_data/AB_original/AB_Mar27_normalized.csv \
    ../../combinatorial_antibiotics/data/screening_data/AB_original/For_gen_AB_DRH.csv \
    ../../combinatorial_antibiotics/data/chembl/chembl_antibacterial_antibiotic.csv \
    ../../combinatorial_antibiotics/data/screening_data/AB_2560_hits.csv \
    ../../combinatorial_antibiotics/data/screening_data/AB_Mar27_hits.csv \
    ../../combinatorial_antibiotics/data/screening_data/For_gen_AB_DRH_hits.csv \
    --max_molecules 2000 \
    --data_names AB_2560_normalized AB_Mar27_normalized For_gen_AB_DRH chembl_antibiotic \
     AB_2560_normalized_hits AB_Mar27_normalized_hits For_gen_AB_DRH_hits  \
    --highlight_data_names AB_2560_normalized_hits AB_Mar27_normalized_hits For_gen_AB_DRH_hits \
    --smiles_columns SMILES SMILES SMILES canonical_smiles smiles smiles smiles \
    --save_path ../../combinatorial_antibiotics/plots/paper/tsne/train_vs_train_hits_vs_chembl.pdf
```

### REAL plots

Visualize REAL reactions using [chem_utils](https://github.com/swansonk14/chem_utils)..
```bash
python visualize_reactions.py \
    --data_path ../../combinatorial_antibiotics/data/real_reaction_smarts.csv \
    --save_dir ../../combinatorial_antibiotics/plots/paper/real_reactions \
    --name_column reaction_id
```

Plot REAL reaction and reactant counts.
TODO: count reactants and plot
```bash
python plot_real_counts.py \
    --reaction_counts_path data/reaction_counts_REAL_space.csv \
    --save_dir plots/paper/real_counts
```

Plot molecular weight distribution of REAL and train molecules.
TODO: need actual sample of REAL space, not REAL database
```bash
python property_distribution.py \
    --data_paths ../../combinatorial_antibiotics/data/Enamine_REAL_SMILES_sampled.csv \
    ../../combinatorial_antibiotics/data/screening_data/AB_combined.csv \
    --property mol_weight \
    --max_value 1000 \
    --save_path ../../combinatorial_antibiotics/plots/paper/mol_weight_train_vs_real.pdf
```

Plot logP distribution of REAL and train molecules.
TODO: need actual sample of REAL space, not REAL database
```bash
python property_distribution.py \
    --data_paths ../../combinatorial_antibiotics/data/Enamine_REAL_SMILES_sampled.csv \
    ../../combinatorial_antibiotics/data/screening_data/AB_combined.csv \
    --property logp \
    --min_value -10 \
    --max_value 10 \
    --save_path ../../combinatorial_antibiotics/plots/paper/logp_train_vs_real.pdf
```

TODO: coverage of train and train hits by fragments
Do Tversky similarity between training set and fragments with fragment as reference but look for greatest similarity for each fragment across training molecules to show the greatest degree to which each fragment appears in a training molecule

Plot t-SNE of training data and REAL space sample using [chem_utils](https://github.com/swansonk14/chem_utils).
TODO: need actual sample of REAL space, not REAL database.
```bash
python dimensionality_reduction.py \
    --data_paths ../../combinatorial_antibiotics/data/Enamine_REAL_SMILES_sampled.csv \
    ../../combinatorial_antibiotics/data/screening_data/AB_combined.csv \
    ../../combinatorial_antibiotics/data/screening_data/AB_combined_hits.csv \
    --max_molecules 2000 \
    --data_names Enamine_REAL train train_hits \
    --highlight_data_names train_hits \
    --save_path ../../combinatorial_antibiotics/plots/paper/tsne/train_vs_train_hits_vs_real.pdf
```

### Model on Training Data

TODO: AUC and PRC-AUC curves for random and scaffold splits

Plot model generalization between training sets.
TODO: update/fill in numbers and implement confusion matrix
```bash
python plot_model_generalization.py \
    --save_dir plots/paper/model_generalization
```


### Model on REAL Data

Plot fragment score distribution for each model.
```bash
python plot_fragment_scores.py \
    --fragment_to_score_path ckpt/AB_combined_RF_rdkit/fragments_to_model_scores.json \
    --title "Random Forest Fragment Score Distribution" \
    --save_path plots/paper/fragment_scores/rf_fragment_scores.pdf
```

```bash
python plot_fragment_scores.py \
    --fragment_to_score_path ckpt/AB_combined_chemprop/fragments_to_model_scores.json \
    --title "Chemprop Fragment Score Distribution" \
    --save_path plots/paper/fragment_scores/chemprop_fragment_scores.pdf
```

```bash
python plot_fragment_scores.py \
    --fragment_to_score_path ckpt/AB_combined_chemprop_rdkit/fragments_to_model_scores.json \
    --title "Chemprop RDKit Fragment Score Distribution" \
    --save_path plots/paper/fragment_scores/chemprop_rdkit_fragment_scores.pdf
```

Plot fragment vs full molecule scores for random sample of REAL molecules.
TODO: maybe use actual random sample rather than randomly generated sample?
```bash
python plot_fragment_vs_molecule_scores.py \
    --data_path generations/random_ids_20k/molecules.csv \
    --score_column rf_rdkit_ensemble_preds \
    --fragment_to_score_path ckpt/AB_combined_RF_rdkit/fragments_to_model_scores.json \
    --title "Random Forest Full Molecule vs Average Fragment Scores" \
    --save_path plots/paper/full_vs_fragment_scores/rf_rdkit_full_vs_fragment_scores.pdf
```

```bash
python plot_fragment_vs_molecule_scores.py \
    --data_path generations/random_ids_20k/molecules.csv \
    --score_column chemprop_ensemble_preds \
    --fragment_to_score_path ckpt/AB_combined_chemprop/fragments_to_model_scores.json \
    --title "Chemprop Full Molecule vs Average Fragment Scores" \
    --save_path plots/paper/full_vs_fragment_scores/chemprop_full_vs_fragment_scores.pdf
```

```bash
python plot_fragment_vs_molecule_scores.py \
    --data_path generations/random_ids_20k/molecules.csv \
    --score_column chemprop_rdkit_ensemble_preds \
    --fragment_to_score_path ckpt/AB_combined_chemprop_rdkit/fragments_to_model_scores.json \
    --title "Chemprop RDKit Full Molecule vs Average Fragment Scores" \
    --save_path plots/paper/full_vs_fragment_scores/chemprop_rdkit_full_vs_fragment_scores.pdf
```

### MCTS Analysis

TODO: Fragment counts before and after fragment diversity (need to do a 20k run without fragment diversity)

Score of molecules binned by rollout.
TODO: make this look nicer
```bash
python plot_mcts_over_time.py \
    --data_path generations/mcts_AB_combined_rf_rdkit_ids_20k/molecules.csv \
    --save_path plots/paper/mcts_over_time/mcts_over_time_rf_rdkit_line.pdf \
    --model_name "Random Forest" \
    --plot_type line \
    --increment 500
```

```bash
python plot_mcts_over_time.py \
    --data_path generations/mcts_AB_combined_chemprop_ids_20k/molecules.csv \
    --save_path plots/paper/mcts_over_time/mcts_over_time_chemprop_line.pdf \
    --model_name "Chemprop" \
    --plot_type line \
    --increment 500
```

```bash
python plot_mcts_over_time.py \
    --data_path generations/mcts_AB_combined_chemprop_rdkit_ids_20k/molecules.csv \
    --save_path plots/paper/mcts_over_time/mcts_over_time_chemprop_rdkit_line.pdf \
    --model_name "Chemprop RDKit" \
    --plot_type line \
    --increment 500
```

TODO: show how often MCTS finds molecules with full molecule score higher than fragment score and see if this is more often than random chance
But maybe not since the numbers don't support this

### Generated Sets