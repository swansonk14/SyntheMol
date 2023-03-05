# Plots

This file contains instructions for generating plots analyzing the data and results.

TODO: redo some names here

- [Data plots](#data-plots)
  * [Data distribution](#data-distribution)
  * [t-SNE of training data and ChEMBL data](#t-sne-of-training-data-and-chembl-data)
- [Property prediction model](#property-prediction-model)
  * [Model performance on training data](#model-performance-on-training-data)
  * [Model generalization](#model-generalization)
- [REAL data](#real-data)
  * [REAL reactions](#real-reactions)
  * [REAL reaction and reactant counts](#real-reaction-and-reactant-counts)
  * [REAL vs train molecular properties](#real-vs-train-molecular-properties)
  * [t-SNE of REAL vs train](#t-sne-of-real-vs-train)
- [Model on REAL Data](#model-on-real-data)
  * [Building block scores](#building-block-scores)
  * [Building block vs full molecule scores](#building-block-vs-full-molecule-scores)
  * [Assess REAL molecules](#assess-real-molecules)
- [MCTS Analysis](#mcts-analysis)
  * [Building block diversity](#building-block-diversity)
  * [Score by rollout](#score-by-rollout)
- [Generated Sets](#generated-sets)
  * [Generate set characteristics](#generate-set-characteristics)
  * [Images of generated molecules](#images-of-generated-molecules)
  * [t-SNE of the filtering steps](#t-sne-of-the-filtering-steps)
  * [t-SNE of final generated sets](#t-sne-of-final-generated-sets)
  * [ClinTox predictions](#clintox-predictions)


## Data plots

### Data distribution

Plot data values for each training set.
```bash
python plots/plot_regression_values.py \
    --data_path data/screening_data/AB_original/AB_2560_normalized.csv \
    --rep1_column Rep_1 \
    --rep2_column Rep_2 \
    --save_dir plots/paper/AB_2560_normalized
```

```bash
python plots/plot_regression_values.py \
    --data_path data/screening_data/AB_original/AB_Mar27_normalized.csv \
    --rep1_column Rep_1 \
    --rep2_column Rep_2 \
    --save_dir plots/paper/AB_Mar27_normalized
```

```bash
python plots/plot_regression_values.py \
    --data_path data/screening_data/AB_original/For_gen_AB_DRH.csv \
    --rep1_column Rep_1 \
    --rep2_column Rep_2 \
    --save_dir plots/paper/For_gen_AB_DRH
```

### t-SNE of training data and ChEMBL data

Plot t-SNE of training data and ChEMBL antibiotics using [chem_utils](https://github.com/swansonk14/chem_utils).
```bash
python dimensionality_reduction.py \
    --data_paths ../../combinatorial_antibiotics/data/chembl/chembl_antibacterial_antibiotic.csv \
    ../../combinatorial_antibiotics/data/screening_data/AB_original/AB_2560_normalized.csv \
    ../../combinatorial_antibiotics/data/screening_data/AB_original/AB_Mar27_normalized.csv \
    ../../combinatorial_antibiotics/data/screening_data/AB_original/For_gen_AB_DRH.csv \
    ../../combinatorial_antibiotics/data/screening_data/AB_2560_hits.csv \
    ../../combinatorial_antibiotics/data/screening_data/AB_Mar27_hits.csv \
    ../../combinatorial_antibiotics/data/screening_data/For_gen_AB_DRH_hits.csv \
    --max_molecules 2000 \
    --colors red orange blue purple orange blue purple \
    --data_names chembl_antibiotic AB_2560_normalized AB_Mar27_normalized For_gen_AB_DRH \
     AB_2560_normalized_hits AB_Mar27_normalized_hits For_gen_AB_DRH_hits  \
    --highlight_data_names AB_2560_normalized_hits AB_Mar27_normalized_hits For_gen_AB_DRH_hits \
    --smiles_columns smiles SMILES SMILES SMILES smiles smiles smiles \
    --save_dir ../../combinatorial_antibiotics/plots/paper/tsne/train_vs_train_hits_vs_chembl
```

## Property prediction model

### Model performance on training data

Plot ROC-AUC and PRC-AUC curves for each model. (Replace model paths and names as needed and curve type with ROC or PRC.)
```bash
python plots/plot_auc.py \
    --data_dir ckpt/AB_combined_RF_rdkit \
    --save_dir plots/paper/auc \
    --model_name "Random Forest" \
    --curve_type ROC
```

### Model generalization

Plot model generalization between training sets.

Separately process each training set.
```bash
#!/bin/bash

for DATA_NAME in AB_2560_normalized AB_Mar27_normalized For_gen_AB_DRH
do
python process_data.py \
    --data_paths data/screening_data/AB_original/${DATA_NAME}.csv \
    --save_path data/screening_data/AB_original/${DATA_NAME}_binarized.csv
done
```

Results of separate data processing.
```
AB_2560_normalized
Data size = 2,371
Mean activity = 0.9759337199493885
Std activity = 0.2673771821337539
Activity threshold of mean - 2 std = 0.4411793556818807
Number of hits = 130
Number of non-hits = 2,241

Full data size = 2,371
Data size after dropping non-conflicting duplicates = 2,371
Data size after dropping conflicting duplicates = 2,371

Final data size = 2,371
Number of hits = 130
Number of non-hits = 2,241

AB_Mar27_normalized
Data size = 5,376
Mean activity = 0.9980530249533108
Std activity = 0.13303517604898704
Activity threshold of mean - 2 std = 0.7319826728553367
Number of hits = 112
Number of non-hits = 5,264

Full data size = 5,376
Data size after dropping non-conflicting duplicates = 5,376
Data size after dropping conflicting duplicates = 5,376

Final data size = 5,376
Number of hits = 112
Number of non-hits = 5,264

For_gen_AB_DRH
Data size = 6,680
Mean activity = 0.9701813713618264
Std activity = 0.17623388953330232
Activity threshold of mean - 2 std = 0.6177135922952217
Number of hits = 294
Number of non-hits = 6,386

Full data size = 6,680
Data size after dropping non-conflicting duplicates = 6,648
Data size after dropping conflicting duplicates = 6,646

Final data size = 6,646
Number of hits = 292
Number of non-hits = 6,354
```

Train models on each training set.
```bash
#!/bin/bash

for DATA_NAME in AB_2560_normalized AB_Mar27_normalized For_gen_AB_DRH
do
python train_model.py \
    --data_path data/screening_data/AB_original/${DATA_NAME}_binarized.csv \
    --save_dir ckpt/${DATA_NAME}_RF_rdkit \
    --model_type rf \
    --fingerprint_type rdkit \
    --num_models 10

python train_model.py \
    --data_path data/screening_data/AB_original/${DATA_NAME}_binarized.csv \
    --save_dir ckpt/${DATA_NAME}_chemprop \
    --model_type chemprop \
    --num_models 10

python train_model.py \
    --data_path data/screening_data/AB_original/${DATA_NAME}_binarized.csv \
    --save_dir ckpt/${DATA_NAME}_chemprop_rdkit \
    --model_type chemprop \
    --fingerprint_type rdkit \
    --num_models 10
done
```

Make predictions on other training sets.
```bash
#!/bin/bash

for DATA1_NAME in AB_2560_normalized AB_Mar27_normalized For_gen_AB_DRH
do
    for DATA2_NAME in AB_2560_normalized AB_Mar27_normalized For_gen_AB_DRH
    do
        if [ "$DATA1_NAME" != "$DATA2_NAME" ]
        then
            echo $DATA1_NAME $DATA2_NAME
            python predict_model.py \
                --data_path data/screening_data/AB_original/${DATA2_NAME}_binarized.csv \
                --model_path ckpt/${DATA1_NAME}_RF_rdkit \
                --save_path ckpt/${DATA1_NAME}_RF_rdkit/${DATA2_NAME}_preds.csv \
                --model_type rf \
                --fingerprint_type rdkit \
                --average_preds

            python compute_auc.py \
                --data_path ckpt/${DATA1_NAME}_RF_rdkit/${DATA2_NAME}_preds.csv \
                --pred_column rf_rdkit_ensemble_preds \
                --true_column activity

            python predict_model.py \
                --data_path data/screening_data/AB_original/${DATA2_NAME}_binarized.csv \
                --model_path ckpt/${DATA1_NAME}_chemprop \
                --save_path ckpt/${DATA1_NAME}_chemprop/${DATA2_NAME}_preds.csv \
                --model_type chemprop \
                --average_preds

            python compute_auc.py \
                --data_path ckpt/${DATA1_NAME}_chemprop/${DATA2_NAME}_preds.csv \
                --pred_column chemprop_ensemble_preds \
                --true_column activity

            python predict_model.py \
                --data_path data/screening_data/AB_original/${DATA2_NAME}_binarized.csv \
                --model_path ckpt/${DATA1_NAME}_chemprop_rdkit \
                --save_path ckpt/${DATA1_NAME}_chemprop_rdkit/${DATA2_NAME}_preds.csv \
                --model_type chemprop \
                --fingerprint_type rdkit \
                --average_preds

            python compute_auc.py \
                --data_path ckpt/${DATA1_NAME}_chemprop_rdkit/${DATA2_NAME}_preds.csv \
                --pred_column chemprop_rdkit_ensemble_preds \
                --true_column activity
        fi
    done
done
```

Plot generalization across training sets as a confusion matrix.

```bash
python plots/plot_model_generalization.py \
    --save_dir plots/paper/model_generalization
```


## REAL data

### REAL reactions

Visualize REAL reactions using [chem_utils](https://github.com/swansonk14/chem_utils).
```bash
python visualize_reactions.py \
    --data_path ../../combinatorial_antibiotics/data/real_reaction_smarts.csv \
    --save_dir ../../combinatorial_antibiotics/plots/paper/real_reactions \
    --name_column reaction_id
```

### REAL reaction and reactant counts

Plot REAL reaction and reactant counts.
```bash
python plots/plot_real_counts.py \
    --reaction_counts_path data/Enamine_REAL_space_counts/real_space_reaction_counts.csv \
    --building_block_counts_path data/Enamine_REAL_space_counts/real_space_building_block_counts.csv \
    --save_dir plots/paper/real_counts
```

### REAL vs train molecular properties

Compute properties of REAL building blocks using [chem_utils](https://github.com/swansonk14/chem_utils).
```bash
python compute_properties.py \
    --data_path ../../combinatorial_antibiotics/data/2021q3-4_Enamine_REAL_reagents_SMILES_no_salts_unique.csv \
    --properties logp mol_weight \
    --save_path ../../combinatorial_antibiotics/data/2021q3-4_Enamine_REAL_reagents_SMILES_no_salts_unique_with_properties.csv
```

Compute properties of REAL molecules using [chem_utils](https://github.com/swansonk14/chem_utils).
```bash
python compute_properties.py \
    --data_path ../../combinatorial_antibiotics/data/Enamine_REAL_space_sampled_25k.csv \
    --properties logp mol_weight \
    --save_path ../../combinatorial_antibiotics/data/Enamine_REAL_space_sampled_25k_with_properties.csv
```

Compute properties of train molecules using [chem_utils](https://github.com/swansonk14/chem_utils).
```bash
python compute_properties.py \
    --data_path ../../combinatorial_antibiotics/data/screening_data/AB_combined.csv \
    --properties logp mol_weight qed sa_score \
    --save_path ../../combinatorial_antibiotics/data/screening_data/AB_combined_with_properties.csv
```

Plot logP distributions using [chem_utils](https://github.com/swansonk14/chem_utils).
```bash
python plots/plot_property_distribution.py \
    --data_paths ../../combinatorial_antibiotics/data/2021q3-4_Enamine_REAL_reagents_SMILES_no_salts_unique_with_properties.csv \
    ../../combinatorial_antibiotics/data/Enamine_REAL_space_sampled_25k_with_properties.csv \
    ../../combinatorial_antibiotics/data/screening_data/AB_combined_with_properties.csv \
    --property_column logp \
    --save_dir ../../combinatorial_antibiotics/plots/paper/properties/logp_train_vs_real \
    --min_value -10 \
    --max_value 10
```

Plot molecular weight distributions using [chem_utils](https://github.com/swansonk14/chem_utils).
```bash
python plots/plot_property_distribution.py \
    --data_paths ../../combinatorial_antibiotics/data/2021q3-4_Enamine_REAL_reagents_SMILES_no_salts_unique_with_properties.csv \
    ../../combinatorial_antibiotics/data/Enamine_REAL_space_sampled_25k_with_properties.csv \
    ../../combinatorial_antibiotics/data/screening_data/AB_combined_with_properties.csv \
    --property_column mol_weight \
    --save_dir ../../combinatorial_antibiotics/plots/paper/properties/mol_weight_train_vs_real \
    --max_value 1000
```


### t-SNE of REAL vs train

Plot t-SNE of training data and REAL space sample using [chem_utils](https://github.com/swansonk14/chem_utils).
```bash
python dimensionality_reduction.py \
    --data_paths ../../combinatorial_antibiotics/data/Enamine_REAL_space_sampled_25k.csv \
    ../../combinatorial_antibiotics/data/2021q3-4_Enamine_REAL_reagents_SMILES_no_salts_unique.csv \
    ../../combinatorial_antibiotics/data/screening_data/AB_combined.csv \
    ../../combinatorial_antibiotics/data/screening_data/AB_combined_hits.csv \
    --max_molecules 7500 1000 1000 500 \
    --data_names REAL_molecules REAL_building_blocks train train_hits \
    --highlight_data_names train_hits \
    --save_dir ../../combinatorial_antibiotics/plots/paper/tsne/train_vs_train_hits_vs_real_vs_real_building_blocks
```


## Model on REAL Data

### Building block scores

Plot building block score distribution for each model.
```bash
python plots/plot_building_block_scores.py \
    --building_block_to_score_path ckpt/AB_combined_RF_rdkit/building_block_to_model_scores.json \
    --title "Random Forest Building Block Score Distribution" \
    --save_dir plots/paper/building_block_scores/rf_building_block_scores
```

```bash
python plots/plot_building_block_scores.py \
    --building_block_to_score_path ckpt/AB_combined_chemprop/building_block_to_model_scores.json \
    --title "Chemprop Building Block Score Distribution" \
    --save_dir plots/paper/building_block_scores/chemprop_building_block_scores
```

```bash
python plots/plot_building_block_scores.py \
    --building_block_to_score_path ckpt/AB_combined_chemprop_rdkit/building_block_to_model_scores.json \
    --title "Chemprop RDKit Building Block Score Distribution" \
    --save_dir plots/paper/building_block_scores/chemprop_rdkit_building_block_scores
```

### Building block vs full molecule scores

Plot building block vs full molecule scores for random sample of REAL molecules.

First, make predictions on a random sample of REAL molecules using each model.
```bash
#!/bin/bash

python predict_model.py \
    --data_path data/Enamine_REAL_space_sampled_25k.csv \
    --model_path ckpt/AB_combined_RF_rdkit \
    --model_type rf \
    --fingerprint_type rdkit \
    --average_preds

python predict_model.py \
    --data_path data/Enamine_REAL_space_sampled_25k.csv \
    --model_path ckpt/AB_combined_chemprop \
    --model_type chemprop \
    --average_preds

python predict_model.py \
    --data_path data/Enamine_REAL_space_sampled_25k.csv \
    --model_path ckpt/AB_combined_chemprop_rdkit \
    --model_type chemprop \
    --fingerprint_type rdkit \
    --average_preds
```

Then, plot the building block vs full molecule scores. (Note: Only 24,276 out of 25,000 molecules have all required building block SMILES.)
```bash
python plots/plot_building_block_vs_molecule_scores.py \
    --data_path data/Enamine_REAL_space_sampled_25k.csv \
    --score_column rf_rdkit_ensemble_preds \
    --building_block_path data/2021q3-4_Enamine_REAL_reagents_SMILES_no_salts.csv \
    --building_block_to_score_path ckpt/AB_combined_RF_rdkit/building_block_to_model_scores.json \
    --title "Random Forest Full Molecule vs Average Building Block Scores" \
    --save_dir plots/paper/full_vs_building_block_scores/rf_rdkit_full_vs_building_block_scores
```

```bash
python plots/plot_building_block_vs_molecule_scores.py \
    --data_path data/Enamine_REAL_space_sampled_25k.csv \
    --score_column chemprop_ensemble_preds \
    --building_block_path data/2021q3-4_Enamine_REAL_reagents_SMILES_no_salts.csv \
    --building_block_to_score_path ckpt/AB_combined_chemprop/building_block_to_model_scores.json \
    --title "Chemprop Full Molecule vs Average Building Block Scores" \
    --save_dir plots/paper/full_vs_building_block_scores/chemprop_full_vs_building_block_scores
```

```bash
python plots/plot_building_block_vs_molecule_scores.py \
    --data_path data/Enamine_REAL_space_sampled_25k.csv \
    --score_column chemprop_rdkit_ensemble_preds \
    --building_block_path data/2021q3-4_Enamine_REAL_reagents_SMILES_no_salts.csv \
    --building_block_to_score_path ckpt/AB_combined_chemprop_rdkit/building_block_to_model_scores.json \
    --title "Chemprop RDKit Full Molecule vs Average Building Block Scores" \
    --save_dir plots/paper/full_vs_building_block_scores/chemprop_rdkit_full_vs_building_block_scores
```

### Assess REAL molecules

Assess scores and similarity distributions of REAL molecules.
```bash
python assess_real_molecules.py \
    --data_path data/Enamine_REAL_space_sampled_25k.csv \
    --save_dir plots/paper/real_analysis \
    --train_hits_path data/screening_data/AB_combined_hits.csv \
    --score_columns rf_rdkit_ensemble_preds chemprop_ensemble_preds chemprop_rdkit_ensemble_preds
```


## MCTS Analysis

### Building block diversity

Building block counts before and after building block diversity. Run `tree_search.py` using same settings as above but without `--fragment_diversity` and run assess_generated_molecules.py and look at `fragment_counts.pdf`.

### Score by rollout

Score of molecules binned by rollout.
```bash
python plots/plot_mcts_over_time.py \
    --data_path generations/mcts_AB_combined_rf_rdkit_ids_20k/molecules.csv \
    --save_dir plots/paper/mcts_over_time/mcts_over_time_rf_rdkit_violin \
    --model_name "Random Forest" \
    --plot_type violin \
    --increment 2000
```

```bash
python plots/plot_mcts_over_time.py \
    --data_path generations/mcts_AB_combined_chemprop_ids_20k/molecules.csv \
    --save_dir plots/paper/mcts_over_time/mcts_over_time_chemprop_violin \
    --model_name "Chemprop" \
    --plot_type violin \
    --increment 2000
```

```bash
python plots/plot_mcts_over_time.py \
    --data_path generations/mcts_AB_combined_chemprop_rdkit_ids_20k/molecules.csv \
    --save_dir plots/paper/mcts_over_time/mcts_over_time_chemprop_rdkit_violin \
    --model_name "Chemprop RDKit" \
    --plot_type violin \
    --increment 2000
```


## Generated Sets

### Generate set characteristics

Assess generated molecules (both original 20k and final 50) for each model. This has already been done using assess_generated_molecules.py above.

### Images of generated molecules

TODO: Images of generated molecules with indications of model, synthesis success, and experimental scores.

### t-SNE of the filtering steps

t-SNE for each step of filtering using [chem_utils](https://github.com/swansonk14/chem_utils).
Note: Need to uncomment hack in dimensionality_reduction.py.
```bash
python dimensionality_reduction.py \
    --data_paths ../../combinatorial_antibiotics/data/chembl/chembl_antibacterial_antibiotic.csv \
    ../../combinatorial_antibiotics/data/screening_data/AB_combined.csv \
    ../../combinatorial_antibiotics/data/screening_data/AB_combined_hits.csv \
    ../../combinatorial_antibiotics/generations/mcts_AB_combined_RF_rdkit_ids_20k/molecules.csv \
    ../../combinatorial_antibiotics/generations/mcts_AB_combined_RF_rdkit_ids_20k/molecules_train_sim_below_0.5.csv \
    ../../combinatorial_antibiotics/generations/mcts_AB_combined_RF_rdkit_ids_20k/molecules_train_sim_below_0.5_chembl_sim_below_0.5.csv \
    ../../combinatorial_antibiotics/generations/mcts_AB_combined_RF_rdkit_ids_20k/molecules_train_sim_below_0.5_chembl_sim_below_0.5_top_20_percent.csv \
    ../../combinatorial_antibiotics/generations/mcts_AB_combined_RF_rdkit_ids_20k/molecules_train_sim_below_0.5_chembl_sim_below_0.5_top_20_percent_selected_50.csv \
    --data_names chembl_antibiotic train train_hits random_forest random_forest_train_sim random_forest_train_sim_chembl_sim \
    random_forest_train_sim_chembl_sim_top_score random_forest_train_sim_chembl_sim_top_score_selected \
    --colors orange blue purple brown brown brown brown red \
    --highlight_data_names random_forest_train_sim_chembl_sim_top_score_selected \
    --save_dir ../../combinatorial_antibiotics/plots/paper/tsne/random_forest_filter
```

```bash
python dimensionality_reduction.py \
    --data_paths ../../combinatorial_antibiotics/data/chembl/chembl_antibacterial_antibiotic.csv \
    ../../combinatorial_antibiotics/data/screening_data/AB_combined.csv \
    ../../combinatorial_antibiotics/data/screening_data/AB_combined_hits.csv \
    ../../combinatorial_antibiotics/generations/mcts_AB_combined_chemprop_ids_20k/molecules.csv \
    ../../combinatorial_antibiotics/generations/mcts_AB_combined_chemprop_ids_20k/molecules_train_sim_below_0.5.csv \
    ../../combinatorial_antibiotics/generations/mcts_AB_combined_chemprop_ids_20k/molecules_train_sim_below_0.5_chembl_sim_below_0.5.csv \
    ../../combinatorial_antibiotics/generations/mcts_AB_combined_chemprop_ids_20k/molecules_train_sim_below_0.5_chembl_sim_below_0.5_top_20_percent.csv \
    ../../combinatorial_antibiotics/generations/mcts_AB_combined_chemprop_ids_20k/molecules_train_sim_below_0.5_chembl_sim_below_0.5_top_20_percent_selected_50.csv \
    --data_names chembl_antibiotic train train_hits chemprop chemprop_train_sim chemprop_train_sim_chembl_sim \
    chemprop_train_sim_chembl_sim_top_score chemprop_train_sim_chembl_sim_top_score_selected \
    --colors orange blue purple brown brown brown brown red \
    --highlight_data_names chemprop_train_sim_chembl_sim_top_score_selected \
    --save_dir ../../combinatorial_antibiotics/plots/paper/tsne/chemprop_filter
```

```bash
python dimensionality_reduction.py \
    --data_paths ../../combinatorial_antibiotics/data/chembl/chembl_antibacterial_antibiotic.csv \
    ../../combinatorial_antibiotics/data/screening_data/AB_combined.csv \
    ../../combinatorial_antibiotics/data/screening_data/AB_combined_hits.csv \
    ../../combinatorial_antibiotics/generations/mcts_AB_combined_chemprop_rdkit_ids_20k/molecules.csv \
    ../../combinatorial_antibiotics/generations/mcts_AB_combined_chemprop_rdkit_ids_20k/molecules_train_sim_below_0.5.csv \
    ../../combinatorial_antibiotics/generations/mcts_AB_combined_chemprop_rdkit_ids_20k/molecules_train_sim_below_0.5_chembl_sim_below_0.5.csv \
    ../../combinatorial_antibiotics/generations/mcts_AB_combined_chemprop_rdkit_ids_20k/molecules_train_sim_below_0.5_chembl_sim_below_0.5_top_20_percent.csv \
    ../../combinatorial_antibiotics/generations/mcts_AB_combined_chemprop_rdkit_ids_20k/molecules_train_sim_below_0.5_chembl_sim_below_0.5_top_20_percent_selected_50.csv \
    --data_names chembl_antibiotic train train_hits chemprop_rdkit chemprop_rdkit_train_sim chemprop_rdkit_train_sim_chembl_sim \
    chemprop_rdkit_train_sim_chembl_sim_top_score chemprop_rdkit_train_sim_chembl_sim_top_score_selected \
    --colors orange blue purple brown brown brown brown red \
    --highlight_data_names chemprop_rdkit_train_sim_chembl_sim_top_score_selected \
    --save_dir ../../combinatorial_antibiotics/plots/paper/tsne/chemprop_rdkit_filter
```

### t-SNE of final generated sets

t-SNE for final generated sets using [chem_utils](https://github.com/swansonk14/chem_utils).
```bash
python dimensionality_reduction.py \
    --data_paths ../../combinatorial_antibiotics/data/screening_data/AB_combined.csv \
    ../../combinatorial_antibiotics/generations/mcts_AB_combined_RF_rdkit_ids_20k/molecules.csv \
    ../../combinatorial_antibiotics/generations/mcts_AB_combined_chemprop_ids_20k/molecules.csv \
    ../../combinatorial_antibiotics/generations/mcts_AB_combined_chemprop_rdkit_ids_20k/molecules.csv \
    ../../combinatorial_antibiotics/data/screening_data/AB_combined_hits.csv \
    ../../combinatorial_antibiotics/generations/mcts_AB_combined_RF_rdkit_ids_20k/molecules_train_sim_below_0.5_chembl_sim_below_0.5_top_20_percent_selected_50.csv \
    ../../combinatorial_antibiotics/generations/mcts_AB_combined_chemprop_ids_20k/molecules_train_sim_below_0.5_chembl_sim_below_0.5_top_20_percent_selected_50.csv \
    ../../combinatorial_antibiotics/generations/mcts_AB_combined_chemprop_rdkit_ids_20k/molecules_train_sim_below_0.5_chembl_sim_below_0.5_top_20_percent_selected_50.csv \
    --data_names train random_forest_full chemprop_full chemprop_rdkit_full train_hits random_forest chemprop chemprop_rdkit \
    --max_molecules 2000 \
    --highlight_data_names train_hits random_forest chemprop chemprop_rdkit \
    --colors blue black red orange blue black red orange \
    --save_dir ../../combinatorial_antibiotics/plots/paper/tsne/train_vs_train_hits_vs_generated_selected
```

### ClinTox predictions

Show potent molecules vs test molecules.
```bash
python plots/plot_toxicity.py \
    --test_dir paper/models/clintox_chemprop_rdkit \
    --generated_path paper/data/8_synthesized/potent.csv \
    --save_dir paper/plots/toxicity
```

```
Found 10 test set predictions
Size of test set predictions = 1,478
Size of generated predictions = 6
Generated preds = [0.054, 0.078, 0.123, 0.200, 0.240, 0.247]
Best F1 threshold = 0.416
F1 = 0.475
Precision = 0.547
Recall = 0.420
Percentiles of generated preds among toxic molecules = [14.3, 21.4, 31.2, 45.5, 47.3, 47.3]
Percentiles of generated preds among non-toxic molecules = [76.2, 80.1, 86.4, 91.7, 93.1, 93.2]
```
