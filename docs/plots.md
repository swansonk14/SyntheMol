# Plots

Instructions for producing plots analyzing the data and results. Assumes relevant data has already been downloaded (see [docs/README.md](README.md)).

- [Data](#data)
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
  * [Chemprop vs SyntheMol](#chemprop-vs-synthemol)
- [MCTS Analysis](#mcts-analysis)
  * [Building block diversity](#building-block-diversity)
  * [Score by rollout](#score-by-rollout)
- [Generated Sets](#generated-sets)
  * [Generate set characteristics](#generate-set-characteristics)
  * [Images of generated molecules](#images-of-generated-molecules)
  * [t-SNE of final generated sets](#t-sne-of-final-generated-sets)
  * [ClinTox predictions](#clintox-predictions)


## Data

### Data distribution

Plot data values for each training set.
```bash
for LIBRARY in library_1 library_2 library_3
do
python scripts/plot/plot_regression_values.py \
    --data_path data/Data/1_training_data/${LIBRARY}.csv \
    --rep1_column Rep_1 \
    --rep2_column Rep_2 \
    --save_dir plots/${LIBRARY}
done
```

### t-SNE of training data and ChEMBL data

Plot t-SNE of training data and ChEMBL antibiotics.
```bash
chemfunc dimensionality_reduction \
    --data_paths data/Data/2_chembl/chembl_antibacterial_antibiotic.csv \
    data/Data/1_training_data/library_1.csv \
    data/Data/1_training_data/library_2.csv \
    data/Data/1_training_data/library_3.csv \
    data/Data/1_training_data/library_1_hits.csv \
    data/Data/1_training_data/library_2_hits.csv \
    data/Data/1_training_data/library_3_hits.csv \
    --max_molecules 2000 \
    --colors red orange blue purple orange blue purple \
    --data_names chembl_antibiotic library_1 library_2 library_3 \
     library_1_hits library_2_hits library_3_hits  \
    --highlight_data_names library_1_hits library_2_hits library_3_hits \
    --smiles_columns smiles SMILES SMILES SMILES smiles smiles smiles \
    --save_path plots/tsne/train_vs_train_hits_vs_chembl.pdf
```

## Property prediction model

### Model performance on training data

Plot ROC-AUC and PRC-AUC curves for each model. (Replace model paths and names as needed and curve type with ROC or PRC.)
```bash
python scripts/plot/plot_auc.py \
    --data_dir data/Models/antibiotic_random_forest \
    --save_dir plots/auc \
    --model_name "Random Forest" \
    --curve_type ROC
```

### Model generalization

Plot model generalization between training sets.

Separately process each training set.
```bash
for LIBRARY in library_1 library_2 library_3
do
python scripts/data/process_data.py \
    --data_paths data/Data/1_training_data/${LIBRARY}.csv \
    --save_path data/Data/1_training_data/${LIBRARY}_binarized.csv
done
```

Results of separate data processing.
```
library_1
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

library_2
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

library_3
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
```

Train models on each training set.
```bash
for LIBRARY in library_1 library_2 library_3
do
python train_model.py \
    --data_path data/Data/1_training_data/${LIBRARY}_binarized.csv \
    --save_dir data/Models/${LIBRARY}_random_forest \
    --model_type random_forest \
    --fingerprint_type rdkit \
    --num_models 10

python train_model.py \
    --data_path data/Data/1_training_data/${LIBRARY}_binarized.csv \
    --save_dir data/Models/${LIBRARY}_chemprop \
    --model_type chemprop \
    --num_modeldata/D 10

python train_model.py \
    --data_path data/Data/1_training_data/${LIBRARY}_binarized.csv \
    --save_dir data/Models/${LIBRARY}_chemprop_rdkit \
    --model_type chemprop \
    --fingerprint_type rdkit \
    --num_models 10
done
```

Make predictions on other training sets.
```bash
for LIBRARY_A in library_1 library_2 library_3
do
    for LIBRARY_B in library_1 library_2 library_3
    do
        if [ "$LIBRARY_A" != "$LIBRARY_B" ]
        then
            echo $LIBRARY_A $LIBRARY_B
            python scripts/models/predict.py \
                --data_path data/Data/1_training_data/${LIBRARY_B}_binarized.csv \
                --model_path data/Models/${LIBRARY_A}_random_forest \
                --save_path data/Models/${LIBRARY_A}_random_forest/${LIBRARY_B}_preds.csv \
                --model_type random_forest \
                --fingerprint_type rdkit \
                --average_preds

            python scripts/models/compute_auc.py \
                --data_path data/Models/${LIBRARY_A}_random_forest/${LIBRARY_B}_preds.csv \
                --pred_column random_forest_ensemble_preds \
                --true_column activity

            python scripts/models/predict.py \
                --data_path data/Data/1_training_data/${LIBRARY_B}_binarized.csv \
                --model_path data/Models/${LIBRARY_A}_chemprop \
                --save_path data/Models/${LIBRARY_A}_chemprop/${LIBRARY_B}_preds.csv \
                --model_type chemprop \
                --average_preds

            python scripts/models/compute_auc.py \
                --data_path data/Models/${LIBRARY_A}_chemprop/${LIBRARY_B}_preds.csv \
                --pred_column chemprop_ensemble_preds \
                --true_column activity

            python scripts/models/predict.py \
                --data_path data/Data/1_training_data/${LIBRARY_B}_binarized.csv \
                --model_path data/Models/${LIBRARY_A}_chemprop_rdkit \
                --save_path data/Models/${LIBRARY_A}_chemprop_rdkit/${LIBRARY_B}_preds.csv \
                --model_type chemprop \
                --fingerprint_type rdkit \
                --average_preds

            python scripts/models/compute_auc.py \
                --data_path data/Models/${LIBRARY_A}_chemprop_rdkit/${DATA2_NAME}_preds.csv \
                --pred_column chemprop_rdkit_ensemble_preds \
                --true_column activity
        fi
    done
done
```

Plot generalization across training sets as a confusion matrix.

```bash
python scripts/plot/plot_model_generalization.py \
    --save_dir plots/model_generalization
```


## REAL data

### REAL reactions

Visualize REAL reactions.
```bash
chemfunc visualize_reactions \
    --data_path data/Data/4_real_space/reactions.csv \
    --save_dir plots/real_reactions \
    --smarts_column reaction_smarts \
    --name_column reaction_id
```

### REAL reaction and reactant counts

Plot REAL reaction and reactant counts.
```bash
python scripts/plot/plot_real_counts.py \
    --reaction_counts_path data/Data/4_real_space/reaction_counts.csv \
    --building_block_counts_path data/Data/4_real_space/building_block_counts.csv \
    --save_dir plots/real_counts
```

### REAL vs train molecular properties

Deduplicate building blocks by SMILES.
```bash
chemfunc deduplicate_smiles \
    --data_path data/Data/4_real_space/building_blocks.csv \
    --save_path data/Data/4_real_space/building_blocks_unique.csv
```

This leaves 132,479 unique building block molecules.

Compute properties of REAL building blocks, REAL molecules, and train molecules.
```bash
for DATA in 4_real_space/building_blocks_unique 4_real_space/random_real 1_training_data/antibiotics
do
chemfunc compute_properties \
    --data_path data/Data/${DATA}.csv \
    --properties logp mol_weight
done
```

Plot logP distributions.
```bash
chemfunc plot_property_distribution \
    --data_paths data/Data/4_real_space/building_blocks_unique.csv \
    data/Data/4_real_space/random_real.csv \
    data/Data/1_training_data/antibiotics.csv \
    --property_column logp \
    --save_dir plots/properties/logp_train_vs_real \
    --min_value -10 \
    --max_value 10
```

Plot molecular weight distributions.
```bash
chemfunc plot_property_distribution \
    --data_paths data/Data/4_real_space/building_blocks_unique.csv \
    data/Data/4_real_space/random_real.csv \
    data/Data/1_training_data/antibiotics.csv \
    --property_column mol_weight \
    --save_dir plots/properties/mol_weight_train_vs_real \
    --max_value 1000
```


### t-SNE of REAL vs train

Plot t-SNE of training data and REAL space sample.
```bash
chemfunc dimensionality_reduction \
    --data_paths data/Data/4_real_space/random_real.csv \
    data/Data/4_real_space/building_blocks_unique.csv \
    data/Data/1_training_data/antibiotics.csv \
    data/Data/1_training_data/antibiotics_hits.csv \
    --max_molecules 7500 1000 1000 500 \
    --data_names REAL_molecules REAL_building_blocks train train_hits \
    --highlight_data_names train_hits \
    --save_path plots/tsne/train_vs_train_hits_vs_real_vs_real_building_blocks.pdf
```


## Model on REAL Data

### Building block scores

Plot building block score distribution for each model.
```bash
python scripts/plot/plot_building_block_scores.py \
    --building_blocks_path data/Models/antibiotic_random_forest/building_block_scores.csv \
    --title "Random Forest Building Block Score Distribution" \
    --save_dir plots/building_block_scores/random_forest_building_block_scores
```

```bash
python scripts/plot/plot_building_block_scores.py \
    --building_blocks_path data/Models/antibiotic_chemprop/building_block_scores.csv \
    --title "Chemprop Building Block Score Distribution" \
    --save_dir plots/building_block_scores/chemprop_building_block_scores
```

```bash
python scripts/plot/plot_building_block_scores.py \
    --building_blocks_path data/Models/antibiotic_chemprop_rdkit/building_block_scores.csv \
    --title "Chemprop RDKit Building Block Score Distribution" \
    --save_dir plots/building_block_scores/chemprop_rdkit_building_block_scores
```

### Building block vs full molecule scores

Plot building block vs full molecule scores for random sample of REAL molecules.

First, make predictions on a random sample of REAL molecules using each model.
```bash
python predict_model.py \
    --data_path data/Data/4_real_space/random_real.csv \
    --model_path data/Models/antibiotic_random_forest \
    --model_type random_forest \
    --fingerprint_type rdkit \
    --average_preds

python predict_model.py \
    --data_path data/Data/4_real_space/random_real.csv \
    --model_path data/Models/antibiotic_chemprop \
    --model_type chemprop \
    --average_preds

python predict_model.py \
    --data_path data/Data/4_real_space/random_real.csv \
    --model_path data/Models/antibiotic_chemprop_rdkit \
    --model_type chemprop \
    --fingerprint_type rdkit \
    --average_preds
```

Then, plot the building block vs full molecule scores. (Note: Only 24,276 out of 25,000 molecules have all required building block SMILES.)
```bash
python scripts/plot/plot_building_block_vs_molecule_scores.py \
    --data_path data/Data/4_real_space/random_real.csv \
    --building_blocks_path data/Models/random_forest/building_block_scores.csv \
    --building_blocks_score_column random_forest_rdkit_ensemble_preds \
    --title "Random Forest Full Molecule vs Average Building Block Scores" \
    --save_dir plots/full_vs_building_block_scores/random_forest_full_vs_building_block_scores
```

```bash
python scripts/plot/plot_building_block_vs_molecule_scores.py \
    --data_path data/Data/4_real_space/random_real.csv \
    --building_blocks_path data/Models/chemprop/building_block_scores.csv \
    --building_blocks_score_column chemprop_ensemble_preds \
    --title "Chemprop Full Molecule vs Average Building Block Scores" \
    --save_dir plots/full_vs_building_block_scores/chemprop_full_vs_building_block_scores
```

```bash
python scripts/plots/plot_building_block_vs_molecule_scores.py \
    --data_path data/Data/4_real_space/random_real.csv \
    --building_blocks_path data/Models/chemprop_rdkit/building_block_scores.csv \
    --building_blocks_score_column chemprop_rdkit_ensemble_preds \
    --title "Chemprop RDKit Full Molecule vs Average Building Block Scores" \
    --save_dir plots/full_vs_building_block_scores/chemprop_rdkit_full_vs_building_block_scores
```

### Assess REAL molecules

Assess scores and similarity distributions of REAL molecules.
```bash
python scripts/plot/plot_molecule_analysis.py \
    --data_path data/Data/4_real_space/random_real.csv \
    --save_dir plots/real_analysis \
    --train_hits_path data/Data/1_training_data/antibiotics_hits.csv \
    --score_columns random_forest_rdkit_ensemble_preds chemprop_ensemble_preds chemprop_rdkit_ensemble_preds
```


### Chemprop vs SyntheMol

Make predictions on a random sample of 10 million REAL molecules using the Chemprop model for a time-based comparison against SyntheMol with Chemprop. Note that is uses a GPU and 16 parallel data loaders.

First, unzip the random sample of 10 million REAL molecules if necessary.
```bash
gunzip data/Data/4_real_space/random_real_10m.csv.gz
```

Now, make Chemprop predictions.
```bash
python scripts/models/predict.py \
    --data_path data/Data/4_real_space/random_real_10m.csv \
    --model_path data/Models/clogp_chemprop_30_epochs \
    --model_type chemprop \
    --average_preds \
    --num_workers 16 \
    --no_cache \
    --use_gpu \
    --save_path data/Data/4_real_space/random_real_10m.csv
```


## MCTS Analysis

### Building block diversity

Building block counts before and after building block diversity. Run `SyntheMol.generate` as before but with the `--no_building_block_diversity` flag. Then run `SyntheMol.plot.plot_generated_molecule_analysis` and look at `building_block_counts.pdf`.

### Score by rollout

Score of molecules binned by rollout.
```bash
python scripts/plot/plot_mcts_over_time.py \
    --data_path data/Data/6_generations_chemprop/molecules.csv \
    --save_dir plots/mcts_over_time/chemprop \
    --model_name "Chemprop" \
    --increment 2000
```

```bash
python scripts/plot/plot_mcts_over_time.py \
    --data_path data/Data/7_generations_chemprop_rdkit/molecules.csv \
    --save_dir plots/mcts_over_time/chemprop_rdkit \
    --model_name "Chemprop RDKit" \
    --increment 2000
```

```bash
python scripts/plot/plot_mcts_over_time.py \
    --data_path data/Data/8_generations_random_forest/molecules.csv \
    --save_dir plots/mcts_over_time/random_forest \
    --model_name "Random Forest" \
    --increment 2000
```


## Generated Sets

### Generate set characteristics

Assess generated molecules for novelty, score, and diversity for each model.
```bash
for NAME in 6_generations_chemprop 7_generations_chemprop_rdkit 8_generations_random_forest
do
python scripts/plot/plot_generated_molecule_analysis.py \
    --data_path data/Data/${NAME}/molecules.csv \
    --save_dir data/Data/${NAME} \
    --reference_paths data/Data/1_training_data/antibiotics_hits.csv data/2_chembl/chembl.csv
done
```

Assess selected molecules for novelty, score, and diversity for each model.
```bash
for NAME in chemprop chemprop_rdkit random_forest
do
python scripts/plot/plot_generated_molecule_analysis.py \
    --data_path data/Data/${NAME}/molecules_train_sim_below_0.5_chembl_sim_below_0.5_top_20_percent_selected_50.csv \
    --save_dir data/Data/${NAME}/analysis_molecules_train_sim_below_0.5_chembl_sim_below_0.5_top_20_percent_selected_50 \
    --reference_paths data/Data/1_training_data/antibiotics_hits.csv data/2_chembl/chembl.csv
done
```

### Images of generated molecules

Visualize the selected molecules.
```bash
for NAME in chemprop chemprop_rdkit random_forest
do
chemfunc visualize_molecules \
    --data_path data/Data/${NAME}/molecules_train_sim_below_0.5_chembl_sim_below_0.5_top_20_percent_selected_50.csv \
    --save_dir data/Data/${NAME}
done
```


### t-SNE of final generated sets

t-SNE for final generated sets.
```bash
chemfunc dimensionality_reduction \
    --data_paths data/Data/1_training_data/antibiotics.csv \
    data/Data/6_generations_chemprop/molecules.csv \
    data/Data/7_generations_chemprop_rdkit/molecules.csv \
    data/Data/8_generations_random_forest/molecules.csv \
    data/Data/1_training_data/antibiotics_hits.csv \
    data/Data/6_generations_chemprop/molecules_train_sim_below_0.5_chembl_sim_below_0.5_top_20_percent_selected_50.csv \
    data/Data/7_generations_chemprop_rdkit/molecules_train_sim_below_0.5_chembl_sim_below_0.5_top_20_percent_selected_50.csv \
    data/Data/8_generations_random_forest/molecules_train_sim_below_0.5_chembl_sim_below_0.5_top_20_percent_selected_50.csv \
    --data_names train random_forest_full chemprop_full chemprop_rdkit_full train_hits random_forest chemprop chemprop_rdkit \
    --max_molecules 2000 \
    --highlight_data_names train_hits chemprop chemprop_rdkit random_forest \
    --colors blue black red orange blue black red orange \
    --save_path plots/tsne/train_vs_train_hits_vs_generated_selected.pdf
```

### ClinTox predictions

Show potent molecules vs test molecules.
```bash
python scripts/plot/plot_toxicity.py \
    --test_dir data/Models/clintox_chemprop_rdkit \
    --generated_path data/Data/8_synthesized/potent.csv \
    --save_dir plots/toxicity
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
