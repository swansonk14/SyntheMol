# Ablation Experiments

Instructions for reproducing the ablation experiments in the paper to determine the importance of various components of
SyntheMol-RL.

Note: All SyntheMol experiments should be run with five random seeds (0-4) using the `--rng_seed <seed>` flag to measure variability.

## Property predictors - scaffold split

### Train models

For each model type, train 10 models using 10-fold cross-validation.

Chemprop for _S. aureus_

```bash
chemprop_train \
    --data_path rl/data/s_aureus/s_aureus.csv \
    --dataset_type classification \
    --target_column s_aureus_activity \
    --num_folds 10 \
    --split_type scaffold_balanced \
    --metric prc-auc \
    --extra_metrics auc \
    --save_dir rl/models/s_aureus_chemprop_scaffold \
    --save_preds \
    --quiet
```

Chemprop for solubility

```bash
chemprop_train \
    --data_path rl/data/solubility/solubility.csv \
    --dataset_type regression \
    --target_column solubility \
    --num_folds 10 \
    --split_type scaffold_balanced \
    --metric mae \
    --extra_metrics r2 \
    --save_dir rl/models/solubility_chemprop_scaffold \
    --save_preds \
    --quiet
```

Chemprop-RDKit for _S. aureus_

```bash
chemprop_train \
    --data_path rl/data/s_aureus/s_aureus.csv \
    --dataset_type classification \
    --target_column s_aureus_activity \
    --features_path rl/data/s_aureus/s_aureus.npz \
    --no_features_scaling \
    --num_folds 10 \
    --split_type scaffold_balanced \
    --metric prc-auc \
    --extra_metrics auc \
    --save_dir rl/models/s_aureus_chemprop_rdkit_scaffold \
    --save_preds \
    --quiet
```

Chemprop-RDKit for solubility

```bash
chemprop_train \
    --data_path rl/data/solubility/solubility.csv \
    --dataset_type regression \
    --target_column solubility \
    --features_path rl/data/solubility/solubility.npz \
    --no_features_scaling \
    --num_folds 10 \
    --split_type scaffold_balanced \
    --metric mae \
    --extra_metrics r2 \
    --save_dir rl/models/solubility_chemprop_rdkit_scaffold \
    --save_preds \
    --quiet
```

MLP-RDKit for _S. aureus_

```bash
chemprop_train \
    --data_path rl/data/s_aureus/s_aureus.csv \
    --dataset_type classification \
    --target_column s_aureus_activity \
    --features_path rl/data/s_aureus/s_aureus.npz \
    --no_features_scaling \
    --features_only \
    --num_folds 10 \
    --split_type scaffold_balanced \
    --metric prc-auc \
    --extra_metrics auc \
    --save_dir rl/models/s_aureus_mlp_rdkit_scaffold \
    --save_preds \
    --quiet
```

MLP-RDKit for solubility

```bash
chemprop_train \
    --data_path rl/data/solubility/solubility.csv \
    --dataset_type regression \
    --target_column solubility \
    --features_path rl/data/solubility/solubility.npz \
    --no_features_scaling \
    --features_only \
    --num_folds 10 \
    --split_type scaffold_balanced \
    --metric mae \
    --extra_metrics r2 \
    --save_dir rl/models/solubility_mlp_rdkit_scaffold \
    --save_preds \
    --quiet
```

Results for _S. aureus_ (10-fold cross-validation, 8-core, 1-GPU machine):

| Model          | ROC-AUC         | PRC-AUC         | Time     |
|----------------|-----------------|-----------------|----------|
| Chemprop       | 0.831 +/- 0.019 | 0.477 +/- 0.033 | 76m, 42s |
| Chemprop-RDKit | 0.848 +/- 0.024 | 0.536 +/- 0.056 | 57m, 37s |
| MLP-RDKit      | 0.836 +/- 0.027 | 0.527 +/- 0.058 | 44m, 13s |

Results for solubility (10-fold cross-validation, 8-core, 1-GPU machine):

| Model          | MAE             | R^2             | Time     |
|----------------|-----------------|-----------------|----------|
| Chemprop       | 0.804 +/- 0.044 | 0.765 +/- 0.037 | 39m, 23s |
| Chemprop-RDKit | 0.760 +/- 0.038 | 0.788 +/- 0.034 | 39, 05s  |
| MLP-RDKit      | 0.800 +/- 0.041 | 0.768 +/- 0.032 | 31m, 32s |

## Chemical space

Analyze REAL versus WuXi from final SyntheMol generations.

## Multiparameter

Final (dynamic weights) versus the fixed weights below.

RL Chemprop-RDKit

```bash
for S_AUREUS_WEIGHT in 0.00 0.86 0.90 0.92 0.94 0.96 0.98 1.00
do
SOLUBILITY_WEIGHT="0$(echo "1.0 - S_AUREUS_WEIGHT" | bc)"
synthemol \
    --score_model_paths rl/models/s_aureus_chemprop_rdkit rl/models/solubility_chemprop_rdkit \
    --score_types chemprop chemprop \
    --score_fingerprint_types rdkit rdkit \
    --score_names 'S. aureus' 'Solubility' \
    --base_score_weights ${S_AUREUS_WEIGHT} ${SOLUBILITY_WEIGHT} \
    --chemical_spaces real wuxi \
    --building_blocks_paths rl/data/real/building_blocks.csv rl/data/wuxi/building_blocks.csv \
    --building_blocks_score_columns s_aureus_activity solubility \
    --reaction_to_building_blocks_paths rl/data/real/reaction_to_building_blocks.pkl rl/data/wuxi/reaction_to_building_blocks.pkl \
    --save_dir rl/generations/rl_chemprop_rdkit_s_aureus_${S_AUREUS_WEIGHT}_solubility_${SOLUBILITY_WEIGHT}_real_wuxi \
    --n_rollout 10000 \
    --search_type rl \
    --rl_model_type chemprop \
    --rl_model_fingerprint_type rdkit \
    --rl_model_paths rl/models/s_aureus_chemprop_rdkit/fold_0/model_0/model.pt rl/models/solubility_chemprop_rdkit/fold_0/model_0/model.pt \
    --rl_prediction_type regression \
    --use_gpu \
    --num_workers 8 \
    --wandb_project_name synthemol_rl \
    --wandb_run_name rl_chemprop_rdkit_s_aureus_${S_AUREUS_WEIGHT}_solubility_${SOLUBILITY_WEIGHT}_real_wuxi \
    --wandb_log
done
```

RL MLP-RDKit

```bash
for S_AUREUS_WEIGHT in 0.00 0.86 0.90 0.92 0.94 0.96 0.98 1.00
do
SOLUBILITY_WEIGHT="0$(echo "1.0 - S_AUREUS_WEIGHT" | bc)"
synthemol \
    --score_model_paths rl/models/s_aureus_chemprop_rdkit rl/models/solubility_chemprop_rdkit \
    --score_types chemprop chemprop \
    --score_fingerprint_types rdkit rdkit \
    --score_names 'S. aureus' 'Solubility' \
    --base_score_weights ${S_AUREUS_WEIGHT} ${SOLUBILITY_WEIGHT} \
    --chemical_spaces real wuxi \
    --building_blocks_paths rl/data/real/building_blocks.csv rl/data/wuxi/building_blocks.csv \
    --building_blocks_score_columns s_aureus_activity solubility \
    --reaction_to_building_blocks_paths rl/data/real/reaction_to_building_blocks.pkl rl/data/wuxi/reaction_to_building_blocks.pkl \
    --save_dir rl/generations/rl_mlp_rdkit_s_aureus_${S_AUREUS_WEIGHT}_solubility_${SOLUBILITY_WEIGHT}_real_wuxi \
    --n_rollout 10000 \
    --search_type rl \
    --rl_model_type mlp \
    --rl_model_fingerprint_type rdkit \
    --rl_model_paths rl/models/s_aureus_mlp_rdkit/fold_0/model_0/model.pt rl/models/solubility_mlp_rdkit/fold_0/model_0/model.pt \
    --rl_prediction_type regression \
    --wandb_project_name synthemol_rl \
    --wandb_run_name rl_mlp_rdkit_s_aureus_${S_AUREUS_WEIGHT}_solubility_${SOLUBILITY_WEIGHT}_real_wuxi \
    --wandb_log
done
```

## Dynamic temperature

Final (target similarity of 0.6) versus target similarities of 0.4, 0.5, 0.7, 0.8.

RL Chemprop-RDKit

```bash
for SIMILARITY_TARGET in 0.4 0.5 0.7 0.8
do
synthemol \
    --score_model_paths rl/models/s_aureus_chemprop_rdkit rl/models/solubility_chemprop_rdkit \
    --score_types chemprop chemprop \
    --score_fingerprint_types rdkit rdkit \
    --score_names 'S. aureus' 'Solubility' \
    --success_thresholds '>=0.5' '>=-4' \
    --chemical_spaces real wuxi \
    --building_blocks_paths rl/data/real/building_blocks.csv rl/data/wuxi/building_blocks.csv \
    --building_blocks_score_columns s_aureus_activity solubility \
    --reaction_to_building_blocks_paths rl/data/real/reaction_to_building_blocks.pkl rl/data/wuxi/reaction_to_building_blocks.pkl \
    --save_dir rl/generations/rl_chemprop_rdkit_s_aureus_solubility_dynamic_weights_real_wuxi_similarity_target_${SIMILARITY_TARGET} \
    --n_rollout 10000 \
    --search_type rl \
    --rl_model_type chemprop \
    --rl_model_fingerprint_type rdkit \
    --rl_model_paths rl/models/s_aureus_chemprop_rdkit/fold_0/model_0/model.pt rl/models/solubility_chemprop_rdkit/fold_0/model_0/model.pt \
    --rl_prediction_type regression \
    --rl_temperature_similarity_target ${SIMILARITY_TARGET} \
    --use_gpu \
    --num_workers 8 \
    --wandb_project_name synthemol_rl \
    --wandb_run_name rl_chemprop_rdkit_s_aureus_solubility_dynamic_weights_real_wuxi_similarity_target_${SIMILARITY_TARGET} \
    --wandb_log
done
```

RL MLP-RDKit

```bash
for SIMILARITY_TARGET in 0.4 0.5 0.7 0.8
do
synthemol \
    --score_model_paths rl/models/s_aureus_chemprop_rdkit rl/models/solubility_chemprop_rdkit \
    --score_types chemprop chemprop \
    --score_fingerprint_types rdkit rdkit \
    --score_names 'S. aureus' 'Solubility' \
    --success_thresholds '>=0.5' '>=-4' \
    --chemical_spaces real wuxi \
    --building_blocks_paths rl/data/real/building_blocks.csv rl/data/wuxi/building_blocks.csv \
    --building_blocks_score_columns s_aureus_activity solubility \
    --reaction_to_building_blocks_paths rl/data/real/reaction_to_building_blocks.pkl rl/data/wuxi/reaction_to_building_blocks.pkl \
    --save_dir rl/generations/rl_mlp_rdkit_s_aureus_solubility_dynamic_weights_real_wuxi_similarity_target_${SIMILARITY_TARGET} \
    --n_rollout 10000 \
    --search_type rl \
    --rl_model_type mlp \
    --rl_model_fingerprint_type rdkit \
    --rl_model_paths rl/models/s_aureus_mlp_rdkit/fold_0/model_0/model.pt rl/models/solubility_mlp_rdkit/fold_0/model_0/model.pt \
    --rl_prediction_type regression \
    --rl_temperature_similarity_target ${SIMILARITY_TARGET} \
    --wandb_project_name synthemol_rl \
    --wandb_run_name rl_mlp_rdkit_s_aureus_solubility_dynamic_weights_real_wuxi_similarity_target_${SIMILARITY_TARGET} \
    --wandb_log
done
```

## Exploration parameters (RL temperature vs MCTS explore weight)

RL with fixed temperatures of 0.01, 0.05, 0.1, 0.5, 1.0 versus MCTS with fixed explore weights of 0.5, 1.0, 5.0, 10.0,
50.0 (note: 10.0 is covered by the final model).

RL Chemprop-RDKit

```bash
for TEMPERATURE in 0.01 0.05 0.1 0.5 1.0
do
synthemol \
    --score_model_paths rl/models/s_aureus_chemprop_rdkit rl/models/solubility_chemprop_rdkit \
    --score_types chemprop chemprop \
    --score_fingerprint_types rdkit rdkit \
    --score_names 'S. aureus' 'Solubility' \
    --success_thresholds '>=0.5' '>=-4' \
    --chemical_spaces real wuxi \
    --building_blocks_paths rl/data/real/building_blocks.csv rl/data/wuxi/building_blocks.csv \
    --building_blocks_score_columns s_aureus_activity solubility \
    --reaction_to_building_blocks_paths rl/data/real/reaction_to_building_blocks.pkl rl/data/wuxi/reaction_to_building_blocks.pkl \
    --save_dir rl/generations/rl_chemprop_rdkit_s_aureus_solubility_dynamic_weights_real_wuxi_temperature_${TEMPERATURE} \
    --n_rollout 10000 \
    --search_type rl \
    --rl_model_type chemprop \
    --rl_model_fingerprint_type rdkit \
    --rl_model_paths rl/models/s_aureus_chemprop_rdkit/fold_0/model_0/model.pt rl/models/solubility_chemprop_rdkit/fold_0/model_0/model.pt \
    --rl_prediction_type regression \
    --rl_temperature_similarity_target -1 \
    --rl_base_temperature ${TEMPERATURE} \
    --use_gpu \
    --num_workers 8 \
    --wandb_project_name synthemol_rl \
    --wandb_run_name rl_chemprop_rdkit_s_aureus_solubility_dynamic_weights_real_wuxi_temperature_${TEMPERATURE} \
    --wandb_log
done
```

RL MLP-RDKit

```bash
for TEMPERATURE in 0.01 0.05 0.1 0.5 1.0
do
synthemol \
    --score_model_paths rl/models/s_aureus_chemprop_rdkit rl/models/solubility_chemprop_rdkit \
    --score_types chemprop chemprop \
    --score_fingerprint_types rdkit rdkit \
    --score_names 'S. aureus' 'Solubility' \
    --success_thresholds '>=0.5' '>=-4' \
    --chemical_spaces real wuxi \
    --building_blocks_paths rl/data/real/building_blocks.csv rl/data/wuxi/building_blocks.csv \
    --building_blocks_score_columns s_aureus_activity solubility \
    --reaction_to_building_blocks_paths rl/data/real/reaction_to_building_blocks.pkl rl/data/wuxi/reaction_to_building_blocks.pkl \
    --save_dir rl/generations/rl_mlp_rdkit_s_aureus_solubility_dynamic_weights_real_wuxi_temperature_${TEMPERATURE} \
    --n_rollout 10000 \
    --search_type rl \
    --rl_model_type mlp \
    --rl_model_fingerprint_type rdkit \
    --rl_model_paths rl/models/s_aureus_mlp_rdkit/fold_0/model_0/model.pt rl/models/solubility_mlp_rdkit/fold_0/model_0/model.pt \
    --rl_prediction_type regression \
    --rl_temperature_similarity_target -1 \
    --rl_base_temperature ${TEMPERATURE} \
    --wandb_project_name synthemol_rl \
    --wandb_run_name rl_mlp_rdkit_s_aureus_solubility_dynamic_weights_real_wuxi_temperature_${TEMPERATURE} \
    --wandb_log
done
```

MCTS

```bash
for EXPLORE_WEIGHT in 0.5 1.0 5.0 50.0
do
synthemol \
    --score_model_paths rl/models/s_aureus_chemprop_rdkit rl/models/solubility_chemprop_rdkit \
    --score_types chemprop chemprop \
    --score_fingerprint_types rdkit rdkit \
    --score_names 'S. aureus' 'Solubility' \
    --success_thresholds '>=0.5' '>=-4' \
    --chemical_spaces real wuxi \
    --building_blocks_paths rl/data/real/building_blocks.csv rl/data/wuxi/building_blocks.csv \
    --building_blocks_score_columns s_aureus_activity solubility \
    --reaction_to_building_blocks_paths rl/data/real/reaction_to_building_blocks.pkl rl/data/wuxi/reaction_to_building_blocks.pkl \
    --save_dir rl/generations/mcts_s_aureus_solubility_dynamic_weights_real_wuxi_explore_weight_${EXPLORE_WEIGHT} \
    --n_rollout 10000 \
    --search_type mcts \
    --explore_weight ${EXPLORE_WEIGHT} \
    --wandb_project_name synthemol_rl \
    --wandb_run_name mcts_s_aureus_solubility_dynamic_weights_real_wuxi_explore_weight_${EXPLORE_WEIGHT} \
    --wandb_log
done
```

## RL Training

RL Chemprop-RDKit (pretrained models, no RL training)

```bash
synthemol \
    --score_model_paths rl/models/s_aureus_chemprop_rdkit rl/models/solubility_chemprop_rdkit \
    --score_types chemprop chemprop \
    --score_fingerprint_types rdkit rdkit \
    --score_names 'S. aureus' 'Solubility' \
    --success_thresholds '>=0.5' '>=-4' \
    --chemical_spaces real wuxi \
    --building_blocks_paths rl/data/real/building_blocks.csv rl/data/wuxi/building_blocks.csv \
    --building_blocks_score_columns s_aureus_activity solubility \
    --reaction_to_building_blocks_paths rl/data/real/reaction_to_building_blocks.pkl rl/data/wuxi/reaction_to_building_blocks.pkl \
    --save_dir rl/generations/rl_chemprop_rdkit_s_aureus_solubility_dynamic_weights_real_wuxi_no_rl_training \
    --n_rollout 10000 \
    --search_type rl \
    --rl_model_type chemprop \
    --rl_model_fingerprint_type rdkit \
    --rl_model_paths rl/models/s_aureus_chemprop_rdkit/fold_0/model_0/model.pt rl/models/solubility_chemprop_rdkit/fold_0/model_0/model.pt \
    --rl_prediction_type regression \
    --rl_train_epochs 0 \
    --use_gpu \
    --num_workers 8 \
    --wandb_project_name synthemol_rl \
    --wandb_run_name rl_chemprop_rdkit_s_aureus_solubility_dynamic_weights_real_wuxi_no_rl_training \
    --wandb_log
```

RL Chemprop-RDKit (from scratch models, with RL training)

```bash
synthemol \
    --score_model_paths rl/models/s_aureus_chemprop_rdkit rl/models/solubility_chemprop_rdkit \
    --score_types chemprop chemprop \
    --score_fingerprint_types rdkit rdkit \
    --score_names 'S. aureus' 'Solubility' \
    --success_thresholds '>=0.5' '>=-4' \
    --chemical_spaces real wuxi \
    --building_blocks_paths rl/data/real/building_blocks.csv rl/data/wuxi/building_blocks.csv \
    --building_blocks_score_columns s_aureus_activity solubility \
    --reaction_to_building_blocks_paths rl/data/real/reaction_to_building_blocks.pkl rl/data/wuxi/reaction_to_building_blocks.pkl \
    --save_dir rl/generations/rl_chemprop_rdkit_s_aureus_solubility_dynamic_weights_real_wuxi_no_pretraining \
    --n_rollout 10000 \
    --search_type rl \
    --rl_model_type chemprop \
    --rl_model_fingerprint_type rdkit \
    --rl_prediction_type regression \
    --use_gpu \
    --num_workers 8 \
    --wandb_project_name synthemol_rl \
    --wandb_run_name rl_chemprop_rdkit_s_aureus_solubility_dynamic_weights_real_wuxi_no_pretraining \
    --wandb_log
```

RL MLP-RDKit (pretrained models, no RL training)

```bash
synthemol \
    --score_model_paths rl/models/s_aureus_chemprop_rdkit rl/models/solubility_chemprop_rdkit \
    --score_types chemprop chemprop \
    --score_fingerprint_types rdkit rdkit \
    --score_names 'S. aureus' 'Solubility' \
    --success_thresholds '>=0.5' '>=-4' \
    --chemical_spaces real wuxi \
    --building_blocks_paths rl/data/real/building_blocks.csv rl/data/wuxi/building_blocks.csv \
    --building_blocks_score_columns s_aureus_activity solubility \
    --reaction_to_building_blocks_paths rl/data/real/reaction_to_building_blocks.pkl rl/data/wuxi/reaction_to_building_blocks.pkl \
    --save_dir rl/generations/rl_mlp_rdkit_s_aureus_solubility_dynamic_weights_real_wuxi_no_rl_training \
    --n_rollout 10000 \
    --search_type rl \
    --rl_model_type mlp \
    --rl_model_fingerprint_type rdkit \
    --rl_model_paths rl/models/s_aureus_mlp_rdkit/fold_0/model_0/model.pt rl/models/solubility_mlp_rdkit/fold_0/model_0/model.pt \
    --rl_prediction_type regression \
    --rl_train_epochs 0 \
    --wandb_project_name synthemol_rl \
    --wandb_run_name rl_mlp_rdkit_s_aureus_solubility_dynamic_weights_real_wuxi_no_rl_training \
    --wandb_log
```

RL MLP-RDKit (from scratch models, with RL training)

```bash
synthemol \
    --score_model_paths rl/models/s_aureus_chemprop_rdkit rl/models/solubility_chemprop_rdkit \
    --score_types chemprop chemprop \
    --score_fingerprint_types rdkit rdkit \
    --score_names 'S. aureus' 'Solubility' \
    --success_thresholds '>=0.5' '>=-4' \
    --chemical_spaces real wuxi \
    --building_blocks_paths rl/data/real/building_blocks.csv rl/data/wuxi/building_blocks.csv \
    --building_blocks_score_columns s_aureus_activity solubility \
    --reaction_to_building_blocks_paths rl/data/real/reaction_to_building_blocks.pkl rl/data/wuxi/reaction_to_building_blocks.pkl \
    --save_dir rl/generations/rl_mlp_rdkit_s_aureus_solubility_dynamic_weights_real_wuxi_no_pretraining \
    --n_rollout 10000 \
    --search_type rl \
    --rl_model_type mlp \
    --rl_model_fingerprint_type rdkit \
    --rl_prediction_type regression \
    --wandb_project_name synthemol_rl \
    --wandb_run_name rl_mlp_rdkit_s_aureus_solubility_dynamic_weights_real_wuxi_no_pretraining \
    --wandb_log
```


## Analyze ablation experiments

Analyze each ablation experiment by determining the number of compounds that pass all the filtering steps.

```bash
for DIR in rl/generations/*/
do
chemfunc nearest_neighbor \
    --data_path ${DIR}/molecules.csv \
    --reference_data_path rl/data/s_aureus/s_aureus_hits.csv \
    --reference_name train_hits \
    --metric tversky

chemfunc nearest_neighbor \
    --data_path ${DIR}/molecules.csv \
    --reference_data_path rl/data/chembl/chembl.csv \
    --reference_name chembl \
    --metric tversky

python scripts/data/select_molecules.py \
    --data_path ${DIR}/molecules.csv \
    --save_molecules_path ${DIR}/hits.csv \
    --save_analysis_path ${DIR}/analysis.csv \
    --score_columns "S. aureus" "Solubility" \
    --score_comparators ">=0.5" ">=-4" \
    --novelty_thresholds 0.6 0.6 \
    --similarity_threshold 0.6 \
    --sort_column "S. aureus" \
    --descending
done
```
