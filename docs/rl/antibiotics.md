# Generating Novel Antibiotics with SyntheMol-RL

Instructions for generating antibiotic candidates for _Acinetobacter baumannii_ using SyntheMol-RL from the paper [TODO](TODO).

This includes instructions for processing antibiotics data, training antibacterial activity prediction models, generating molecules with SyntheMol, and selecting candidates. Assumes relevant data has already been downloaded (see [docs/README.md](README.md)).

- [Process antibiotics training data](#process-antibiotics-training-data)
- [Process ChEMBL antibacterials](#process-chembl-antibacterials)
- [Build bioactivity prediction models](#build-bioactivity-prediction-models)
  * [Train models](#train-models)
  * [Compute model scores for building blocks](#compute-model-scores-for-building-blocks)
- [Generate molecules with SyntheMol-RL](#generate-molecules-with-synthemol-rl)
- [Filter generated molecules](#filter-generated-molecules)
  * [Novelty](#novelty)
  * [Bioactivity](#bioactivity)
  * [Diversity](#diversity)
- [Map molecules to REAL IDs](#map-molecules-to-real-ids)
- [Predict toxicity](#predict-toxicity)


## Process S. aureus training data

The training data consists of molecules tested for inhibitory activity against _Staphylococcus aureus_. The following command processes the data to compute binary activity labels based on the inhibition values based on the mean (of two replicates) normalized 16-hour optical density (OD) values.

```bash
python scripts/data/process_data.py \
    --data_paths rl/data/s_aureus/s_aureus_raw.csv \
    --value_column "16h Normalized Mean" \
    --save_path rl/data/s_aureus/s_aureus.csv \
    --save_hits_path rl/data/s_aureus/s_aureus_hits.csv \
    --activity_column s_aureus_activity
```

Output:
```
Data size = 10,716
Mean activity = 0.8938596024636059
Std activity = 0.3028107500370331
Activity threshold of mean - 2 std = 0.2882381023895396
Number of hits = 1,149
Number of non-hits = 9,567

Full data size = 10,716
Data size after dropping non-conflicting duplicates = 10,660
Data size after dropping conflicting duplicates = 10,658

Final data size = 10,658
Number of hits = 1,137
Number of non-hits = 9,521
```


## Process ChEMBL antibacterials

Download lists of known antibiotic-related compounds from ChEMBL using the following search terms. For each, click the CSV download button, unzip the downloaded file, and rename the CSV file appropriately.

- [antibacterial](https://www.ebi.ac.uk/chembl/g/#search_results/compounds/query=antibacterial)
- [antibiotic](https://www.ebi.ac.uk/chembl/g/#search_results/compounds/query=antibiotic)

The results of the queries on November 8, 2023, are available in the Google Drive folder in the `data/chembl` subfolder.

The files are:

- `chembl_antibacterial.csv` with 611 molecules (590 with SMILES)
- `chembl_antibiotic.csv` with 636 molecules (591 with SMILES)

Merge these two files to form a single collection of antibiotic-related compounds.
```bash
python scripts/data/merge_chembl_downloads.py \
    --data_paths rl/data/chembl/chembl_antibacterial.csv rl/data/chembl/chembl_antibiotic.csv \
    --labels antibacterial antibiotic \
    --save_path rl/data/chembl/chembl.csv
```

The file `chembl.csv` contains 1,007 molecules.


## Build bioactivity prediction models


Here, we build three binary classification bioactivity prediction models to predict antibiotic activity against _S. aureus_. The three models are:

1. Chemprop: a graph neural network model
2. Chemprop-RDKit: a graph neural network model augmented with 200 RDKit features
3. MLP-RDKit: a multilayer perceptron using 200 RDKit features


### Compute RDKit features

Pre-compute the 200 RDKit features for the training data and the building blocks.

Training data
```bash
chemfunc save_fingerprints \
    --data_path rl/data/s_aureus/s_aureus.csv \
    --fingerprint_type rdkit \
    --save_path rl/data/s_aureus/s_aureus.npz
```

Time: 1 minute, 28 seconds with an 8-core machine.

Building blocks
```bash
for CHEMICAL_SPACE in real wuxi
do
chemfunc save_fingerprints \
    --data_path rl/data/${CHEMICAL_SPACE}/building_blocks.csv \
    --fingerprint_type rdkit \
    --save_path rl/data/${CHEMICAL_SPACE}/building_blocks.npz
done
```

Time: REAL = 10 minutes, 7 seconds; WuXi = 2 minutes, 3 seconds with an 8-core machine.


### Train models

For each model type, train 10 models using 10-fold cross-validation.

Chemprop
```bash
chemprop_train \
    --data_path rl/data/s_aureus/s_aureus.csv \
    --dataset_type classification \
    --target_column s_aureus_activity \
    --num_folds 10 \
    --split_type cv \
    --metric prc-auc \
    --extra_metrics auc \
    --save_dir rl/models/s_aureus_chemprop \
    --quiet
```

Time: 52 minutes, 19 seconds with an 8-core, 1-GPU machine.

Chemprop-RDKit
```bash
chemprop_train \
    --data_path rl/data/s_aureus/s_aureus.csv \
    --dataset_type classification \
    --target_column s_aureus_activity \
    --features_path rl/data/s_aureus/s_aureus.npz \
    --no_features_scaling \
    --num_folds 10 \
    --split_type cv \
    --metric prc-auc \
    --extra_metrics auc \
    --save_dir rl/models/s_aureus_chemprop_rdkit \
    --quiet
```

Time: 51 minutes, 46 seconds with an 8-core, 1-GPU machine.

MLP-RDKit
```bash
chemprop_train \
    --data_path rl/data/s_aureus/s_aureus.csv \
    --dataset_type classification \
    --target_column s_aureus_activity \
    --features_path rl/data/s_aureus/s_aureus.npz \
    --no_features_scaling \
    --features_only \
    --num_folds 10 \
    --split_type cv \
    --metric prc-auc \
    --extra_metrics auc \
    --save_dir rl/models/s_aureus_mlp_rdkit \
    --quiet
```

Time: 42 minutes, 34 seconds with an 8-core, 1-GPU machine.

| Model          | ROC-AUC         | PRC-AUC         | Time     |
|----------------|-----------------|-----------------|----------|
| Chemprop       | 0.861 +/- 0.013 | 0.531 +/- 0.042 | 52m, 19s |
| Chemprop-RDKit | 0.874 +/- 0.017 | 0.575 +/- 0.046 | 51m, 46s |
| MLP-RDKit      | 0.873 +/- 0.019 | 0.554 +/- 0.043 | 42m, 34s |


### Compute model scores for building blocks

After training, use the model to pre-compute scores of building blocks.

Chemprop
```bash
for CHEMICAL_SPACE in real wuxi
do
chemprop_predict \
    --test_path rl/data/${CHEMICAL_SPACE}/building_blocks.csv \
    --checkpoint_dir rl/models/s_aureus_chemprop \
    --preds_path rl/models/s_aureus_chemprop/${CHEMICAL_SPACE}_building_blocks.csv
done
```

Time: REAL = 12 minutes, 21 seconds; WuXi = 6 minutes, 29 seconds with an 8-core, 1-GPU machine.

Chemprop-RDKit
```bash
for CHEMICAL_SPACE in real wuxi
do
chemprop_predict \
    --test_path rl/data/${CHEMICAL_SPACE}/building_blocks.csv \
    --checkpoint_dir rl/models/s_aureus_chemprop_rdkit \
    --preds_path rl/models/s_aureus_chemprop_rdkit/${CHEMICAL_SPACE}_building_blocks.csv \
    --features_path rl/data/${CHEMICAL_SPACE}/building_blocks.npz \
    --no_features_scaling
done
```

Time: REAL = 12 minutes, 28 seconds; WuXi = 1 minute, 58 seconds with an 8-core, 1-GPU machine.

MLP-RDKit
```bash
for CHEMICAL_SPACE in real wuxi
do
chemprop_predict \
    --test_path rl/data/${CHEMICAL_SPACE}/building_blocks.csv \
    --checkpoint_dir rl/models/s_aureus_mlp_rdkit \
    --preds_path rl/models/s_aureus_mlp_rdkit/${CHEMICAL_SPACE}_building_blocks.csv \
    --features_path rl/data/${CHEMICAL_SPACE}/building_blocks.npz \
    --no_features_scaling
done
```

Time: REAL = 12 minutes, 27 seconds; WuXi = 2 minutes, 29 seconds with an 8-core, 1-GPU machine.


TODO: solubility predictions (or just all ADMET with ADMET-AI)


## Generate molecules with SyntheMol-RL

Generate molecules with SyntheMol-RL.


### Final generations

RL models for _S. aureus_ and solubility dynamic multiparameter REAL & WuXi

TODO: Use ADMET-AI solubility model?
TODO: implement weight loading for MLP

```bash
for RL_MODEL_TYPE in mlp_rdkit chemprop_rdkit
do
synthemol \
    --model_paths rl/models/s_aureus_chemprop_rdkit TODO:solubility \
    --model_types chemprop TODO:solublity_model \
    --fingerprint_types rdkit rdkit \
    --model_names 'S. aureus' 'Solubility' \
    --success_thresholds '>=0.5' '>=-4' \
    --chemical_spaces real wuxi \
    --building_blocks_paths rl/models/s_aureus_chemprop_rdkit/real_building_blocks.csv rl/models/s_aureus_chemprop_rdkit/wuxi_building_blocks.csv \
    --building_blocks_score_columns activity TODO:solubility \
    --save_dir rl/generations/rl_${RL_MODEL_TYPE}_s_aureus_solubility_dynamic_weights_real_wuxi \
    --n_rollout TODO:rollouts \
    --search_type rl \
    --rl_model_type ${RL_MODEL_TYPE} \
    --rl_model_paths rl/models/s_aureus_${RL_MODEL_TYPE}/fold_0/model_0/model.pt TODO:solubility \
    --rl_prediction_type regression \
    --wandb_project_name synthemol_rl \
    --wandb_run_name rl_${RL_MODEL_TYPE}_s_aureus_solubility_dynamic_weights_real_wuxi \
    --wandb_log
done
```


## Filter generated molecules


### Novelty


### Bioactivity


### Diversity


## Map molecules to REAL IDs


## Predict toxicity

