# Generating Novel Antibiotics with SyntheMol-RL

Instructions for generating antibiotic candidates for _Acinetobacter baumannii_ using SyntheMol-RL from the
paper [TODO](TODO).

This includes instructions for processing antibiotics data, training antibacterial activity prediction models,
generating molecules with SyntheMol, and selecting candidates. Assumes relevant data has already been downloaded (
see [docs/README.md](README.md)).

* [Data](#data)
    + [Process S. aureus training data](#process-s-aureus-training-data)
    + [Process ChEMBL antibacterials](#process-chembl-antibacterials)
    + [Solubility data](#solubility-data)
* [Build bioactivity prediction models](#build-bioactivity-prediction-models)
    + [Compute RDKit features](#compute-rdkit-features)
    + [Train models](#train-models)
    + [Add labels to test preds](#add-labels-to-test-preds)
    + [Compute model scores for building blocks](#compute-model-scores-for-building-blocks)
* [Generate molecules with SyntheMol-RL](#generate-molecules-with-synthemol-rl)
* [Screen molecules with Chemprop-RDKit](#screen-molecules-with-chemprop-rdkit)
* [Select generated molecules](#select-generated-molecules)
* [Select screened molecules](#select-screened-molecules)
* [Visualize selected molecules](#visualize-selected-molecules)

## Data

### Process S. aureus training data

The training data consists of molecules tested for inhibitory activity against _Staphylococcus aureus_. The following
command processes the data to compute binary activity labels based on the inhibition values based on the mean (of two
replicates) normalized 16-hour optical density (OD) values.

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

### Process ChEMBL antibacterials

Download lists of known antibiotic-related compounds from ChEMBL using the following search terms. For each, click the
CSV download button, unzip the downloaded file, and rename the CSV file appropriately.

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

### Solubility data

Aqueous solubility data was obtained from [ADMET-AI](TODO), which preprocessed data from
the [Therapeutics Data Commons](https://tdcommons.ai/). This dataset contains 9,982 molecules with aqueous solubility
measurements in units of log mol/L. The data is saved to `data/solubility/solubility.csv`.

## Build bioactivity prediction models

Here, we build three binary classification bioactivity prediction models to predict antibiotic activity against _S.
aureus_ and to predict aqueous solubility. The three models are:

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

Solubility data

```bash
chemfunc save_fingerprints \
    --data_path rl/data/solubility/solubility.csv \
    --fingerprint_type rdkit \
    --save_path rl/data/solubility/solubility.npz
```

Time: 2 minutes, 0 seconds with an 8-core machine.

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

Chemprop for _S. aureus_

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
    --split_type cv \
    --metric mae \
    --extra_metrics r2 \
    --save_dir rl/models/solubility_chemprop \
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
    --split_type cv \
    --metric prc-auc \
    --extra_metrics auc \
    --save_dir rl/models/s_aureus_chemprop_rdkit \
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
    --split_type cv \
    --metric mae \
    --extra_metrics r2 \
    --save_dir rl/models/solubility_chemprop_rdkit \
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
    --split_type cv \
    --metric prc-auc \
    --extra_metrics auc \
    --save_dir rl/models/s_aureus_mlp_rdkit \
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
    --split_type cv \
    --metric mae \
    --extra_metrics r2 \
    --save_dir rl/models/solubility_mlp_rdkit \
    --save_preds \
    --quiet
```

Results for _S. aureus_ (10-fold cross-validation, 8-core, 1-GPU machine):

| Model          | ROC-AUC         | PRC-AUC         | Time     |
|----------------|-----------------|-----------------|----------|
| Chemprop       | 0.862 +/- 0.013 | 0.527 +/- 0.034 | 68m, 34s |
| Chemprop-RDKit | 0.875 +/- 0.013 | 0.570 +/- 0.050 | 79m, 38s |
| MLP-RDKit      | 0.873 +/- 0.018 | 0.553 +/- 0.040 | 60m, 09s |

Results for solubility (10-fold cross-validation, 8-core, 1-GPU machine):

| Model          | MAE             | R^2             | Time     |
|----------------|-----------------|-----------------|----------|
| Chemprop       | 0.693 +/- 0.027 | 0.806 +/- 0.021 | 35m, 59s |
| Chemprop-RDKit | 0.656 +/- 0.023 | 0.822 +/- 0.021 | 36m, 59s |
| MLP-RDKit      | 0.688 +/- 0.020 | 0.817 +/- 0.013 | 24m, 22s |

The Chemprop-RDKit models are the best for both _S. aureus_ and solubility. Therefore, those models are used to evaluate
building blocks and molecules.

### Add labels to test preds

When saving predictions, Chemprop does not save the labels. This command adds the labels to the predictions files (and
also corrects the formatting of SMILES strings).

```bash
for MODEL in chemprop chemprop_rdkit mlp_rdkit
do
python scripts/data/add_labels_to_test_preds.py \
    --true_path rl/data/s_aureus/s_aureus.csv \
    --preds_path rl/models/s_aureus_${MODEL} \
    --value_column s_aureus_activity

python scripts/data/add_labels_to_test_preds.py \
    --true_path rl/data/solubility/solubility.csv \
    --preds_path rl/models/solubility_${MODEL} \
    --value_column solubility
done
```

### Compute model scores for building blocks

After training, use the Chemprop-RDKit models to pre-compute scores of building blocks.

```bash
for CHEMICAL_SPACE in real wuxi
do
for MODEL in s_aureus solubility
do
chemprop_predict \
    --test_path rl/data/${CHEMICAL_SPACE}/building_blocks.csv \
    --smiles_column smiles \
    --preds_path rl/data/${CHEMICAL_SPACE}/building_blocks.csv \
    --checkpoint_dir rl/models/${MODEL}_chemprop_rdkit \
    --features_path rl/data/${CHEMICAL_SPACE}/building_blocks.npz \
    --no_features_scaling
done
done
```

Time with an 8-core, 1-GPU machine:

| Model      | Chemical Space | Time     |
|------------|----------------|----------|
| S. aureus  | REAL           | 20m, 11s |
| S. aureus  | WuXi           | 2m, 28s  |
| Solubility | REAL           | 20m, 06s |
| Solubility | WuXi           | 2m, 19s  |

## Generate molecules with SyntheMol-RL

Generate molecules with SyntheMol-RL using RL models for _S. aureus_ and solubility with dynamic multiparameter and
exploring REAL & WuXi chemical spaces.

RL Chemprop-RDKit

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
    --save_dir rl/generations/rl_chemprop_rdkit_s_aureus_solubility_dynamic_weights_real_wuxi \
    --n_rollout 10000 \
    --search_type rl \
    --rl_model_type chemprop_rdkit \
    --rl_model_paths rl/models/s_aureus_chemprop_rdkit/fold_0/model_0/model.pt rl/models/solubility_chemprop_rdkit/fold_0/model_0/model.pt \
    --rl_prediction_type regression \
    --use_gpu \
    --num_workers 8 \
    --replicate_rl \
    --wandb_project_name synthemol_rl \
    --wandb_run_name rl_chemprop_rdkit_s_aureus_solubility_dynamic_weights_real_wuxi \
    --wandb_log
```

RL MLP-RDKit

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
    --save_dir rl/generations/rl_mlp_rdkit_s_aureus_solubility_dynamic_weights_real_wuxi \
    --n_rollout 10000 \
    --search_type rl \
    --rl_model_type mlp_rdkit \
    --rl_model_paths rl/models/s_aureus_mlp_rdkit/fold_0/model_0/model.pt rl/models/solubility_mlp_rdkit/fold_0/model_0/model.pt \
    --rl_prediction_type regression \
    --replicate_rl \
    --wandb_project_name synthemol_rl \
    --wandb_run_name rl_mlp_rdkit_s_aureus_solubility_dynamic_weights_real_wuxi \
    --wandb_log
```

MCTS

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
    --save_dir rl/generations/mcts_s_aureus_solubility_dynamic_weights_real_wuxi \
    --n_rollout 10000 \
    --search_type mcts \
    --replicate_rl \
    --wandb_project_name synthemol_rl \
    --wandb_run_name mcts_s_aureus_solubility_dynamic_weights_real_wuxi \
    --wandb_log
```

## Screen molecules with Chemprop-RDKit

Chemprop-RDKit on random REAL 14 million and 100 and WuXi 7 million and 50.

```bash
for SPACE in real wuxi
do

if [ "${SPACE}" = "real" ]; then
  SIZES=("14m" "100")
else
  SIZES=("7m" "50")
fi

for SIZE in "${SIZES[@]}"
do
FILE_NAME="random_${SPACE}_${SIZE}"

chemfunc save_fingerprints \
    --data_path rl/data/${SPACE}/${FILE_NAME}.csv \
    --fingerprint_type rdkit \
    --save_path rl/data/${SPACE}/${FILE_NAME}.npz

for PROPERTY in s_aureus solubility
do
chemprop_predict \
    --test_path rl/data/${SPACE}/${FILE_NAME}.csv \
    --smiles_column smiles \
    --preds_path rl/screened/chemprop_rdkit_${FILE_NAME}_${PROPERTY}.csv \
    --checkpoint_dir rl/models/${PROPERTY}_chemprop_rdkit \
    --features_path rl/data/${SPACE}/${FILE_NAME}.npz \
    --no_features_scaling \
    --no_cache_mol
done

python -c "import pandas as pd
data = pd.read_csv('rl/screened/chemprop_rdkit_${FILE_NAME}_s_aureus.csv')
sol = pd.read_csv('rl/screened/chemprop_rdkit_${FILE_NAME}_solubility.csv')
data['solubility'] = sol['solubility']
data.to_csv('rl/screened/chemprop_rdkit_${FILE_NAME}.csv', index=False)"

done
done
```

## Select generated molecules

Compute similarity to training hits and ChEMBL antibiotics.

```bash
for MODEL in rl_chemprop_rdkit rl_mlp_rdkit mcts
do
chemfunc nearest_neighbor \
    --data_path rl/generations/${MODEL}_s_aureus_solubility_dynamic_weights_real_wuxi/molecules.csv \
    --reference_data_path rl/data/s_aureus/s_aureus_hits.csv \
    --reference_name train_hits \
    --metric tversky

chemfunc nearest_neighbor \
    --data_path rl/generations/${MODEL}_s_aureus_solubility_dynamic_weights_real_wuxi/molecules.csv \
    --reference_data_path rl/data/chembl/chembl.csv \
    --reference_name chembl \
    --metric tversky
done
```

Select the top 150 diverse, novel hit molecules. Hits are defined as _S. aureus_ >= 0.5 and solubility >= -4. Novelty is
defined as maximum 0.6 Tversky similarity to training hits and ChEMBL antibiotics. Diversity is defined as maximum 0.6
Tanimoto similarity to other selected molecules (maximum independent set). Final selection is the top 150 diverse, novel
hits molecules sorted by _S. aureus_ score.

```bash
for MODEL in rl_chemprop_rdkit rl_mlp_rdkit mcts
do
python scripts/data/select_molecules.py \
    --data_path rl/generations/${MODEL}_s_aureus_solubility_dynamic_weights_real_wuxi/molecules.csv \
    --save_molecules_path rl/selections/${MODEL}/molecules.csv \
    --save_analysis_path rl/selections/${MODEL}/analysis.csv \
    --score_columns "S. aureus" "Solubility" \
    --score_comparators ">=0.5" ">=-4" \
    --novelty_threshold 0.6 \
    --similarity_threshold 0.6 \
    --select_num 150 \
    --sort_column "S. aureus" \
    --descending
done
```

## Select screened molecules

Filter by hits, where hits are defined as _S. aureus_ >= 0.5 and solubility >= -4.

```bash
for NAME in real_14m wuxi_7m
do
chemfunc filter_molecules \
    --data_path rl/screened/chemprop_rdkit_random_${NAME}.csv \
    --save_path rl/screened/chemprop_rdkit_random_${NAME}_hits.csv \
    --filter_column "s_aureus_activity" \
    --min_value 0.5

chemfunc filter_molecules \
    --data_path rl/screened/chemprop_rdkit_random_${NAME}_hits.csv \
    --save_path rl/screened/chemprop_rdkit_random_${NAME}_hits.csv \
    --filter_column "solubility" \
    --min_value -4
done
```


REAL Output:

```
Original data size = 14,000,000
Data size after filtering with min_value 0.5 = 3,702
Final data size = 3,702

Original data size = 3,702
Data size after filtering with min_value -4.0 = 853
Final data size = 853
```

Note: The REAL random 14M has 3,702 (0.026%) with S. aureus >= 0.5; 11,371,677 (81.23%) with solubility >= -4; and 853 (0.006%) with both.


WuXi Output:

```
Original data size = 7,000,000
Data size after filtering with min_value 0.5 = 105,055
Final data size = 105,055

Original data size = 105,055
Data size after filtering with min_value -4.0 = 437
Final data size = 437
```

Note: The WuXi random 7M has 105,055 (1.50%) with S. aureus >= 0.5; 1,868,378 (26.69%) with solubility >= -4; and 437 (0.006%) with both.


Combine REAL and WuXi hits.

```bash
python -c "import pandas as pd
pd.concat([
    pd.read_csv('rl/screened/chemprop_rdkit_random_real_14m_hits.csv').assign(chemical_space='real'),
    pd.read_csv('rl/screened/chemprop_rdkit_random_wuxi_7m_hits.csv').assign(chemical_space='wuxi')
]).to_csv('rl/screened/chemprop_rdkit_random_real_wuxi_21m_hits.csv', index=False)"
```

Compute similarity to training hits and ChEMBL antibiotics.

```bash
chemfunc nearest_neighbor \
    --data_path rl/screened/chemprop_rdkit_random_real_wuxi_21m_hits.csv \
    --reference_data_path rl/data/s_aureus/s_aureus_hits.csv \
    --reference_name train_hits \
    --metric tversky

chemfunc nearest_neighbor \
    --data_path rl/screened/chemprop_rdkit_random_real_wuxi_21m_hits.csv \
    --reference_data_path rl/data/chembl/chembl.csv \
    --reference_name chembl \
    --metric tversky
```

Select the top 150 diverse, novel hit molecules. Hits are defined as _S. aureus_ >= 0.5 and solubility >= -4. Novelty is
defined as maximum 0.6 Tversky similarity to training hits and ChEMBL antibiotics. Diversity is defined as maximum 0.6
Tanimoto similarity to other selected molecules (maximum independent set). Final selection is the top 150 diverse, novel
hits molecules sorted by _S. aureus_ score.

```bash
python scripts/data/select_molecules.py \
    --data_path rl/screened/chemprop_rdkit_random_real_wuxi_21m_hits.csv \
    --save_molecules_path rl/selections/chemprop_rdkit/molecules.csv \
    --save_analysis_path rl/selections/chemprop_rdkit/analysis.csv \
    --score_columns "s_aureus_activity" "solubility" \
    --score_comparators ">=0.5" ">=-4" \
    --novelty_threshold 0.6 \
    --similarity_threshold 0.6 \
    --select_num 150 \
    --sort_column "s_aureus_activity" \
    --descending
```

Merge random molecules (no filtering).

```bash
python -c "import pandas as pd
pd.concat([
    pd.read_csv('rl/screened/chemprop_rdkit_random_real_100.csv').assign(chemical_space='real'),
    pd.read_csv('rl/screened/chemprop_rdkit_random_wuxi_50.csv').assign(chemical_space='wuxi')
]).to_csv('rl/selections/random/molecules.csv', index=False)"
```

## Visualize selected molecules

Create images of the structures of the selected molecules.

```bash
for NAME in rl_chemprop_rdkit rl_mlp_rdkit mcts chemprop_rdkit random
do
chemfunc visualize_molecules \
    --data_path rl/selections/${NAME}/molecules.csv \
    --save_dir rl/selections/${NAME}
done
```

Plot t-SNE of the selected molecules vs training hits and ChEMBL antibiotics.

```bash
chemfunc plot_tsne \
    --data_paths rl/data/chembl/chembl.csv \
    rl/data/s_aureus/s_aureus_hits.csv \
    rl/selections/random/molecules.csv \
    rl/selections/chemprop_rdkit/molecules.csv \
    rl/selections/mcts/molecules.csv \
    rl/selections/rl_mlp_rdkit/molecules.csv \
    rl/selections/rl_chemprop_rdkit/molecules.csv \
    --data_names chembl train_hits random chemprop_rdkit mcts rl_mlp_rdkit rl_chemprop_rdkit \
    --save_path rl/plots/selections_tsne.pdf \
    --point_size 2000
```

## Compute ADMET properties of selected molecules

Compute ADMET properties of selected molecules using [ADMET-AI](https://github.com/swansonk14/admet_ai).

Note: This requires installing ADMET-AI via `pip install admet_ai`.

```bash
for NAME in rl_chemprop_rdkit rl_mlp_rdkit mcts chemprop_rdkit random
do
admet_predict \
    --data_path rl/selections/${NAME}/molecules.csv
done
```

## Select final candidates based on availability and ClinTox score

Select 50 molecules from each model from the available molecules ranked by ClinTox score (lowest to highest). Note that this uses the `available.csv` file in the `rl/quotes` folder that is created based on the availability of the molecules in quotes from Enamine and WuXi

```bash
mkdir rl/purchase

python -c "import pandas as pd
available = pd.read_csv('rl/quotes/available.csv')
for name in ['rl_chemprop_rdkit', 'rl_mlp_rdkit', 'mcts', 'chemprop_rdkit', 'random']:
    data = pd.read_csv(f'rl/selections/{name}/molecules.csv')
    data = data[data['smiles'].isin(available['smiles'])]
    data = data.sort_values('ClinTox', ascending=True).head(50)
    data.to_csv(f'rl/purchase/{name}.csv', index=False)"
```
