# SyntheMol

SyntheMol is a generative AI method for designing easily synthesizable and structurally novel drug candidates. SyntheMol consists of a Monte Carlo tree search (MCTS) guided by a molecular property prediction model that searches the [Enamine REAL Space](https://enamine.net/compound-collections/real-compounds/real-space-navigator) for molecules with a desired property. REAL Space molecules are constructed by combining off-the-shelf molecular building blocks with known chemical reactions. SyntheMol efficiently searches the space of combinations of building blocks and reactions to find molecules with a desired property. Here, we apply SyntheMol to design novel antibiotics for the Gram-negative bacterium _Acinetobacter baumannii_.

TODO: Change all instances of fragments or reagents in the code to building blocks. Change rf to random_forest as model type.

- [Installation](#installation)
- [Data](#data)
- [Process REAL data](#process-real-data)
  * [Process building blocks](#process-building-blocks)
    + [Convert building blocks from SDF to SMILES](#convert-building-blocks-from-sdf-to-smiles)
    + [Remove salts from building blocks](#remove-salts-from-building-blocks)
  * [Process REAL Space molecules](#process-real-space-molecules)
    + [Download REAL Space](#download-real-space)
    + [Map REAL reactions to building blocks](#map-real-reactions-to-building-blocks)
    + [Count REAL reactions and building blocks](#count-real-reactions-and-building-blocks)
    + [Sample REAL molecules](#sample-real-molecules)
    + [Count REAL molecules limited to selected reactions and building blocks](#count-real-molecules-limited-to-selected-reactions-and-building-blocks)
- [Process antibiotics training data](#process-antibiotics-training-data)
- [Process ChEMBL antibacterials](#process-chembl-antibacterials)
- [Build property prediction models](#build-property-prediction-models)
  * [Train models](#train-models)
  * [Map building blocks to model scores](#map-building-blocks-to-model-scores)
- [Generate molecules with SyntheMol](#generate-molecules-with-synthemol)
- [Assess generated molecules](#assess-generated-molecules)
- [Select generated molecules](#select-generated-molecules)
  * [Filter by novelty](#filter-by-novelty)
  * [Filter by model score](#filter-by-model-score)
  * [Filter by diversity](#filter-by-diversity)
- [Assess selected molecules](#assess-selected-molecules)
- [Visualize molecules](#visualize-molecules)
- [Map molecules to REAL IDs](#map-molecules-to-real-ids)
- [Predict toxicity](#predict-toxicity)
- [Plots](#plots)

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

Install [chem_utils](https://github.com/swansonk14/chem_utils), which contains useful scripts and functions for working with small molecules.
```bash
cd ..
git clone git@github.com:swansonk14/chem_utils.git
cd chem_utils
pip install -e .
cd ../combinatorial_antibiotics
```

[//]: # (&#40;Note: Had some issues with pip so commented out the `install_requires` lines in `setup.py`.&#41;)

Install [chemprop](https://github.com/chemprop/chemprop), which contain a graph neural network for molecular property prediction.
```bash
cd ..
git clone git@github.com:chemprop/chemprop.git
cd chemprop
pip install -e .
cd ../combinatorial_antibiotics
```

[//]: # (&#40;Note: Had some issues with pip so commented out the `install_requires` and `extras_require` lines in `setup.py`.&#41;)


## Data

All raw and processed data is available in this Google Drive folder: https://drive.google.com/drive/folders/1VLPPUbY_FTKMjlXgRm09bPSSms206Dce?usp=share_link.


## Process REAL data

Process the Enamine REAL Space data consisting of building blocks, reactions, and molecules. The sections below detail the steps for processing the building blocks and molecules. A PDF with the 169 reactions is provided in the Google Drive folder, and SMARTS of the 13 reactions used here are provided in `real_reactions.py`.


### Process building blocks

Process the REAL building blocks.


#### Convert building blocks from SDF to SMILES

The `building_blocks.sdf` file (from the Google Drive folder) is the 2021 q3-4 version containing 138,085 molecules.

Convert the building blocks from SDF to SMILES using [chem_utils](https://github.com/swansonk14/chem_utils).
```bash
python sdf_to_smiles.py \
    --data_path ../../combinatorial_antibiotics/data/building_blocks.sdf \
    --save_path ../../combinatorial_antibiotics/data/building_blocks_raw.csv \
    --properties Reagent_ID Catalog_ID
```

All 138,085 molecules were successfully converted from SDF to SMILES, and among those 134,609 are unique molecules.

Note: The SMILES are likely not all unique because they do not include stereochemistry even though it is annotated with the `CFG` tag in the SDF file.


#### Remove salts from building blocks

Remove the salts from the building blocks using [chem_utils](https://github.com/swansonk14/chem_utils). This will also canonicalize the SMILES using RDKit's canonicalization method.
```bash
python canonicalize_smiles.py \
    --data_path ../../combinatorial_antibiotics/data/building_blocks_raw.csv \
    --save_path ../../combinatorial_antibiotics/data/building_blocks.csv \
    --remove_salts \
    --delete_disconnected_mols
```

This removes 25 molecules whose salts cannot be stripped, leaving 138,060 molecules (unique IDs, duplicate SMILES).

Note: This step is crucial to prevent errors in running reactions. Salts can cause reactions to create products that are the same as the reactants, leading to undesired infinite loops during molecule generation.

Deduplicate by SMILES using [chem_utils](https://github.com/swansonk14/chem_utils).
```bash
python deduplicate_smiles.py \
    --data_path ../../combinatorial_antibiotics/data/building_blocks.csv \
    --save_path ../../combinatorial_antibiotics/data/building_blocks_unique.csv
```

This leaves 132,479 unique molecules.


### Process REAL Space molecules

Process the REAL Space molecules.


#### Download REAL Space

Download the building block and reaction IDs used for the full REAL Space of 31 billion compounds (2022 q1-2 version downloaded on August 30, 2022).
```bash
lftp -c "open -p 21 -u username,password ftp-rdb-fr.chem-space.com; mirror -c --parallel=16 . data/real_space"
```

In the above command, replace `username` and `password` with the appropriate values provided by Enamine.


#### Map REAL reactions to building blocks

Determine which building blocks are valid in which REAL reactions.
```bash
python map_real_reactions_to_reagents.py \
    --data_dir data/real_space \
    --save_path data/reaction_to_building_blocks.json \
    --parallel
```

Total number of molecules = 31,507,987,117


#### Count REAL reactions and building blocks

Determine which reactions and building blocks are most common in REAL space.
```bash
python count_real_space.py \
    --data_dir data/real_space \
    --save_dir data/real_space_counts
```


#### Sample REAL molecules

Randomly sample 25,000 REAL Space molecules. Used for analysis of a representative sample of REAL Space molecules.
```bash
python sample_real_space.py \
    --data_dir data/real_space \
    --save_path data/random_real.csv \
    --num_molecules 25000 \
    --parallel
```


#### Count REAL molecules limited to selected reactions and building blocks

Count REAL space when limiting to the selected 13 reactions and the building blocks that we have after processing.
```bash
python count_real_space.py \
    --data_dir data/Enamine_REAL_space \
    --save_dir data/Enamine_REAL_space_counts_with_selected_fragments_and_reactions \
    --fragment_path data/2021q3-4_Enamine_REAL_reagents_SMILES_no_salts.csv \
    --only_selected_reactions
```

Number of files = 2,993
Number of building blocks = 138,060 (132,479 unique molecules)
Number of selected reactions = 13
Total number of molecules = 31,507,987,117
Total number of molecules with selected building blocks/reactions = 29,575,293,692


## Process antibiotics training data

The antibiotics training data consists of three libraries of molecules tested against the Gram-negative bacterium _Acinetobacter baumannii_. Here, we merge the three libraries into a single file and determine which molecules are hits (i.e., active against _A. baumannii_). The data referred to here is in the `Data/1. Training Data` subfolder from the Google Drive folder.

Merge the three libraries into a single file and determine which molecules are hits.
```bash
python process_data.py \
    --data_paths data/library_1.csv data/library_2.csv data/library_3.csv \
    --save_path data/antibiotics.csv \
    --save_hits_path data/antibiotics_hits.csv
```

Output:
```
library_1
Data size = 2,371
Mean activity = 0.9759337199493885
Std activity = 0.2673771821337539
Activity threshold of mean - 2 std = 0.4411793556818807
Number of hits = 130
Number of non-hits = 2,241

library_2
Data size = 6,680
Mean activity = 0.9701813713618264
Std activity = 0.17623388953330232
Activity threshold of mean - 2 std = 0.6177135922952217
Number of hits = 294
Number of non-hits = 6,386

library_3
Data size = 5,376
Mean activity = 0.9980530249533108
Std activity = 0.13303517604898704
Activity threshold of mean - 2 std = 0.7319826728553367
Number of hits = 112
Number of non-hits = 5,264

Full data size = 14,427
Data size after dropping non-conflicting duplicates = 13,594
Data size after dropping conflicting duplicates = 13,524

Final data size = 13,524
Number of hits = 470
Number of non-hits = 13,054
```


## Process ChEMBL antibacterials

Download lists of known antibiotic-related compounds from ChEMBL using the following search terms. For each, click the CSV download button, unzip the downloaded file, and rename the CSV file appropriately.

- [antibacterial](https://www.ebi.ac.uk/chembl/g/#search_results/compounds/query=antibacterial)
- [antibiotic](https://www.ebi.ac.uk/chembl/g/#search_results/compounds/query=antibiotic)

The results of the queries on August 8, 2022, are available in the Google Drive folder in the `Data/2. ChEMBL` subfolder.

The files are:

- `chembl_antibacterial.csv` with 604 molecules (589 with SMILES)
- `chembl_antibiotic.csv` with 636 molecules (587 with SMILES)

Merge these two files to form a single collection of antibiotic-related compounds.
```bash
python merge_chembl_downloads.py \
    --data_paths data/chembl/chembl_antibacterial.csv data/chembl/chembl_antibiotic.csv \
    --labels antibacterial antibiotic \
    --save_path data/chembl/chembl.csv
```

The file `chembl.csv` contains 1,005 molecules.


## Build property prediction models

Here, we build three binary classification molecular property prediction models to predict antibiotic activity against _A. baumannii_. The three models are:

1. Chemprop: a graph neural network model
2. Chemprop-RDKit: a graph neural network model augmented with 200 RDKit features
3. Random forest: a random forest model using 200 RDKit features


### Train models

For each model type, we trained 10 models using 10-fold cross-validation on the training data. Each ensemble of 10 models took less than 90 minutes to train on a 16-CPU machine. Trained models are available in the `Models` subfolder of the Google Drive folder.

Chemprop
```bash
python train_model.py \
    --data_path data/antibiotics.csv \
    --save_dir models/chemprop \
    --dataset_type classification \
    --model_type chemprop \
    --property_column activity \
    --num_models 10
```

Chemprop-RDKit
```bash
python train_model.py \
    --data_path data/antibiotics.csv \
    --save_dir models/chemprop_rdkit \
    --model_type chemprop \
    --dataset_type classification \
    --fingerprint_type rdkit \
    --property_column activity \
    --num_models 10
```

Random forest
```bash
python train_model.py \
    --data_path data/antibiotics.csv \
    --save_dir models/random_forest \
    --model_type rf \
    --dataset_type classification \
    --fingerprint_type rdkit \
    --property_column activity \
    --num_models 10
```

| Model          | ROC-AUC         | PRC-AUC         |
|----------------|-----------------|-----------------|
| Chemprop       | 0.803 +/- 0.036 | 0.354 +/- 0.086 |
| Chemprop-RDKit | 0.827 +/- 0.028 | 0.388 +/- 0.078 |
| Random forest  | 0.835 +/- 0.035 | 0.401 +/- 0.099 |


### Map building blocks to model scores

In order to speed up the generative model, we pre-compute the scores of all building blocks for each model.

Chemprop
```bash
python map_fragments_to_model_scores.py \
    --fragment_path data/building_blocks.csv \
    --model_path models/chemprop \
    --save_path models/chemprop/building_block_scores.json \
    --model_type chemprop
```

Chemprop-RDKit
```bash
python map_fragments_to_model_scores.py \
    --fragment_path data/building_blocks.csv \
    --model_path models/chemprop_rdkit \
    --save_path models/chemprop_rdkit/building_block_scores.json \
    --model_type chemprop \
    --fingerprint_type rdkit
```

Random forest with RDKit
```bash
python map_fragments_to_model_scores.py \
    --fragment_path data/building_blocks.csv \
    --model_path models/random_forest \
    --save_path models/random_forest/building_block_scores.json \
    --model_type rf \
    --fingerprint_type rdkit
```


## Generate molecules with SyntheMol

Here, we apply SyntheMol to generate molecules using a Monte Carlo tree search (MCTS). We limit the search space to molecules that can be formed using a single reaction (as opposed to multi-step reactions) to ensure easy, quick, and cheap chemical synthesis.

Prior to applying SyntheMol to generate antibiotics, we first run a simulation study with the property cLogP, which can be computed on the generated molecules, thereby providing an _in silico_ measure of the efficacy of SyntheMol. Instructions for running the simulation study and analyzing the results are in `simulation.md`.

Below, we run SyntheMol to generate antibiotics. We run SyntheMol three times, each time using a different property prediction model to score molecules and guide the search.

Chemprop
```bash
python tree_search.py \
    --model_path models/chemprop \
    --fragment_path data/building_blocks.csv \
    --reaction_to_reagents_path data/reaction_to_building_blocks.json \
    --fragment_to_model_score_path models/chemprop/building_block_scores.json \
    --save_dir generations/chemprop \
    --search_type mcts \
    --model_type chemprop \
    --n_rollout 20000 \
    --fragment_diversity \
    --max_reactions 1
```

Chemprop-RDKit
```bash
python tree_search.py \
    --model_path models/chemprop_rdkit \
    --fragment_path data/building_blocks.csv \
    --reaction_to_reagents_path data/reaction_to_building_blocks.json \
    --fragment_to_model_score_path models/chemprop_rdkit/building_block_scores.json \
    --save_dir generations/chemprop_rdkit \
    --search_type mcts \
    --model_type chemprop \
    --fingerprint_type rdkit \
    --n_rollout 20000 \
    --fragment_diversity \
    --max_reactions 1
```

Random forest
```bash
python tree_search.py \
    --model_path models/random_forest \
    --fragment_path data/building_blocks.csv \
    --reaction_to_reagents_path data/reaction_to_building_blocks.json \
    --fragment_to_model_score_path models/random_forest/building_block_scores.json \
    --save_dir generations/random_forest \
    --search_type mcts \
    --model_type rf \
    --fingerprint_type rdkit \
    --n_rollout 20000 \
    --fragment_diversity \
    --max_reactions 1
```


## Assess generated molecules

Optionally, assess various quantities about the generated molecules such as novelty, score, and diversity.
```bash
#!/bin/bash

for NAME in chemprop chemprop_rdkit random_forest
do
python assess_generated_molecules.py \
    --data_path generations/${NAME}/molecules.csv \
    --save_dir generations/${NAME} \
    --reference_paths data/antibiotics_hits.csv data/chembl.csv
done
```


## Select generated molecules

Run three filtering steps to select molecules for experimental testing based on novelty, score, and diversity.


### Filter by novelty

Filter for molecules that are structurally novel, i.e., dissimilar from the training hits and the ChEMBL antibacterials.

First, Compute the nearest neighbor Tversky similarity between the generated molecules and the training hits or ChEMBL antibacterials to determine the novelty of each generated molecule. Tversky similarity is used to determine what proportion of the functional groups of known antibacterials are contained within the generated molecules. These functions are from [chem_utils](https://github.com/swansonk14/chem_utils).

Compute Tversky similarity to train hits.
```bash
#!/bin/bash

for NAME in chemprop chemprop_rdkit random_forest
do
python nearest_neighbor.py \
    --data_path ../../combinatorial_antibiotics/generations/${NAME}/molecules.csv \
    --reference_data_path ../../combinatorial_antibiotics/data/antibiotics_hits.csv \
    --reference_name antibiotics_hits \
    --metrics tversky
done
```

Compute Tversky similarity to ChEMBL antibacterials.
```bash
#!/bin/bash

for NAME in chemprop chemprop_rdkit random_forest
do
python nearest_neighbor.py \
    --data_path ../../combinatorial_antibiotics/generations/${NAME}/molecules.csv \
    --reference_data_path ../../combinatorial_antibiotics/data/chembl.csv \
    --reference_name chembl \
    --metrics tversky
done
```

Now, filter to only keep molecules with nearest neighbor Tverksy similarity to the train hits or ChEMBL antibacterials <= 0.5 using [chem_utils](https://github.com/swansonk14/chem_utils).

Filter by Tversky similarity to train hits.
```bash
#!/bin/bash

for NAME in chemprop chemprop_rdkit random_forest
do
python filter_molecules.py \
    --data_path ../../combinatorial_antibiotics/generations/${NAME}/molecules.csv \
    --save_path ../../combinatorial_antibiotics/generations/${NAME}/molecules_train_sim_below_0.5.csv \
    --filter_column antibiotics_hits_tversky_nearest_neighbor_similarity \
    --max_value 0.5
done
```

Filter by Tversky similarity to ChEMBL antibacterials.
```bash
#!/bin/bash

for NAME in chemprop chemprop_rdkit random_forest
do
python filter_molecules.py \
    --data_path ../../combinatorial_antibiotics/generations/${NAME}/molecules_train_sim_below_0.5.csv \
    --save_path ../../combinatorial_antibiotics/generations/${NAME}/molecules_train_sim_below_0.5_chembl_sim_below_0.5.csv \
    --filter_column chembl_tversky_nearest_neighbor_similarity \
    --max_value 0.5
done
```


### Filter by model score

Filter for high-scoring molecules by only keeping molecules with the top 20% of model scores using [chem_utils](https://github.com/swansonk14/chem_utils).
```bash
#!/bin/bash

for NAME in chemprop chemprop_rdkit random_forest
do
python filter_molecules.py \
    --data_path ../../combinatorial_antibiotics/generations/${NAME}/molecules_train_sim_below_0.5_chembl_sim_below_0.5.csv \
    --save_path ../../combinatorial_antibiotics/generations/${NAME}/molecules_train_sim_below_0.5_chembl_sim_below_0.5_top_20_percent.csv \
    --filter_column score \
    --top_proportion 0.2
done
```


### Filter by diversity

Filter for diverse molecules by clustering molecules based on their Morgan fingerprint and only keeping the top scoring molecule from each cluster using [chem_utils](https://github.com/swansonk14/chem_utils).

Cluster molecules based on their Morgan fingerprint.
```bash
#!/bin/bash

for NAME in chemprop chemprop_rdkit random_forest
do
python cluster_molecules.py \
    --data_path ../../combinatorial_antibiotics/generations/${NAME}/molecules_train_sim_below_0.5_chembl_sim_below_0.5_top_20_percent.csv \
    --num_clusters 50
done
```

Select the top scoring molecule from each cluster.
```bash
#!/bin/bash

for NAME in chemprop chemprop_rdkit random_forest
do
python select_from_clusters.py \
    --data_path ../../combinatorial_antibiotics/generations/${NAME}/molecules_train_sim_below_0.5_chembl_sim_below_0.5_top_20_percent.csv \
    --save_path ../../combinatorial_antibiotics/generations/${NAME}/molecules_train_sim_below_0.5_chembl_sim_below_0.5_top_20_percent_selected_50.csv \
    --value_column score
done
```

We now have 50 molecules selected from each model's generations that meet our desired novelty, score, and diversity criteria.


## Assess selected molecules

Optionally, assess various quantities about the selected molecules.

```bash
#!/bin/bash

for NAME in chemprop chemprop_rdkit random_forest
do
python assess_generated_molecules.py \
    --data_path generations/${NAME}/molecules_train_sim_below_0.5_chembl_sim_below_0.5_top_20_percent_selected_50.csv \
    --save_dir generations/${NAME}/analysis_molecules_train_sim_below_0.5_chembl_sim_below_0.5_top_20_percent_selected_50 \
    --reference_paths data/antibiotics_hits.csv data/chembl.csv
done
```


## Visualize molecules

Optionally, visualize the selected molecules using [chem_utils](https://github.com/swansonk14/chem_utils).
```bash
#!/bin/bash

for NAME in chemprop chemprop_rdkit random_forest
do
python visualize_molecules.py \
    --data_path ../../combinatorial_antibiotics/generations/${NAME}/molecules_train_sim_below_0.5_chembl_sim_below_0.5_top_20_percent_selected_50.csv \
    --save_dir ../../combinatorial_antibiotics/generations/${NAME}
done
```


## Map molecules to REAL IDs

TODO: change to use all possible IDs for each SMILES

Map generated molecules to REAL IDs in the format expected by Enamine to enable a lookup in the Enamine REAL Space database.
```bash
#!/bin/bash

for NAME in chemprop chemprop_rdkit random_forest
do
python map_generated_molecules_to_real_ids.py \
    --data_path generations/${NAME}/molecules_train_sim_below_0.5_chembl_sim_below_0.5_top_20_percent_selected_50.csv \
    --save_dir generations/${NAME}/molecules_train_sim_below_0.5_chembl_sim_below_0.5_top_20_percent_selected_50_real_ids
done
```


## Predict toxicity

Predict the toxicity of the synthesized molecules by training a toxicity prediction model on the ClinTox dataset from [MoleculeNet](https://moleculenet.org/). ClinTox is a clinical trial toxicity binary classification dataset with 1,478 molecules with 112 (7.58%) toxic.

Train toxicity model using the `CT_TOX` property of the ClinTox dataset. (Note that the ClinTox properties `CT_TOX` and `FDA_APPROVED` are nearly always identical, so we only use one.)
```bash
python train_model.py \
    --data_path data/clintox.csv \
    --save_dir models/clintox \
    --model_type chemprop \
    --dataset_type classification \
    --fingerprint_type rdkit \
    --property_column CT_TOX \
    --num_models 10
```

The model has an ROC-AUC of 0.881 +/- 0.045 and a PRC-AUC of 0.514 +/- 0.141 across 10-fold cross-validation.

Make toxicity predictions on the synthesized molecules. The list of successfully synthesized generated molecules is in the `Data/8. Synthesized` subfolder of the Google Drive folder.
```bash
python predict_model.py \
    --data_path generations/synthesized.csv \
    --model_path models/clintox \
    --model_type chemprop \
    --fingerprint_type rdkit \
    --average_preds
```

## Plots

Instructions for generating plots analysing the data and results are in `plots/README.md`.
