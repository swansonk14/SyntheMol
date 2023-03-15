# SyntheMol: Generative AI for Drug Discovery

SyntheMol is a generative AI method for designing structurally novel and diverse drug candidates with predicted bioactivity that are easy to synthesize.

SyntheMol consists of a Monte Carlo tree search (MCTS) that explores a combinatorial chemical space consisting of molecular building blocks and chemical reactions. The MCTS is guided by a bioactivity prediction AI model, such as a graph neural network or a random forest. Currently, SyntheMol is designed to use 139,444 building blocks and 13 chemical reactions from the [Enamine REAL Space](https://enamine.net/compound-collections/real-compounds/real-space-navigator), which can produce over 30 billion molecules. However, SyntheMol can be easily adapted to use any set of building blocks and reactions.

SyntheMol is described in the paper [TODO](TODO), where we applied SyntheMol to design novel antibiotic candidates for the Gram-negative bacterium _Acinetobacter baumannii_. Full details for reproducing the results in the paper are provided in the [docs](docs) directory.

If you use SyntheMol in your research, please cite:

```
TODO
```

## Table of contents

* [Installation](#installation)
* [Combinatorial chemical space](#combinatorial-chemical-space)
* [Bioactivity prediction model](#bioactivity-prediction-model)
  + [Train model](#train-model)
  + [Pre-compute building block scores](#pre-compute-building-block-scores)
* [Generate molecules](#generate-molecules)
* [Filter generated molecules](#filter-generated-molecules)
  + [Novelty](#novelty)
  + [Bioactivity](#bioactivity)
  + [Diversity](#diversity)


## Installation

TODO: put SyntheMol (and chem_utils) on pip

TODO: remove environment.yml?

TODO: rdkit?

Create conda environment.
```bash
conda create --name SyntheMol python=3.10
conda activate SyntheMol
```

Install SyntheMol via pip.
```bash
pip install SyntheMol
```

Alternatively, clone the repo and install SyntheMol locally.
```bash
git clone https://github.com/swansonk14/SyntheMol.git
cd SyntheMol
pip install -e .
```

Download the necessary data files.
```bash
gdown "https://drive.google.com/drive/folders/1LLLwxe_nQAnsRSQpIRq_ngyCm1txS-Sq" -O /path/to/SyntheMol/files --folder
```

**Note:** Replace `/path/to/SyntheMol` with the path to the SyntheMol package. The path to SyntheMol can be found by running `python -c "import SyntheMol; print(SyntheMol.__path__)"`.

**Note:** If you get the issue `ImportError: libXrender.so.1: cannot open shared object file: No such file or directory`, run `conda install -c conda-forge xorg-libxrender`.


## Combinatorial chemical space

SyntheMol is currently designed to use 139,444 unique building blocks and 13 chemical reactions from the [Enamine REAL Space](https://enamine.net/compound-collections/real-compounds/real-space-navigator), which can produce over 30 billion molecules. However, an alternate combinatorial chemical space can optionally be used by replacing the building blocks and chemical reactions as follows.

TODO: check numbers with new building blocks

**Building blocks:** Replace `data/building_blocks.csv` with a custom file containing the building blocks. The file should be a CSV file with a header row and two columns: `smiles` and `ID`. The `smiles` column should contain the SMILES string for each building block, and the `ID` column should contain a unique ID for each building block.

**Chemical reactions:** In `SyntheMol/reactions/custom.py`, set `CUSTOM_REACTIONS` to a list of `Reaction` objects similar to the `REAL_REACTIONS` list in `SyntheMol/reactions/real.py`. If `CUSTOM_REACTIONS` is defined (i.e., not `None`), then it will automatically be used instead of `REAL_REACTIONS`.


## Bioactivity prediction model

SyntheMol requires a bioactivity prediction model to guide its generative process. SyntheMol is designed to use one of three types of models:

1. **Chemprop:** a message passing neural network from https://github.com/chemprop/chemprop
2. **Chemprop-RDKit:** Chemprop augmented with 200 RDKit molecular features
3. **Random forest:** a scikit-learn random forest model trained on 200 RDKit molecular features


### Train model

All three model types can be trained using [Chemprop](https://github.com/chemprop/chemprop), which is installed along with SyntheMol. All three model types can be trained on either regression or binary classification bioactivities. Full details are provided in the [Chemprop](https://github.com/chemprop/chemprop) README. Below is an example for training a Chemprop model on a binary classification task.

TODO: enable random forest RDKit in chemprop

TODO: enable regression models

Data file
```bash
# data/data.csv
smiles,activity
Br.CC(Cc1ccc(O)cc1)NCC(O)c1cc(O)cc(O)c1,0
CC[Hg]Sc1ccccc1C(=O)[O-].[Na+],1
O=C(O)CCc1ccc(NCc2cccc(Oc3ccccc3)c2)cc1,0
...
```

TODO: make sure training and predict commands work like this in chemprop

Train Chemprop
```bash
python -m chemprop.train \
    --data_path data/data.csv \
    --dataset_type classification \
    --save_dir models/chemprop
```


### Pre-compute building block scores

After training, use the model to pre-compute scores of building blocks to accelerate the SyntheMol generation process. Below is an example using the trained Chemprop model.

```bash
python -m chemprop.predict \
    --test_path data/building_blocks.csv \
    --preds_path models/chemprop/building_block_scores.csv \
    --checkpoint_dir models/chemprop
```


## Generate molecules

SyntheMol uses the bioactivity prediction model within a Monte Carlo tree search to generate molecules. Below is an example for generating molecules with a trained Chemprop model using 20,000 MCTS rollouts.

```bash
python -m SyntheMol.generate \
    --model_path models/chemprop \
    --model_type chemprop \
    --save_dir generations/chemprop \
    --building_blocks_path models/chemprop/building_block_scores.csv \
    --bulding_blocks_score_column activity \
    --n_rollout 20000
```


## Filter generated molecules

Optionally, the generated molecules can be filtered for structural novelty, predicted bioactivity, and structural diversity. These filtering steps use [chem_utils](https://github.com/swansonk14/chem_utils), which is installed along with SyntheMol. Below is an example for filtering the generated molecules.

### Novelty

Filter for novelty by comparing the generated molecules to a set of active molecules (hits) from the training set or literature and removing similar generated molecules.

Hits file
```bash
# data/hits.csv
smiles,activity
CC[Hg]Sc1ccccc1C(=O)[O-].[Na+],1
O=C(NNc1ccccc1)c1ccncc1,1
...
```

Compute Tversky similarity between generated molecules and hits
```bash
python -m chem_utils.nearest_neighbor \
    --data_path generations/chemprop/molecules.csv \
    --reference_data_path data/hits.csv \
    --reference_name hits \
    --metrics tversky
```

Filter by similarity, only keeping molecules with a nearest neighbor similarity to hits of at most 0.5
```bash
python -m chem_utils.filter_molecules \
    --data_path generations/chemprop/molecules.csv \
    --save_path generations/chemprop/molecules_novel.csv \
    --filter_column hits_tversky_nearest_neighbor_similarity \
    --max_value 0.5
```


### Bioactivity

Filter for predicted bioactivity by keeping the molecules with the top 20% highest predicted bioactivity.

```bash
python -m chem_utils.filter_molecules \
    --data_path generations/chemprop/molecules_novel.csv \
    --save_path generations/chemprop/molecules_novel_bioactive.csv \
    --filter_column score \
    --top_proportion 0.2
```


### Diversity

Filter for diversity by clustering molecules based on their Morgan fingerprint and only keeping the top scoring molecule from each cluster.

Cluster molecules into 50 clusters
```bash
python -m chem_utils.cluster_molecules \
    --data_path generations/chemprop/molecules_novel_bioactive.csv \
    --num_clusters 50
```

Select the top scoring molecule from each cluster
```bash
python -m chem_utils.select_from_clusters \
    --data_path generations/molecules_novel_bioactive.csv \
    --save_path generations/molecules_novel_bioactive_diverse.csv \
    --value_column score
```
