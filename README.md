# SyntheMol: Generative AI for Drug Discovery

[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/synthemol)](https://badge.fury.io/py/synthemol)
[![PyPI version](https://badge.fury.io/py/synthemol.svg)](https://badge.fury.io/py/synthemol)
[![Downloads](https://pepy.tech/badge/synthemol)](https://pepy.tech/project/synthemol)
[![license](https://img.shields.io/github/license/swansonk14/synthemol.svg)](https://github.com/swansonk14/SyntheMol/blob/main/LICENSE.txt)

SyntheMol is a generative AI method for designing structurally novel and diverse drug candidates with predicted
bioactivity that are easy to synthesize.

SyntheMol consists of either a reinforcement learning model (SyntheMol-RL) or a Monte Carlo tree search (SyntheMol-MCTS) to explore a combinatorial chemical space consisting of
molecular building blocks and chemical reactions. In both cases, SyntheMol is guided by a bioactivity prediction AI model, such as a
graph neural network or multilayer perceptron. SyntheMol uses two chemical spaces which jointly contain over 46 billion molecules:

[Enamine REAL Space](https://enamine.net/compound-collections/real-compounds/real-space-navigator)
* 137,656 molecular building blocks
* 70 chemical reactions
* 30,330,025,259 molecules

[WuXi GalaXi](https://www.biosolveit.de/wp-content/uploads/2021/08/Xu.pdf)
* 14,977 molecular building blocks
* 36 chemical reactions
* 16,146,071,436 molecules

Notably, SyntheMol can be easily adapted to use any set of building blocks and reactions.

SyntheMol-RL will be described in a forthcoming paper.

SyntheMol-MCTS is described in the following paper, where we applied SyntheMol to design novel antibiotic candidates for the Gram-negative bacterium _Acinetobacter baumannii_.

Swanson, K., Liu, G., Catacutan, D. B., Arnold, A., Zou, J., Stokes, J. M. [Generative AI for designing and validating easily synthesizable and structurally novel antibiotics](https://www.nature.com/articles/s42256-024-00809-7). _Nature Machine Intelligence_, 2024.

Full details for reproducing the SyntheMol results in both papers are provided in the [docs](docs) directory.

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

SyntheMol can be installed in < 3 minutes on any operating system using pip (optionally within a conda environment).
SyntheMol can be run on a standard laptop (e.g., 16 GB memory and 8-16 CPUs), although a GPU is useful for faster
training and prediction of the underlying bioactivity prediction model (Chemprop).

Optionally, create a conda environment.

```bash
conda create -y -n synthemol python=3.12
conda activate synthemol
```

Install SyntheMol via pip.

```bash
pip install synthemol
```

Alternatively, clone the repo and install SyntheMol locally.

```bash
git clone https://github.com/swansonk14/SyntheMol.git
cd SyntheMol
pip install -e .
```

If there are version issues with the required packages, create a conda environment with specific working versions of the
packages as follows.

SyntheMol-RL

```bash
pip install -r requirements_rl.txt
pip install -e .
```

SyntheMol-MCTS

```bash
pip install -r requirements_mcts.txt
pip install -e .
```

**Note:** If you get the
issue `ImportError: libXrender.so.1: cannot open shared object file: No such file or directory`,
run `conda install -c conda-forge xorg-libxrender`.

## Combinatorial chemical space

An alternate combinatorial chemical space can optionally be used by replacing the building blocks and chemical reactions as follows.

**Building blocks:** Create a building blocks file similar to `synthemol/resources/real/building_blocks.csv` with a custom file containing the building blocks. The file
should be a CSV file with a header row and two columns: `smiles` and `reagent_id`. The `smiles` column should contain the SMILES
string for each building block, and the `reagent_id` column should contain a unique ID for each building block.

**Chemical reactions:** In `SyntheMol/reactions/custom.py`, set `CUSTOM_REACTIONS` to a list of `Reaction` objects
similar to the `REAL_REACTIONS` list in `SyntheMol/reactions/real.py`. Then specify `--chemical_spaces custom` when running SyntheMol.

## Bioactivity prediction model

SyntheMol requires a bioactivity prediction model to guide its generative process. SyntheMol is designed to use one of
these types of models:

1. **Chemprop:** a message passing neural network from https://github.com/chemprop/chemprop
2. **Chemprop-RDKit:** Chemprop augmented with 200 RDKit molecular features
3. **MLP-RDKit:** a feed-forward neural network using 200 RDKit molecular features
4. **Random forest:** a scikit-learn random forest model trained on 200 RDKit molecular features (SyntheMol-MCTS only)

### Train model

All model types can be trained using [Chemprop](https://github.com/chemprop/chemprop), which is installed along
with SyntheMol. All model types can be trained on either regression or binary classification bioactivities. Full
details are provided in the [Chemprop](https://github.com/chemprop/chemprop) README. Below is an example for training a
Chemprop model on a binary classification task. By default, training is done on a GPU (if available).

Data file

```bash
# data/data.csv
smiles,activity
Br.CC(Cc1ccc(O)cc1)NCC(O)c1cc(O)cc(O)c1,0
CC[Hg]Sc1ccccc1C(=O)[O-].[Na+],1
O=C(O)CCc1ccc(NCc2cccc(Oc3ccccc3)c2)cc1,0
...
```

Train Chemprop

```bash
chemprop_train \
    --data_path data/data.csv \
    --dataset_type classification \
    --save_dir models/chemprop
```

### Pre-compute building block scores

After training, use the model to pre-compute scores of building blocks to accelerate the SyntheMol generation process.
Below is an example using the trained Chemprop model. By default, prediction is done on a GPU (if available).

```bash
chemprop_predict \
    --test_path "$(python -c 'import synthemol; print(str(synthemol.constants.BUILDING_BLOCKS_PATH))')" \
    --preds_path models/chemprop/building_blocks.csv \
    --checkpoint_dir models/chemprop
```

## Generate molecules

SyntheMol uses the bioactivity prediction model a generative model (RL or MCTS) to generate molecules. SyntheMol-RL can be powered either by a Chemprop-based RL model (RL-Chemprop) or an MLP-based RL model (RL-MLP). Below is an example using SyntheMol-RL (RL-Chemprop version) to generate molecules with a trained Chemprop model for 10,000 rollouts using the Enamine REAL Space (the default).

```bash
synthemol \
    --score_model_paths models/chemprop \
    --score_types chemprop \
    --chemical_spaces real wuxi \
    --building_blocks_score_columns activity \
    --save_dir generations/chemprop \
    --n_rollout 10000 \
    --search_type rl \
    --rl_model_type chemprop \
    --rl_model_fingerprint_type rdkit \
    --rl_model_paths models/chemprop/fold_0/model_0/model.pt \
    --rl_prediction_types classification
```

Note: The `building_blocks_score_columns` must match the column name in the building blocks file that contains the
building block scores. When using `chemprop_train` and `chemprop_predict`, the column name will be the same as the
column that contains target activity/property values in the training data file (e.g., `activity`).

For more complex examples with multiparameter optimization using both the RL-Chemprop and RL-MLP versions of SyntheMol-RL, please see [docs/rl/antibiotics.md](docs/rl/antibiotics.md).

## Filter generated molecules

Optionally, the generated molecules can be filtered for structural novelty, predicted bioactivity, and structural
diversity. These filtering steps use [chemfunc](https://github.com/swansonk14/chemfunc), which is installed along with
SyntheMol. Below is an example for filtering the generated molecules.

### Novelty

Filter for novelty by comparing the generated molecules to a set of active molecules (hits) from the training set or
literature and removing similar generated molecules.

Hits file

```bash
# data/hits.csv
smiles,activity
CC[Hg]Sc1ccccc1C(=O)[O-].[Na+],1
O=C(NNc1ccccc1)c1ccncc1,1
...
```

Compute Tversky similarity between generated molecules and hits.

```bash
chemfunc nearest_neighbor \
    --data_path generations/chemprop/molecules.csv \
    --reference_data_path data/hits.csv \
    --reference_name hits \
    --metric tversky
```

Filter by similarity, only keeping molecules with a nearest neighbor similarity to hits of at most 0.5.

```bash
chemfunc filter_molecules \
    --data_path generations/chemprop/molecules.csv \
    --save_path generations/chemprop/molecules_novel.csv \
    --filter_column hits_tversky_nearest_neighbor_similarity \
    --max_value 0.5
```

### Bioactivity

Filter for predicted bioactivity by keeping the molecules with the top 20% highest predicted bioactivity.

```bash
chemfunc filter_molecules \
    --data_path generations/chemprop/molecules_novel.csv \
    --save_path generations/chemprop/molecules_novel_bioactive.csv \
    --filter_column score \
    --top_proportion 0.2
```

### Diversity

Filter for diversity by clustering molecules based on their Morgan fingerprint and only keeping the top scoring molecule
from each cluster.

Cluster molecules into 50 clusters.

```bash
chemfunc cluster_molecules \
    --data_path generations/chemprop/molecules_novel_bioactive.csv \
    --num_clusters 50
```

Select the top scoring molecule from each cluster.

```bash
chemfunc select_from_clusters \
    --data_path generations/chemprop/molecules_novel_bioactive.csv \
    --save_path generations/chemprop/molecules_novel_bioactive_diverse.csv \
    --value_column score
```
