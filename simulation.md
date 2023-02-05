# Simulation Study

Steps to perform a simulation study of the generative model using a computational molecular property. This assumes that the data setup steps from README.md have already been performed.

TODO: change everything from logP to cLogP


## Compute logP

Compute logP on the training set using [chem_utils](https://github.com/swansonk14/chem_utils).

```bash
python compute_properties.py \
    --data_path ../../combinatorial_antibiotics/data/screening_data/AB_combined.csv \
    --properties logp \
    --save_path ../../combinatorial_antibiotics/data/screening_data/AB_combined_with_properties.csv
```

For comparison purposes, compute logP on a random sample of REAL molecules using [chem_utils](https://github.com/swansonk14/chem_utils).

```bash
python compute_properties.py \
    --data_path ../../combinatorial_antibiotics/data/enamine_REAL_space_sampled_25k.csv \
    --properties logp \
    --save_path ../../combinatorial_antibiotics/data/enamine_REAL_space_sampled_25k_with_properties.csv
```


## Binarize logP

Binarize logP by using 6.5 as a threshold.

Training data.

```python
import pandas as pd

data = pd.read_csv('data/screening_data/AB_combined_with_properties.csv')
data['logp_6.5'] = (data['logp'] > 6.5).astype(int)
print(data['logp_6.5'].value_counts())
data.to_csv('data/screening_data/AB_combined_with_properties.csv', index=False)

hits = data[data['logp_6.5'] == 1]
hits.to_csv('data/screening_data/AB_combined_logp_6.5_hits.csv', index=False)
```

495 positive (3.7%) out of 13,524.

Random sample of REAL molecules.

```python
import pandas as pd

data = pd.read_csv('data/Enamine_REAL_space_sampled_25k_with_properties.csv')
data['logp_6.5'] = (data['logp'] > 6.5).astype(int)
print(data['logp_6.5'].value_counts())
data.to_csv('data/Enamine_REAL_space_sampled_25k_with_properties.csv', index=False)
```

11 positive (0.044%) out of 25,000.


## Train model

Train model on binary logP data.

```bash
python train_model.py \
    --data_path data/screening_data/AB_combined_with_properties.csv \
    --save_dir ckpt/logp_6.5_chemprop \
    --dataset_type classification \
    --model_type chemprop \
    --property_column logp_6.5 \
    --num_models 10 \
    --epochs 1
```

1 epoch: ROC-AUC = 0.859 +/- 0.017, PRC-AUC = 0.198 +/- 0.045
30 epochs: ROC-AUC = 0.973 +/- 0.007, PRC-AUC = 0.736 +/- 0.049


## Map fragments to model scores

Map fragments to model scores.

```bash
python map_fragments_to_model_scores.py \
    --fragment_path data/2021q3-4_Enamine_REAL_reagents_SMILES_no_salts.csv \
    --model_path ckpt/logp_6.5_chemprop \
    --save_path ckpt/logp_6.5_chemprop/fragments_to_model_scores.json \
    --model_type chemprop
```


## Run MCTS

Tree search with MCTS.

```bash
python tree_search.py \
    --model_path ckpt/logp_6.5_chemprop \
    --fragment_path data/2021q3-4_Enamine_REAL_reagents_SMILES_no_salts.csv \
    --reaction_to_reagents_path data/reaction_to_reagents_REAL_space.json \
    --fragment_to_model_score_path ckpt/logp_6.5_chemprop/fragments_to_model_scores.json \
    --save_dir generations/logp_6.5_chemprop \
    --search_type mcts \
    --model_type chemprop \
    --n_rollout 20000 \
    --fragment_diversity \
    --max_reactions 1
```

1 epoch: 27,123 molecules in 7 hours and 19 minutes.
30 epochs: 25,550 generated molecules in 8 hours and 53 minutes.


## Compute true logP

Compute true logP for generated molecules using [chem_utils](https://github.com/swansonk14/chem_utils).

```bash
python compute_properties.py \
    --data_path ../../combinatorial_antibiotics/generations/logp_6.5_chemprop/molecules.csv \
    --properties logp
```

Plot distribution of train vs generated vs REAL logP using [chem_utils](https://github.com/swansonk14/chem_utils).
```bash
python plot_property_distribution.py \
    --data_paths ../../combinatorial_antibiotics/data/screening_data/AB_combined_with_properties.csv \
    ../../combinatorial_antibiotics/data/Enamine_REAL_space_sampled_25k_with_properties.csv \
    ../../combinatorial_antibiotics/generations/logp_6.5_chemprop/molecules.csv \
    --property_column logp \
    --save_dir ../../combinatorial_antibiotics/generations/logp_6.5_chemprop \
    --min_value -7 \
    --max_value 12
```

Binarize true logP for generated molecules

```python
import pandas as pd

data = pd.read_csv('generations/logp_6.5_chemprop/molecules.csv')
data['logp_6.5'] = (data['logp'] > 6.5).astype(int)
print(data['logp_6.5'].value_counts())
data.to_csv('generations/logp_6.5_chemprop/molecules.csv', index=False)
```

Print percent of generated molecules with true logP > 6.5.

```bash
python -c "import pandas as pd
data = pd.read_csv('generations/logp_6.5_chemprop/molecules.csv')
num_pos = sum(data['logp_6.5'])
print(f'{num_pos:,} positive ({100 * num_pos / len(data):.2f}%) out of {len(data):,}')"
```

1 epoch: 3,195 positive (11.78%) out of 27,123 (vs 0.044% for random REAL).
30 epochs: 15,693 positive (61.42%) out of 25,550 (vs 0.044% for random REAL).

Compute metrics.

```bash
python -c "import pandas as pd
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score, average_precision_score
data = pd.read_csv('generations/logp_6.5_chemprop/molecules.csv')
print(spearmanr(data['logp'], data['score']))
print('ROC-AUC: ' + str(roc_auc_score(data['logp_6.5'], data['score'])))
print('PRC-AUC: ' + str(average_precision_score(data['logp_6.5'], data['score'])))"
```

1 epoch:
SpearmanrResult(correlation=0.42308232369949134, pvalue=0.0)
ROC-AUC: 0.6721572712399064
PRC-AUC: 0.18060276214822235 (Note: high b/c 11.8% are positive)

30 epochs:
SpearmanrResult(correlation=0.10148604738350796, pvalue=1.8034441138216467e-59)
ROC-AUC: 0.540647770477802
PRC-AUC: 0.62236162848842 (Note: high b/c 61.4% are positive)

Note: Not very predictive b/c these molecules were selected for high score so model already thinks all are positive rather than thinking some are negative.


## Assess generated molecules

Assess generated molecules

```bash
python assess_generated_molecules.py \
    --data_path generations/logp_6.5_chemprop/molecules.csv \
    --save_dir generations/logp_6.5_chemprop \
    --train_hits_path data/screening_data/AB_combined_logp_6.5_hits.csv
```


## Filter generated molecules

Compute nearest neighbor Tversky similarities to training hits using [chem_utils](https://github.com/swansonk14/chem_utils).

```bash
python nearest_neighbor.py \
    --data_path ../../combinatorial_antibiotics/generations/logp_6.5_chemprop/molecules.csv \
    --reference_data_path ../../combinatorial_antibiotics/data/screening_data/AB_combined_logp_6.5_hits.csv \
    --reference_name train_hits \
    --metrics tversky
```

Filter generated molecules based on nearest neighbor Tversky similarities to training hits using [chem_utils](https://github.com/swansonk14/chem_utils).

```bash
python filter_molecules.py \
    --data_path ../../combinatorial_antibiotics/generations/logp_6.5_chemprop/molecules.csv \
    --save_path ../../combinatorial_antibiotics/generations/logp_6.5_chemprop/molecules_train_sim_below_0.5.csv \
    --filter_column train_hits_tversky_nearest_neighbor_similarity \
    --max_value 0.5
```

```bash
python -c "import pandas as pd
data = pd.read_csv('generations/logp_6.5_chemprop/molecules_train_sim_below_0.5.csv')
num_pos = sum(data['logp_6.5'])
print(f'{num_pos:,} positive ({100 * num_pos / len(data):.2f}%) out of {len(data):,}')"
```

1 epoch: 2,427 positive (10.13%) out of 23,953
30 epochs: 9,736 positive (54.69%) out of 17,803

Filter generated molecules based on logP prediction score using [chem_utils](https://github.com/swansonk14/chem_utils).

```bash
python filter_molecules.py \
    --data_path ../../combinatorial_antibiotics/generations/logp_6.5_chemprop/molecules_train_sim_below_0.5.csv \
    --save_path ../../combinatorial_antibiotics/generations/logp_6.5_chemprop/molecules_train_sim_below_0.5_top_20_percent.csv \
    --filter_column score \
    --top_proportion 0.2
```

```bash
python -c "import pandas as pd
data = pd.read_csv('generations/logp_6.5_chemprop/molecules_train_sim_below_0.5_top_20_percent.csv')
num_pos = sum(data['logp_6.5'])
print(f'{num_pos:,} positive ({100 * num_pos / len(data):.2f}%) out of {len(data):,}')"
```

1 epoch: 782 positive (16.32%) out of 4,791
30 epochs: 1,928 positive (54.14%) out of 3,561

Cluster generated molecules based on Tanimoto similarity using [chem_utils](https://github.com/swansonk14/chem_utils).

```bash
python cluster_molecules.py \
    --data_path ../../combinatorial_antibiotics/generations/logp_6.5_chemprop/molecules_train_sim_below_0.5_top_20_percent.csv \
    --num_clusters 50
```

Select molecules from clusters using [chem_utils](https://github.com/swansonk14/chem_utils).

```bash
python select_from_clusters.py \
    --data_path ../../combinatorial_antibiotics/generations/logp_6.5_chemprop/molecules_train_sim_below_0.5_top_20_percent.csv \
    --save_path ../../combinatorial_antibiotics/generations/logp_6.5_chemprop/molecules_train_sim_below_0.5_top_20_percent_selected_50.csv \
    --value_column score
```

```bash
python -c "import pandas as pd
data = pd.read_csv('generations/logp_6.5_chemprop/molecules_train_sim_below_0.5_top_20_percent_selected_50.csv')
num_pos = sum(data['logp_6.5'])
print(f'{num_pos:,} positive ({100 * num_pos / len(data):.2f}%) out of {len(data):,}')"
```

1 epoch: 6 positive (12.00%) out of 50
30 epochs: 21 positive (42.00%) out of 50


Visualize selected molecules using [chem_utils](https://github.com/swansonk14/chem_utils).

```bash
python visualize_molecules.py \
    --data_path ../../combinatorial_antibiotics/generations/logp_6.5_chemprop/molecules_train_sim_below_0.5_top_20_percent_selected_50.csv \
    --save_dir ../../combinatorial_antibiotics/generations/logp_6.5_chemprop
cd ../../combinatorial_antibiotics/generations/logp_6.5_chemprop
mv 1.png molecules_train_sim_below_0.5_top_20_percent_selected_50.png
```

Assess selected molecules.

```bash
python assess_generated_molecules.py \
    --data_path generations/logp_6.5_chemprop/molecules_train_sim_below_0.5_top_20_percent_selected_50.csv \
    --save_dir generations/logp_6.5_chemprop/analysis_molecules_train_sim_below_0.5_top_20_percent_selected_50 \
    --train_hits_path data/screening_data/AB_combined_logp_6.5_hits.csv
```
