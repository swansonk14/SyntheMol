# Simulation Study

Steps to perform a simulation study of SyntheMol using a computational molecular property. This assumes that the data setup steps from README.md have already been performed.

TODO: convert to python scripts and adjust final steps to use both models and REAL


## Compute cLogP

Compute cLogP on the training set using [chem_utils](https://github.com/swansonk14/chem_utils).
```bash
python compute_properties.py \
    --data_path ../../combinatorial_antibiotics/data/antibiotics.csv \
    --properties clogp
```

For comparison purposes, compute cLogP on a random sample of REAL molecules using [chem_utils](https://github.com/swansonk14/chem_utils).
```bash
python compute_properties.py \
    --data_path ../../combinatorial_antibiotics/data/random_real.csv \
    --properties clogp
```


## Binarize clogP

Binarize cLogP by using 6.5 as a threshold.

Training data
```python
import pandas as pd

data = pd.read_csv('data/antibiotics.csv')
data['clogp_6.5'] = (data['clogp'] > 6.5).astype(int)
print(data['clogp_6.5'].value_counts())
data.to_csv('data/antibiotics.csv', index=False)

hits = data[data['logp_6.5'] == 1]
hits.to_csv('data/antibiotics_hits.csv', index=False)
```

495 positive (3.7%) out of 13,524.

Random sample of REAL molecules
```python
import pandas as pd

data = pd.read_csv('data/random_real.csv')
data['clogp_6.5'] = (data['clogp'] > 6.5).astype(int)
print(data['clogp_6.5'].value_counts())
data.to_csv('data/random_real.csv', index=False)
```

11 positive (0.044%) out of 25,000.


## Train model

Train a Chemprop model on binary cLogP data using 30 epochs (strong model) or 1 epoch (weak model).
```bash
#!/bin/bash

for EPOCHS in 30 1
do
python train_model.py \
    --data_path data/antibiotics.csv \
    --save_dir models/clogp_chemprop_${EPOCHS}_epochs \
    --dataset_type classification \
    --model_type chemprop \
    --property_column clogp_6.5 \
    --num_models 10 \
    --epochs ${EPOCHS}
done
```

30 epochs: ROC-AUC = 0.973 +/- 0.007, PRC-AUC = 0.736 +/- 0.049
1 epoch: ROC-AUC = 0.859 +/- 0.017, PRC-AUC = 0.198 +/- 0.045


## Map building blocks to model scores

Map building blocks to model scores.
```bash
#!/bin/bash

for EPOCHS in 30 1
do
python map_fragments_to_model_scores.py \
    --fragment_path data/building_blocks.csv \
    --model_path models/clogp_chemprop_${EPOCHS}_epochs \
    --save_path models/clogp_chemprop_${EPOCHS}_epochs/building_block_scores.json \
    --model_type chemprop
done
```


## Generate molecules

Apply SyntheMol to generate molecules.
```bash
#!/bin/bash

for EPOCHS in 30 1
do
python tree_search.py \
    --model_path models/clogp_chemprop_${EPOCHS}_epochs \
    --fragment_path data/building_blocks.csv \
    --reaction_to_reagents_path data/reaction_to_building_blocks.json \
    --fragment_to_model_score_path models/clogp_chemprop_${EPOCHS}_epochs/building_block_scores.json \
    --save_dir generations/clogp_chemprop_${EPOCHS}_epochs \
    --search_type mcts \
    --model_type chemprop \
    --n_rollout 20000 \
    --fragment_diversity \
    --max_reactions 1
done
```

1 epoch: 27,123 molecules in 7 hours and 19 minutes.
30 epochs: 25,550 generated molecules in 8 hours and 53 minutes.


## Compute true cLogP

Compute the true cLogP for generated molecules using [chem_utils](https://github.com/swansonk14/chem_utils).
```bash
#!/bin/bash

for EPOCHS in 30 1
do
python compute_properties.py \
    --data_path ../../combinatorial_antibiotics/generations/clogp_chemprop_${EPOCHS}_epochs/molecules.csv \
    --properties clogp
done
```

Plot distribution of train vs generated vs REAL cLogP using [chem_utils](https://github.com/swansonk14/chem_utils).
```bash
python plot_property_distribution.py \
    --data_paths ../../combinatorial_antibiotics/data/screening_data/AB_combined_with_properties.csv \
    ../../combinatorial_antibiotics/data/Enamine_REAL_space_sampled_25k_with_properties.csv \
    ../../combinatorial_antibiotics/generations/clogp_chemprop_30_epochs/molecules.csv \
    ../../combinatorial_antibiotics/generations/clogp_chemprop_1_epochs/molecules.csv \
    --property_column clogp \
    --save_dir ../../combinatorial_antibiotics/plots/clogp_generations \
    --min_value -7 \
    --max_value 12
```

Binarize the true cLogP for generated molecules using the 6.5 threshold.
```python
import pandas as pd

data = pd.read_csv('generations/clogp_chemprop_30_epochs/molecules.csv')
data['clogp_6.5'] = (data['clogp'] > 6.5).astype(int)
print(data['clogp_6.5'].value_counts())
data.to_csv('generations/clogp_chemprop_30_epochs/molecules.csv', index=False)
```

Print the percent of generated molecules with true cLogP > 6.5.
```bash
python -c "import pandas as pd
data = pd.read_csv('generations/clogp_chemprop_30_epochs/molecules.csv')
num_pos = sum(data['clogp_6.5'])
print(f'{num_pos:,} positive ({100 * num_pos / len(data):.2f}%) out of {len(data):,}')"
```

30 epochs: 15,693 positive (61.42%) out of 25,550 (vs 0.044% for random REAL).
1 epoch: 3,195 positive (11.78%) out of 27,123 (vs 0.044% for random REAL).
