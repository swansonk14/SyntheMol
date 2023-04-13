# Simulation Study

This file contains instructions for performing a simulation study of SyntheMol using a computational molecular property, cLogP (the computed octanol-water partition coefficient).

The data and models referred to in this file can be downloaded from the Google Drive folder [here](https://drive.google.com/drive/folders/1VLPPUbY_FTKMjlXgRm09bPSSms206Dce?usp=share_link). Note that the instructions below assume that the relevant data is downloaded to the `data` directory.

TODO: potentially change order of cLogP data in Google Drive and add to supplementary table

* [Compute cLogP](#compute-clogp)
* [Binarize clogP](#binarize-clogp)
* [Train model](#train-model)
* [Map building blocks to model scores](#map-building-blocks-to-model-scores)
* [Generate molecules](#generate-molecules)
* [Compute true cLogP](#compute-true-clogp)


## Compute cLogP

Compute cLogP on the training set.
```bash
python -m chem_utils.compute_properties \
    --data_path data/1_training_data/antibiotics.csv \
    --properties clogp
```

For comparison purposes, compute cLogP on a random sample of REAL molecules.
```bash
python -m chem_utils.compute_properties \
    --data_path data/4_real_space/random_real.csv \
    --properties clogp
```


## Binarize clogP

Binarize cLogP by using 6.5 as a threshold.

```bash
for DATA_NAME in antibiotics random_real
do
python -m chem_utils.regression_to_classification \
    --data_path data/1_training_data/${DATA_NAME}.csv \
    --regression_column clogp \
    --threshold 6.5 \
    --classification_column clogp_6.5
done
```

**Antibiotics:** 495 cLogP positive (3.7%) out of 13,524.

**Random REAL:** 11 cLogP positive (0.044%) out of 25,000.


## Train model

Train a Chemprop model on binary cLogP data using 30 epochs (strong model) or 1 epoch (weak model).
```bash
#!/bin/bash

for EPOCHS in 30 1
do
python -m SyntheMol.models.train \
    --data_path data/1_training_data/antibiotics.csv \
    --save_dir models/clogp_chemprop_${EPOCHS}_epochs \
    --dataset_type classification \
    --model_type chemprop \
    --property_column clogp_6.5 \
    --num_models 10 \
    --epochs ${EPOCHS}
done
```

**30 epochs:** ROC-AUC = 0.973 +/- 0.007, PRC-AUC = 0.736 +/- 0.049

**1 epoch:** ROC-AUC = 0.859 +/- 0.017, PRC-AUC = 0.198 +/- 0.045


## Map building blocks to model scores

Map building blocks to model scores.
```bash
#!/bin/bash

for EPOCHS in 30 1
do
python -m SyntheMol.models.predict \
    --data_path data/4_real_space/building_blocks.csv \
    --model_path models/clogp_chemprop_${EPOCHS}_epochs \
    --model_type chemprop \
    --average_preds
done
```


## Generate molecules

Apply SyntheMol to generate molecules.
```bash
#!/bin/bash

for EPOCHS in 30 1
do
python -m SyntheMol.generate \
    --model_path models/clogp_chemprop_${EPOCHS}_epochs \
    --model_type chemprop \
    --building_blocks_path data/4_real_space/building_blocks.csv \
    --save_dir data/10_generations_clogp/clogp_chemprop_${EPOCHS}_epochs \
    --reaction_to_building_blocks_path data/4_real_space/reaction_to_building_blocks.pkl \
    --max_reactions 1 \
    --n_rollout 20000
done
```

**30 epochs:** 25,550 generated molecules in 8 hours and 53 minutes.

**1 epoch:** 27,123 molecules in 7 hours and 19 minutes.


## Compute true cLogP

Compute the true cLogP for generated molecules.
```bash
#!/bin/bash

for EPOCHS in 30 1
do
python -m chem_utils.compute_properties \
    --data_path data/10_generations_clogp/clogp_chemprop_${EPOCHS}_epochs/molecules.csv \
    --properties clogp
done
```

Plot distribution of train vs generated vs REAL cLogP.
```bash
python -m chem_utils.plot_property_distribution \
    --data_paths data/1_training_data/antibiotics.csv \
    data/4_real_space/random_real.csv \
    data/10_generations_clogp/clogp_chemprop_30_epochs/molecules.csv \
    data/10_generations_clogp/clogp_chemprop_1_epochs/molecules.csv \
    --property_column clogp \
    --save_dir plots/clogp_generations \
    --min_value -7 \
    --max_value 12
```

Binarize the true cLogP for generated molecules using the 6.5 threshold.

```bash
for EPOCHS in 30 1
do
python -m chem_utils.regression_to_classification \
    --data_path data/10_generations_clogp/clogp_chemprop_${EPOCHS}_epochs/molecules.csv \
    --regression_column clogp \
    --threshold 6.5 \
    --classification_column clogp_6.5
done
```

**30 epochs:** 15,693 positive (61.42%) out of 25,550 (vs 0.044% for random REAL).

**1 epoch:** 3,195 positive (11.78%) out of 27,123 (vs 0.044% for random REAL).
