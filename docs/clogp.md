# Generating High cLogP Molecules with SyntheMol

Instructions for performing an _in silico_ study of SyntheMol using a computational molecular property, cLogP, which is the computed octanol-water partition coefficient. Assumes relevant data has already been downloaded (see [docs/README.md](README.md)).

- [Compute cLogP](#compute-clogp)
- [Binarize clogP](#binarize-clogp)
- [Train model](#train-model)
- [Compute model scores for building blocks](#compute-model-scores-for-building-blocks)
- [Generate molecules](#generate-molecules)
- [Compute true cLogP](#compute-true-clogp)


## Compute cLogP

Compute cLogP on the training set.
```bash
chemfunc compute_properties \
    --data_path data/Data/1_training_data/antibiotics.csv \
    --properties clogp
```

For comparison purposes, compute cLogP on a random sample of REAL molecules.
```bash
chemfunc compute_properties \
    --data_path data/Data/4_real_space/random_real.csv \
    --properties clogp
```


## Binarize clogP

Binarize cLogP by using 6.5 as a threshold.

```bash
chemfunc regression_to_classification \
    --data_path data/Data/1_training_data/antibiotics.csv \
    --regression_column clogp \
    --threshold 6.5 \
    --classification_column clogp_6.5
```

**Antibiotics:** 495 cLogP positive (3.7%) out of 13,524.

```bash
chemfunc regression_to_classification \
    --data_path data/Data/4_real_space/random_real.csv \
    --regression_column clogp \
    --threshold 6.5 \
    --classification_column clogp_6.5
```

**Random REAL:** 11 cLogP positive (0.044%) out of 25,000.


## Train model

Train a Chemprop model on binary cLogP data using 30 epochs (strong model) or 1 epoch (weak model).
```bash
for EPOCHS in 1 30
do
python scripts/models/train.py \
    --data_path data/Data/1_training_data/antibiotics.csv \
    --save_dir data/Models/clogp_chemprop_${EPOCHS}_epochs \
    --dataset_type classification \
    --model_type chemprop \
    --property_column clogp_6.5 \
    --num_models 10 \
    --epochs ${EPOCHS}
done
```

**1 epoch:** ROC-AUC = 0.859 +/- 0.017, PRC-AUC = 0.198 +/- 0.045

**30 epochs:** ROC-AUC = 0.973 +/- 0.007, PRC-AUC = 0.736 +/- 0.049


## Compute model scores for building blocks

Compute model scores for building blocks.
```bash
for EPOCHS in 1 30
do
python scripts/models/predict.py \
    --data_path data/Data/4_real_space/building_blocks.csv \
    --save_path data/Models/clogp_chemprop_${EPOCHS}_epochs/building_blocks.csv \
    --model_path data/Models/clogp_chemprop_${EPOCHS}_epochs \
    --model_type chemprop \
    --average_preds
done
```


## Generate molecules

Apply SyntheMol to generate molecules.
```bash
for EPOCHS in 1 30
do
synthemol \
    --model_path data/Models/clogp_chemprop_${EPOCHS}_epochs \
    --model_type chemprop \
    --building_blocks_path data/Models/clogp_chemprop_${EPOCHS}_epochs/building_blocks.csv \
    --building_blocks_score_column chemprop_ensemble_preds \
    --building_blocks_id_column Reagent_ID \
    --reaction_to_building_blocks_path data/Data/4_real_space/reaction_to_building_blocks.pkl \
    --save_dir data/Data/10_generations_clogp/clogp_chemprop_${EPOCHS}_epochs \
    --max_reactions 1 \
    --n_rollout 20000 \
    --replicate
done
```

**1 epoch:** 27,123 molecules in 7 hours and 19 minutes.

**30 epochs:** 25,550 generated molecules in 8 hours and 53 minutes.


## Compute true cLogP

Compute the true cLogP for generated molecules.
```bash
for EPOCHS in 1 30
do
chemfunc compute_properties \
    --data_path data/Data/10_generations_clogp/clogp_chemprop_${EPOCHS}_epochs/molecules.csv \
    --properties clogp
done
```

Plot distribution of train vs generated vs REAL cLogP.
```bash
chemfunc plot_property_distribution \
    --data_paths data/Data/1_training_data/antibiotics.csv \
    data/Data/4_real_space/random_real.csv \
    data/Data/10_generations_clogp/clogp_chemprop_30_epochs/molecules.csv \
    data/Data/10_generations_clogp/clogp_chemprop_1_epochs/molecules.csv \
    --property_column clogp \
    --save_path plots/clogp_generations.pdf \
    --min_value -7 \
    --max_value 12
```

Binarize the true cLogP for generated molecules using the 6.5 threshold.

```bash
for EPOCHS in 1 30
do
chemfunc regression_to_classification \
    --data_path data/Data/10_generations_clogp/clogp_chemprop_${EPOCHS}_epochs/molecules.csv \
    --regression_column clogp \
    --threshold 6.5 \
    --classification_column clogp_6.5
done
```

**1 epoch:** 3,195 positive (11.78%) out of 27,123 (vs 0.044% for random REAL).

**30 epochs:** 15,693 positive (61.42%) out of 25,550 (vs 0.044% for random REAL).
