# Enamine REAL Space Data Processing

Instructions for processing Enamine REAL building blocks, reactions, and molecules. Assumes relevant data has already been downloaded (see [docs/README.md](README.md)).


- [Reactions](#reactions)
- [Building blocks](#building-blocks)
  * [SDF to SMILES](#sdf-to-smiles)
  * [Remove salts](#remove-salts)
- [REAL Space molecules](#real-space-molecules)
  * [Download REAL Space](#download-real-space)
  * [Map REAL reactions to building blocks](#map-real-reactions-to-building-blocks)
  * [Count REAL reactions and building blocks](#count-real-reactions-and-building-blocks)
  * [Sample REAL molecules](#sample-real-molecules)
  * [Count feasible REAL molecules](#count-feasible-real-molecules)


## Reactions

The Enamine REAL Space consists of 169 reactions, which can be found in a PDF in the Google Drive folder. By default, SyntheMol uses 13 of these reactions, which are provided as SMARTS in `synthemol/reactions/real.py`.


## Building blocks

Below, we describe the steps for processing the building blocks.


### SDF to SMILES

The `data/Data/4_real_space/building_blocks.sdf` file is the 2021 q3-4 version containing 138,085 molecules.

Convert the building blocks from SDF to SMILES.
```bash
chemfunc sdf_to_smiles \
    --data_path data/Data/4_real_space/building_blocks_raw.sdf \
    --save_path data/Data/4_real_space/building_blocks_raw.csv \
    --properties Reagent_ID Catalog_ID
```

All 138,085 molecules were successfully converted from SDF to SMILES, and among those 134,609 are unique molecules.

Note: The SMILES are likely not all unique because they do not include stereochemistry even though it is annotated with the `CFG` tag in the SDF file.


### Remove salts

Remove the salts from the building blocks. This will also canonicalize the SMILES using RDKit's canonicalization method.
```bash
chemfunc canonicalize_smiles \
    --data_path data/Data/4_real_space/building_blocks_raw.csv \
    --save_path data/Data/4_real_space/building_blocks.csv \
    --remove_salts \
    --delete_disconnected_mols
```

This removes 25 molecules whose salts cannot be stripped, leaving 138,060 molecules (unique IDs, duplicate SMILES) of which 132,479 are unique.

**Note:** This step is crucial to prevent errors in running reactions. Salts can cause reactions to create products that are the same as the reactants, leading to undesired infinite loops during molecule generation.

## REAL Space molecules

Below, we describe the steps for processing the REAL Space molecules.


### Download REAL Space

Download the building block and reaction IDs used for the full REAL Space of 31 billion compounds (we used the 2022 q1-2 version downloaded on August 30, 2022).
```bash
lftp -c "open -p 21 -u username,password ftp-rdb-fr.chem-space.com; mirror -c --parallel=16 . data/Data/4_real_space/full_real"
```

In the above command, replace `username` and `password` with the appropriate values provided by Enamine.


### Map REAL reactions to building blocks

Determine which building blocks are valid in which REAL reactions.
```bash
python scripts/data/map_real_reactions_to_building_blocks.py \
    --data_dir data/Data/4_real_space/full_real \
    --save_path data/Data/4_real_space/reaction_to_building_blocks.pkl
```

Total number of molecules = 31,507,987,117


### Count REAL reactions and building blocks

Determine which reactions and building blocks are most common in REAL space.
```bash
python scripts/data/count_real_space.py \
    --data_dir data/Data/4_real_space/full_real \
    --save_dir data/Data/4_real_space
```


### Sample REAL molecules

Randomly sample 25,000 REAL Space molecules. This is used for analysis of a representative sample of REAL Space molecules.
```bash
python scripts/data/sample_real_space.py \
    --data_dir data/Data/4_real_space/full_real \
    --save_path data/Data/4_real_space/random_real.csv \
    --num_molecules 25000
```


### Count feasible REAL molecules

Count feasible REAL Space molecules, i.e., those that can be produced when limited to the selected 13 reactions and the building blocks after processing.
```bash
python scripts/data/count_real_space.py \
    --data_dir data/Data/4_real_space/full_real \
    --save_dir data/Data/4_real_space \
    --building_blocks_path data/Data/4_real_space/building_blocks.csv \
    --only_selected_reactions
```

Number of files = 2,993
Number of building blocks = 138,060 (132,479 unique molecules)
Number of selected reactions = 13
Total number of molecules = 31,507,987,117
Total number of molecules with selected building blocks/reactions = 29,575,293,692
