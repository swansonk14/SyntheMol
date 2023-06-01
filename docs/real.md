# Enamine REAL Space Data Processing

This file contains instructions for processing the Enamine REAL Space data in preparation for use with SyntheMol. Here, we use the 2021 q3-4 version of the building blocks and the 2022 q1-2 version of the enumerated REAL Space molecules as described in our antibiotic generation paper [TODO](TODO). However, by default, SyntheMol now downloads the 2022 q1-2 version of the building blocks during installation.

The data referred to in this file can be downloaded from the Google Drive folder [here](https://drive.google.com/drive/folders/1VLPPUbY_FTKMjlXgRm09bPSSms206Dce?usp=share_link). Note that the instructions below assume that the relevant data is downloaded to the `data` directory.


## Reactions

The Enamine REAL Space consists of 169 reactions, which can be found in a PDF in the Google Drive folder. By default, SyntheMol uses 13 of these reactions, which are provided as SMARTS in `SyntheMol/reactions/real.py`.


## Building blocks

Below, we describe the steps for processing the building blocks.


### SDF to SMILES

The `building_blocks.sdf` file (from the Google Drive folder) is the 2021 q3-4 version containing 138,085 molecules.

Convert the building blocks from SDF to SMILES.
```bash
chemfunc sdf_to_smiles \
    --data_path data/4_real_space/building_blocks.sdf \
    --save_path data/4_real_space/building_blocks_raw.csv \
    --properties Reagent_ID Catalog_ID
```

All 138,085 molecules were successfully converted from SDF to SMILES, and among those 134,609 are unique molecules.

Note: The SMILES are likely not all unique because they do not include stereochemistry even though it is annotated with the `CFG` tag in the SDF file.


### Remove salts

Remove the salts from the building blocks. This will also canonicalize the SMILES using RDKit's canonicalization method.
```bash
chemfunc canonicalize_smiles \
    --data_path data/4_real_space/building_blocks_raw.csv \
    --save_path data/4_real_space/building_blocks.csv \
    --remove_salts \
    --delete_disconnected_mols
```

This removes 25 molecules whose salts cannot be stripped, leaving 138,060 molecules (unique IDs, duplicate SMILES) of which 132,479 are unique.

**Note:** This step is crucial to prevent errors in running reactions. Salts can cause reactions to create products that are the same as the reactants, leading to undesired infinite loops during molecule generation.

## REAL Space molecules

Below, we describe the steps for processing the REAL Space molecules.


### Download REAL Space

Download the building block and reaction IDs used for the full REAL Space of 31 billion compounds (2022 q1-2 version downloaded on August 30, 2022).
```bash
lftp -c "open -p 21 -u username,password ftp-rdb-fr.chem-space.com; mirror -c --parallel=16 . data/4_real_space/full_real"
```

In the above command, replace `username` and `password` with the appropriate values provided by Enamine.


### Map REAL reactions to building blocks

Determine which building blocks are valid in which REAL reactions.
```bash
python -m synthemol.data.map_real_reactions_to_building_blocks \
    --data_dir data/4_real_space/full_real \
    --save_path data/4_real_space/reaction_to_building_blocks.pkl
```

Total number of molecules = 31,507,987,117


### Count REAL reactions and building blocks

Determine which reactions and building blocks are most common in REAL space.
```bash
python -m synthemol.data.count_real_space \
    --data_dir data/4_real_space/full_real \
    --save_dir data/4_real_space
```


### Sample REAL molecules

Randomly sample 25,000 REAL Space molecules. This is used for analysis of a representative sample of REAL Space molecules.
```bash
python -m synthemol.data.sample_real_space \
    --data_dir data/4_real_space/full_real \
    --save_path data/4_real_space/random_real.csv \
    --num_molecules 25000
```


### Count feasible REAL molecules

Count feasible REAL Space molecules, i.e., those that can be produced when limited to the selected 13 reactions and the building blocks after processing.
```bash
python -m synthemol.data.count_real_space \
    --data_dir data/4_real_space/full_real \
    --save_dir data/4_real_space \
    --building_blocks_path data/4_real_space/building_blocks.csv \
    --only_selected_reactions
```

Number of files = 2,993
Number of building blocks = 138,060 (132,479 unique molecules)
Number of selected reactions = 13
Total number of molecules = 31,507,987,117
Total number of molecules with selected building blocks/reactions = 29,575,293,692
