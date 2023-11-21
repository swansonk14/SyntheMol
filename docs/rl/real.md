# Enamine REAL Space Data Processing (SyntheMol-RL)

Instructions for processing the Enamine REAL space reactions, building blocks, and molecules using the 2022 REAL release with the 2022 building blocks (SyntheMol-RL paper). These instructions assume that relevant data has already been downloaded (see [docs/README.md](README.md)).


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

The `rl/data/real/building_blocks.sdf` file is the 2022 q1-2 version containing 139,517 molecules.

Convert the building blocks from SDF to SMILES.
```bash
chemfunc sdf_to_smiles \
    --data_path rl/data/real/building_blocks_raw.sdf \
    --save_path rl/data/real/building_blocks_raw.csv \
    --properties reagent_id
```

All 139,517 molecules were successfully converted from SDF to SMILES, and among those 139,444 are unique molecules.


### Remove salts

Remove the salts from the building blocks. This will also canonicalize the SMILES using RDKit's canonicalization method.
```bash
chemfunc canonicalize_smiles \
    --data_path rl/data/real/building_blocks_raw.csv \
    --save_path rl/data/real/building_blocks.csv \
    --remove_salts \
    --delete_disconnected_mols
```

This removes 24 molecules whose salts cannot be stripped, leaving 139,493 molecules (unique IDs, duplicate SMILES) of which 137,656 are unique.

**Note:** This step is crucial to prevent errors in running reactions. Salts can cause reactions to create products that are the same as the reactants, leading to undesired infinite loops during molecule generation.


## REAL Space molecules

Below, we describe the steps for processing the REAL Space molecules.

Note: The enumerated REAL space requires over 1T of storage and a similar amount of RAM in order to process it. However, the following steps are a one-time operation, and only minimal storage and RAM are required to store the processed data and run SyntheMol. Each of the following commands run in <= 24 hours using 10 cores and <= 1 TB of RAM.


### Download REAL Space

Download the building block and reaction IDs used for the full REAL Space of 31,507,987,117 billion compounds (we used the 2022 q1-2 version downloaded on August 30, 2022).
```bash
lftp -c "open -p 21 -u username,password ftp-rdb-fr.chem-space.com; mirror -c --parallel=16 . data/Data/4_real_space/full_real"
```

In the above command, replace `username` and `password` with the appropriate values provided by Enamine.

**Note:** Newer versions of the REAL space can be downloaded using the command below. However, the building block IDs in newer releases may not be compatible with the 2022 q1-2 building blocks used here.

```bash
mkdir rl/data/real/real_space
sftp -oPort=22 username@ftp-rdb-fr.chem-space.com
lcd rl/data/real/real_space
get -r .
```

In the above command, replace `username` with the appropriate value provided by Enamine, and type in the password provided by Enamine when prompted after the first command.


### Map REAL reactions to building blocks

Determine which building blocks are valid in which REAL reactions.
```bash
python scripts/data/map_real_reactions_to_building_blocks.py \
    --data_dir rl/data/real/real_space \
    --save_path rl/data/real/reaction_to_building_blocks_raw.pkl
```

Total number of molecules = 31,507,987,117

Filter out building blocks that do not match the reaction templates and map building block IDs to SMILES.
```bash
python scripts/data/filter_real_reactions_to_building_blocks.py \
    --reaction_to_building_blocks_path rl/data/real/reaction_to_building_blocks_raw.pkl \
    --building_blocks_path rl/data/real/building_blocks.csv \
    --save_path rl/data/real/reaction_to_building_blocks.pkl
```


### Count REAL reactions and building blocks

Determine which reactions and building blocks are most common in REAL space.
```bash
python scripts/data/count_real_space.py \
    --data_dir rl/data/real/real_space \
    --save_dir rl/data/real
```

Number of files = 2,993
Total number of molecules = 31,507,987,117
Total number of molecules with selected building blocks/reactions = 31,507,987,117


### Count feasible REAL molecules

Count feasible REAL Space molecules, i.e., those that can be produced when limited to the selected 13 reactions and the building blocks after processing.
```bash
python scripts/data/count_real_space.py \
    --data_dir rl/data/real/real_space \
    --save_dir rl/data/real \
    --building_blocks_path rl/data/real/building_blocks.csv \
    --only_selected_reactions
```

Number of files = 2,993
Number of building blocks = 139,493
Number of selected reactions = 13
Total number of molecules = 31,507,987,117
Total number of molecules with selected building blocks/reactions = 30,330,025,259


### Sample REAL molecules

Randomly sample 25,000 REAL Space molecules. This is used for analysis of a representative sample of REAL Space molecules.
```bash
python scripts/data/sample_real_space.py \
    --data_dir rl/data/real/real_space \
    --save_path rl/data/real/random_real_25k.csv \
    --num_molecules 25000
```

Randomly sample 25 million REAL Space molecules. This is used for a time-based comparison of Chemprop versus SyntheMol-RL.
```bash
python scripts/data/sample_real_space.py \
    --data_dir rl/data/real/real_space \
    --save_path rl/data/real/random_real_25m.csv \
    --num_molecules 25000000
```
