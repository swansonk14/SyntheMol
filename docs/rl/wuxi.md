# WuXi GalaXi Data Processing

Instructions for processing WuXi GalaXi REAL reactions and building blocks. Assumes relevant data has already been downloaded (see [docs/README.md](README.md)).


## Download

Download the WuXi GalaXi data with the following command.
```bash
gdown "https://drive.google.com/drive/folders/1BCM5P-tDjLFLkyiEJX9p4-6fw0w5-qRg?usp=drive_link" -O $(python -c "import synthemol; from pathlib import Path; print(Path(synthemol.__path__[0]).parent)")/rl/data/wuxi --folder
```

TODO: convert to zip file on google drive


## Reactions

The WuXi GalaXi consists of three phases of reactions. Phase 1 contains reactions with two building blocks, while phases 2 and 3 contain reactions with a core and two building blocks. Altogether, these reactions can be represented by 36 unique SMARTS, which are provided in `synthemol/reactions/wuxi.py`.


## Building blocks

Below, we describe the steps for processing the building blocks.


### SDF to SMILES

The `rl/data/wuxi/building_blocks_sdf` directory contains SDF files with WuXi GalaXi building blocks.

Convert the building blocks from SDF to SMILES.
```bash
for sdf_path in rl/data/wuxi/building_blocks_sdf/*.sdf; do
    sdf_file=$(basename "$sdf_path")
    csv_file=${sdf_file%.*}.csv
    echo "Converting ${csv_file}"
    chemfunc sdf_to_smiles \
        --data_path rl/data/wuxi/building_blocks_sdf/${sdf_file} \
        --save_path rl/data/wuxi/building_blocks_csv/${csv_file} \
        --properties wxid
done
```


### Merge building blocks files

Merge the WuXi building blocks files into a single file.
```bash
python scripts/data/merge_wuxi_building_blocks.py \
    --building_blocks_dir rl/data/wuxi/building_blocks_csv \
    --phase_2_cores_path rl/data/wuxi/phase2_cores.xlsx \
    --phase_3_cores_path rl/data/wuxi/phase3_cores.xlsx \
    --save_path rl/data/wuxi/building_blocks_raw.csv
```

This command prints the following dataset statistics.
```
Data size = 15,488
Data size after removing rows with missing IDs = 15,488
Data size after removing rows with missing SMILES = 15,484
Data size after removing rows with invalid SMILES = 15,484
Number of unique building block IDs = 15,484
Number of unique building block WuXi IDs = 15,399
Number of unique building block subsets = 25
Number of unique building block SMILES = 14,977
```


### Remove salts

Remove the salts from the building blocks. This will also canonicalize the SMILES using RDKit's canonicalization method.
```bash
chemfunc canonicalize_smiles \
    --data_path rl/data/wuxi/building_blocks_raw.csv \
    --save_path rl/data/wuxi/building_blocks.csv \
    --remove_salts \
    --delete_disconnected_mols
```

**Note:** This step is crucial to prevent errors in running reactions. Salts can cause reactions to create products that are the same as the reactants, leading to undesired infinite loops during molecule generation.


### Map WuXi reactions to building blocks

Determine which building blocks are valid in which WuXi reactions.
```bash
python scripts/data/map_wuxi_reactions_to_building_blocks.py \
    --building_blocks_path rl/data/wuxi/building_blocks.csv \
    --save_path rl/data/wuxi/reaction_to_building_blocks.pkl
```

## WuXi GalaXi molecules

Below, we describe the steps for processing the REAL Space molecules.

Note: The enumerated WuXi GalaXi requires over 1 TB of storage and a similar amount of RAM in order to process it. However, the following steps are a one-time operation, and only minimal storage and RAM are required to store the processed data and run SyntheMol. Each of the following commands run in <= 24 hours using 10 cores and <= 1 TB of RAM.


### Download REAL Space

Contact WuXi to request access to download the full enumerated WuXi GalaXi. The following commands assume that CSV file containing compounds from phases 1, 2, and 3 are in a single directory called `rl/data/wuxi/wuxi_galaxi`. In total, the WuXi GalaXi contains 16,146,071,436 molecules with 24,093,422 in phase 1; 6,864,759,544 in phase 2; and 9,257,218,470 in phase 3 from the 12/31/2022 download.

### Sample WuXi molecules

Randomly sample 10,000 WuXi GalaXi molecules. This is used for analysis of a representative sample of WuXi GalaXi molecules.
```bash
python scripts/data/sample_wuxi_galaxi.py \
    --data_dir rl/data/wuxi/wuxi_galaxi \
    --save_path rl/data/wuxi/random_wuxi_10k.csv \
    --num_molecules 10000
```

Randomly sample 7 million WuXi GalaXi molecules for a time-based comparison of Chemprop-RDKit versus SyntheMol-RL.
```bash
python scripts/data/sample_wuxi_galaxi.py \
    --data_dir rl/data/wuxi/wuxi_galaxi \
    --save_path rl/data/wuxi/random_wuxi_7m.csv \
    --num_molecules 7000000
```

Randomly sample 150 WuXi GalaXi molecules for experimental validation.
```bash
python scripts/data/sample_wuxi_galaxi.py \
    --data_dir rl/data/wuxi/wuxi_galaxi \
    --save_path rl/data/wuxi/random_wuxi_150.csv \
    --num_molecules 150
```
