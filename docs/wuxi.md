# WuXi GalaXi Data Processing

Instructions for processing WuXi GalaXi REAL reactions and building blocks. Assumes relevant data has already been downloaded (see [docs/README.md](README.md)).

## Download

Download the WuXi GalaXi data with the following command.
```bash
gdown "https://drive.google.com/drive/folders/1BCM5P-tDjLFLkyiEJX9p4-6fw0w5-qRg?usp=drive_link" -O $(python -c "import synthemol; from pathlib import Path; print(Path(synthemol.__path__[0]).parent)")/data/wuxi --folder
```

TODO: convert to zip file on google drive

## Reactions

The WuXi GalaXi consists of three phases of reactions. Phase 1 contains reactions with two building blocks, while phases 2 and 3 contain reactions with a core and two building blocks. Altogether, these reactions can be represented by 36 unique SMARTS, which are provided in `synthemol/reactions/wuxi.py`.


## Building blocks

Below, we describe the steps for processing the building blocks.

### SDF to SMILES

The `data/wuxi/building_blocks_sdf` directory contains SDF files with WuXi GalaXi building blocks.

Convert the building blocks from SDF to SMILES.
```bash
for sdf_path in data/wuxi/building_blocks_sdf/*.sdf; do
    sdf_file=$(basename "$sdf_path")
    csv_file=${sdf_file%.*}.csv
    echo "Converting ${csv_file}"
    chemfunc sdf_to_smiles \
        --data_path data/wuxi/building_blocks_sdf/${sdf_file} \
        --save_path data/wuxi/building_blocks_csv/${csv_file} \
        --properties wxid
done
```

### Merge building blocks files

Merge the WuXi building blocks files into a single file.
```bash
python scripts/data/merge_wuxi_building_blocks.py \
    --data_dir data/wuxi/building_blocks_csv \
    --save_path data/wuxi/building_blocks.csv
```

This command prints the following dataset statistics.
```
Data size = 14,468
Data size after removing rows with missing SMILES = 14,464
Number of unique building block IDs = 14,460
Number of unique building block SMILES = 14,455
Number of unique building block subsets = 15
```

Note: Unlike the Enamine REAL building blocks, the WuXi building blocks do not have salts, so there is no need to remove them.


### Map WuXi reactions to building blocks

TODO

Determine which building blocks are valid in which WuXi reactions.
```bash
python scripts/data/map_real_reactions_to_building_blocks.py \
    --data_dir data/Data/4_real_space/full_real \
    --save_path data/Data/4_real_space/reaction_to_building_blocks.pkl
```

Filter out building blocks that do not match the reaction templates.
```bash
python scripts/data/filter_real_reactions_to_building_blocks.py \
    --reaction_to_building_blocks_path data/Data/4_real_space/reaction_to_building_blocks.pkl \
    --save_path data/Data/4_real_space/reaction_to_building_blocks_filtered.pkl
```
