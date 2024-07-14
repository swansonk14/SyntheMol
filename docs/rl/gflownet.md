# Generating molecules with GFlowNet

Below are instructions for generating molecules with GFlowNet, a generative flow-based model.


## Installation

Clone the forked repository.

```bash
git clone https://github.com/swansonk14/gflownet
cd gflownet
git checkout antibiotics
```

Install dependencies.

```bash
pip install -e . --find-links https://data.pyg.org/whl/torch-1.13.1+cu117.html
```

If package dependencies seem not to work, you may need to install the exact frozen versions listed `requirements/`, i.e. `pip install -r requirements/main_3.9.txt`.

For antibiotics applications, install additional dependencies:
```bash
pip install chemprop==1.6.1
pip install descriptastorus==2.6.1
pip install typed-argument-parser==1.9.0
```

**Note:** If you get the issue `ImportError: libXrender.so.1: cannot open shared object file: No such file or directory`, run `conda install -c conda-forge xorg-libxrender`.


## Generate molecules 

Run GFlowNet optimizing for S. aureus activity, solubility, and SA score.
```bash
python src/gflownet/tasks/seh_frag_moo.py \
    --objectives s_aureus solubility sa \
    --log_dir logs/s_aureus_solubility_sa
```

Extract the results from the sqlite database to CSV.
```bash
python scripts/extract_results.py \
    --results_path logs/s_aureus_solubility_sa/final/generated_mols_0.db \
    --save_path logs/s_aureus_solubility_sa/final/molecules.csv
```

Rename the columns and rescale solubility and SA score.
```bash
python -c "import pandas as pd;
path = 'logs/s_aureus_solubility_sa/final/molecules.csv';
data = pd.read_csv(path);
data = data.rename(columns={'smi': 'smiles', 'fr_0': 'S. aureus', 'fr_1': 'Solubility', 'fr_2': 'sa_score'});
data['Solubility'] = 14 * data['Solubility'] - 10;
if 'sa_score' in data: data['sa_score'] = -9 * data['sa_score'] + 10;
data.to_csv(path, index=False)"
```

## Evaluate generated molecules

Switch to the SyntheMol conda environment for the following commands.

Compute novelty of the generated molecules.
```bash
chemfunc nearest_neighbor \
    --data_path logs/$s_aureus_solubility_sa/final/molecules.csv \
    --reference_data_path ../SyntheMol/rl/data/s_aureus/s_aureus_hits.csv \
    --reference_name train_hits \
    --metric tversky

chemfunc nearest_neighbor \
    --data_path logs/s_aureus_solubility_sa/final/molecules.csv \
    --reference_data_path ../SyntheMol/rl/data/chembl/chembl.csv \
    --reference_name chembl \
    --metric tversky
```

Select hit molecules that satisfy novelty, diversity, and efficacy thresholds (including synthesizability).
```bash
python ../SyntheMol/scripts/data/select_molecules.py \
    --data_path logs/s_aureus_solubility_sa/final/molecules.csv \
    --save_molecules_path logs/s_aureus_solubility_sa/final/hits.csv \
    --save_analysis_path logs/s_aureus_solubility_sa/final/analysis.csv \
    --score_columns "S. aureus" "Solubility" "sa_score" \
    --score_comparators ">=0.5" ">=-4" "<=4" \
    --novelty_threshold 0.6 0.6 \
    --similarity_threshold 0.6 \
    --select_num 14 \
    --sort_column "S. aureus" \
    --descending
```

Visualize hits.
```bash
chemfunc visualize_molecules \
    --data_path logs/s_aureus_solubility_sa/final/hits.csv \
    --save_dir logs/s_aureus_solubility_sa/final/hits
```
