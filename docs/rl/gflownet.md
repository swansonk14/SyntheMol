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

Run the following experiments to generate molecules.

GFlowNet optimizing for S. aureus activity and solubility.
```bash
python src/gflownet/tasks/seh_frag_moo.py \
    --objectives s_aureus solubility \
    --log_dir logs/s_aureus_solubility
```

GFlowNet optimizing for S. aureus activity, solubility, and SA score.
```bash
python src/gflownet/tasks/seh_frag_moo.py \
    --objectives s_aureus solubility sa \
    --log_dir logs/s_aureus_solubility_sa
```

Next, extract the results from the sqlite database to CSV.
```bash
for MODEL in s_aureus_solubility s_aureus_solubility_sa
do
python scripts/extract_results.py \
    --results_path logs/${MODEL}/final/generated_mols_0.db \
    --save_path logs/${MODEL}/final/molecules.csv
done
```

Rename the columns and rescale solubility and SA score.
```bash
for MODEL in s_aureus_solubility s_aureus_solubility_sa
do
python -c "import pandas as pd;
path = 'logs/${MODEL}/final/molecules.csv';
data = pd.read_csv(path);
data = data.rename(columns={'smi': 'smiles', 'fr_0': 'S. aureus', 'fr_1': 'Solubility', 'fr_2': 'sa_score'});
data['Solubility'] = 14 * data['Solubility'] - 10;
if 'sa_score' in data: data['sa_score'] = -9 * data['sa_score'] + 10;
data.to_csv(path, index=False)"
done
```

## Evaluate generated molecules

Due to Python version incompatibilities with `chemfunc`, use the `synthemol` environment to run the following commands. (Or create a new conda environment with Python version 3.10 and run `pip install chemfunc==1.0.5`.)

Compute novelty of the generated molecules.
```bash
for MODEL in s_aureus_solubility s_aureus_solubility_sa
do
chemfunc nearest_neighbor \
    --data_path logs/${MODEL}/final/molecules.csv \
    --reference_data_path ../SyntheMol/rl/data/s_aureus/s_aureus_hits.csv \
    --reference_name train_hits \
    --metric tversky

chemfunc nearest_neighbor \
    --data_path logs/${MODEL}/final/molecules.csv \
    --reference_data_path ../SyntheMol/rl/data/chembl/chembl.csv \
    --reference_name chembl \
    --metric tversky
done
```

Select hit molecules that satisfy novelty, diversity, and efficacy thresholds (optionally including synthesiability).
```bash
python ../SyntheMol/scripts/data/select_molecules.py \
    --data_path logs/s_aureus_solubility/final/molecules.csv \
    --save_molecules_path logs/s_aureus_solubility/final/hits.csv \
    --save_analysis_path logs/s_aureus_solubility/final/analysis.csv \
    --score_columns "S. aureus" "Solubility" \
    --score_comparators ">=0.5" ">=-4" \
    --novelty_threshold 0.6 0.6 \
    --similarity_threshold 0.6 \
    --select_num 150 \
    --sort_column "S. aureus" \
    --descending

python ../SyntheMol/scripts/data/select_molecules.py \
    --data_path logs/s_aureus_solubility_sa/final/molecules.csv \
    --save_molecules_path logs/s_aureus_solubility_sa/final/hits.csv \
    --save_analysis_path logs/s_aureus_solubility_sa/final/analysis.csv \
    --score_columns "S. aureus" "Solubility" "sa_score" \
    --score_comparators ">=0.5" ">=-4" "<=4" \
    --novelty_threshold 0.6 0.6 \
    --similarity_threshold 0.6 \
    --select_num 150 \
    --sort_column "S. aureus" \
    --descending
```

Visualize hits.
```bash
for MODEL in s_aureus_solubility s_aureus_solubility_sa
do
chemfunc visualize_molecules \
    --data_path logs/${MODEL}/final/hits.csv \
    --save_dir logs/${MODEL}/final/hits
done
```
