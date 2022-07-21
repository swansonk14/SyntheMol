import time

from pathlib import Path
from typing import Union, Optional
from importlib.util import find_spec
from rdkit import Chem
from rdkit.Chem import rdFMCS
import pandas as pd 

from tap import Tap

class Args(Tap):
    data_path: Path  # Path to CSV file containing generated molecules.
    save_path: Path  # Path to directory where plots and results will be saved.
    # train_hits_path: Optional[Path] = None  # Path to CSV file containing hits from the training set for computing novelty. #Might implement at a later time depending on computational limits
    smiles_column: str = 'smiles'  # The name of the column containing SMILES in data_path.
    min_numAtoms: Optional[int] = 0


#TODO: Potential scoring function depending on numAtoms in MCS (currently just binary value over threshold)

#Created using list from following link (critically important) https://academic.oup.com/cid/article/63/8/1087/2389125 + Pubchem to find SMILES
start_time = time.time()
reference_antibiotics = {'aminoglycosides': ("gentamicin", "CC(C1CCC(C(O1)OC2C(CC(C(C2O)OC3C(C(C(CO3)(C)O)NC)O)N)N)N)NC"),
'ansamycins': ('rifampin', 'CC1C=CC=C(C(=O)NC2=C(C(=C3C(=C2O)C(=C(C4=C3C(=O)C(O4)(OC=CC(C(C(C(C(C(C1O)C)O)C)OC(=O)C)C)OC)C)C)O)O)C=NN5CCN(CC5)C)C'),
'carbapenems': ('meropenem', 'CC1C2C(C(=O)N2C(=C1SC3CC(NC3)C(=O)N(C)C)C(=O)O)C(C)O'),
'cephalosporins': ('ceftrioxone', 'CON=C(c1csc(n1)N)C(=O)NC1C(=O)N2C1SCC(=C2C(=O)O)CSc1nc(=O)c(=O)[nH]n1C' ), #from https://www.guidetopharmacology.org/GRAC/LigandDisplayForward?ligandId=5326
'phosphonic_acid_derivatives': ('fosfomycin', 'CC1C(O1)P(=O)(O)O'),
'glycopeptides' : ('vancomycin', 'CC1C(C(CC(O1)OC2C(C(C(OC2OC3=C4C=C5C=C3OC6=C(C=C(C=C6)C(C(C(=O)NC(C(=O)NC5C(=O)NC7C8=CC(=C(C=C8)O)C9=C(C=C(C=C9O)O)C(NC(=O)C(C(C1=CC(=C(O4)C=C1)Cl)O)NC7=O)C(=O)O)CC(=O)N)NC(=O)C(CC(C)C)NC)O)Cl)CO)O)O)(C)N)O'),
'glycylcyclines': ('tigecycline', 'CC(C)(C)NCC(=O)NC1=CC(=C2CC3CC4C(C(=O)C(=C(C4(C(=O)C3=C(C2=C1O)O)O)O)C(=O)N)N(C)C)N(C)C'),
'lipopeptides' : ('daptomycin', 'CCCCCCCCCC(=O)NC(CC1=CNC2=CC=CC=C21)C(=O)NC(CC(=O)N)C(=O)NC(CC(=O)O)C(=O)NC3C(OC(=O)C(NC(=O)C(NC(=O)C(NC(=O)CNC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)CNC3=O)CCCN)CC(=O)O)C)CC(=O)O)CO)C(C)CC(=O)O)CC(=O)C4=CC=CC=C4N)C' ),
'macrolides' : ('erythromycin', 'CCC1C(C(C(C(=O)C(CC(C(C(C(C(C(=O)O1)C)OC2CC(C(C(O2)C)O)(C)OC)C)OC3C(C(CC(O3)C)N(C)C)O)(C)O)C)C)O)(C)O'),
'monobactams': ('aztreonam', 'CC1C(C(=O)N1S(=O)(=O)O)NC(=O)C(=NOC(C)(C)C(=O)O)C2=CSC(=N2)N'),
'oxazolidinones': ('linezolid', 'CC(=O)NCC1CN(C(=O)O1)C2=CC(=C(C=C2)N3CCOCC3)F'),
'penicillins' : ('ampicillin', 'CC1(C(N2C(S1)C(C2=O)NC(=O)C(C3=CC=CC=C3)N)C(=O)O)C'),
'quinolones': ('ciprofloxacin', 'C1CC1N2C=C(C(=O)C3=CC(=C(C=C32)N4CCNCC4)F)C(=O)O'),
'misc' : ('isoniazid', 'C1=CN=CC=C1C(=O)NN'),
 }

reference_antibiotics_classes = list(reference_antibiotics.keys())

#Uses reference antibiotics and target molecule and returns numAtoms, antibiotic class it shares and the SMARTs string of similarity (last only if specified)
def find_most_similar_class(reference_antibiotics = reference_antibiotics, mol = Union[str, Chem.Mol], min_numAtoms = int) -> tuple[int, str, str]: 
    if type(mol) == str: #Creates Mol object if input is SMILES str
        mol = Chem.MolFromSmiles(mol)
    mcs_paired_list = [] #List of MCS generated objects. Compared to reference antibiotics
    for antibiotic_class in reference_antibiotics:
        mcs_paired_list.append(rdFMCS.FindMCS([mol, Chem.MolFromSmiles(reference_antibiotics[antibiotic_class][1])],
        ringMatchesRingOnly = True, #Atoms must both be in (or not be in) rings to count towards MCS
        completeRingsOnly= True, #Must fully include ring in MCS
        ))
    mcs_num_atoms = [element.numAtoms for element in mcs_paired_list]
    max_mcs_num_atoms = int(max(mcs_num_atoms))
    return([max_mcs_num_atoms if max_mcs_num_atoms >= min_numAtoms else 0 , reference_antibiotics_classes[mcs_num_atoms.index(max_mcs_num_atoms)], mcs_paired_list[mcs_num_atoms.index(max_mcs_num_atoms)].smartsString])


def save_file(save_path = Path, data = pd.DataFrame):
    save_path.parent.mkdir(parents=True, exist_ok=True)
    data.to_csv(save_path, index=False)

def find_mcs_score(args: Args) -> None:
    data = pd.read_csv(args.data_path)
    # data.drop_duplicates(subset=reference_smiles_column, inplace=True)
    # data.sort_values(by=reference_smiles_column, ignore_index=True, inplace=True)
    print('Calculating MCS scores')
    mcs_scores = [find_most_similar_class(mol = element, min_numAtoms= args.min_numAtoms) for element in data[args.smiles_column]]
    mcs_scores = pd.DataFrame(mcs_scores, columns = ['numAtoms', 'antibiotic_class', 'mcs_smarts']) 
    data = pd.concat([data, mcs_scores], axis = 1)
    print(data)

    print('Saving')
    save_file(save_path=args.save_path, data=data)
    

    


if __name__ == '__main__':
    find_mcs_score(Args().parse_args())
print( "%s seconds" %(time.time()-start_time))
#Create reference path of handpicked antibiotic compounds from each class

    #kcluster reference path
    #retrieve one from each cluster
#MCS to each cluster (What is the efficiency of this? How hard can we push this clustering?)