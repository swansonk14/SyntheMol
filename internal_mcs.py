#Add in closest molecule
import time
import multiprocessing
from pathlib import Path
from typing import Union, Optional
from importlib.util import find_spec
from xmlrpc.client import Boolean
from rdkit import Chem
from rdkit.Chem import rdFMCS
import pandas as pd 
from functools import partial
from contextlib import closing

import numpy as np
from tqdm import tqdm
from tap import Tap

from multiprocessing import Pool
#Also deduplicates list

class Args(Tap):
    data_path: Path  # Path to CSV file containing generated molecules.
    save_path: Path  # Path to directory where plots and results will be saved.
    # train_hits_path: Optional[Path] = None  # Path to CSV file containing hits from the training set for computing novelty. #Might implement at a later time depending on computational limits
    smiles_column: str = 'smiles'  # The name of the column containing SMILES in data_path.
    deduplicate: Boolean = False


#TODO: Potential scoring function depending on numAtoms in MCS (currently just binary value over threshold)



#Uses reference antibiotics and target molecule and returns numAtoms, antibiotic class it shares and the SMARTs string of similarity (last only if specified)
def find_largest_MCS(target_mol =  Chem.Mol, list_mols = list) -> tuple[int, str]: 
    for mol in list_mols:
        mcs_track = 0
        mcs = rdFMCS.FindMCS(
            [target_mol, mol],
            ringMatchesRingOnly = True,
            completeRingsOnly = True)
        mcs_size = mcs.numAtoms
        if mcs_size>mcs_track:
            mcs_track = mcs_size
            mcs_smile = mcs.smartsString
    
    return(mcs_track, mcs_smile)


    # mcs_paired_list = [] #List of MCS generated objects. Compared to reference antibiotics
    # for antibiotic_class in reference_antibiotics:
    #     mcs_paired_list.append(rdFMCS.FindMCS([mol, Chem.MolFromSmiles(reference_antibiotics[antibiotic_class][1])],
    #     ringMatchesRingOnly = True, #Atoms must both be in (or not be in) rings to count towards MCS
    #     completeRingsOnly= True, #Must fully include ring in MCS
    #     ))
    # mcs_num_atoms = [element.numAtoms for element in mcs_paired_list]
    # max_mcs_num_atoms = int(max(mcs_num_atoms))
    # return([max_mcs_num_atoms if max_mcs_num_atoms >= min_numAtoms else 0 , reference_antibiotics_classes[mcs_num_atoms.index(max_mcs_num_atoms)], mcs_paired_list[mcs_num_atoms.index(max_mcs_num_atoms)].smartsString])



def save_file(save_path = Path, data = pd.DataFrame):
    save_path.parent.mkdir(parents=True, exist_ok=True)
    data.to_csv(save_path, index=False, mode= 'w+')

def internal_mcs(args: Args) -> None:
    data = pd.read_csv(args.data_path)
    if args.deduplicate:
        data.drop_duplicates(subset=args.smiles_column, inplace=True)
        data.sort_values(by=args.smiles_column, ignore_index=True, inplace=True)
    print('Converting smiles to Mols')
    list_mols = []
    for smile in tqdm(data[args.smiles_column]):
        list_mols.append(Chem.MolFromSmiles(smile))
    
    time0 = time.time()
    print('Calculating MCS scores')
    mcs_generator = partial(find_largest_MCS, list_mols = list_mols)
    print(multiprocessing.cpu_count())
    with closing(Pool(4)) as pool:
        mcs_scores = np.array(list(tqdm(pool.map(mcs_generator, list_mols), total = len(list_mols), desc = 'Finding largest MCS')))

    print(mcs_scores)
    mcs_scores = pd.DataFrame(mcs_scores)
    # mcs_scores = [find_largest_MCS(target_mol = element, list_mols = list_mols) for element in tqdm(list_mols)]
    # mcs_scores = pd.DataFrame(mcs_scores, columns = ['numAtoms','mcs_smarts']) 
    data = pd.concat([data, mcs_scores], axis = 1)
    print(data)

    print('Saving')
    save_file(save_path=args.save_path, data=data)
    print(time.time() - time0)
    

    


if __name__ == '__main__':
    internal_mcs(Args().parse_args())


#Create reference path of handpicked antibiotic compounds from each class

    #kcluster reference path
    #retrieve one from each cluster
#MCS to each cluster (What is the efficiency of this? How hard can we push this clustering?)