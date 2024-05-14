from synthemol.models.bengio2021flow import load_original_model, mol2graph
from typing import List, Optional
import torch
import torch_geometric.data as gd
import numpy as np
from rdkit.Chem import AllChem
from rdkit import Chem
from pydantic import Field
from pydantic.dataclasses import dataclass
import os


@dataclass
class Parameters:
    model = load_original_model()
    device: str = "cuda"
    log_dir: str = "./logs_seh"
    beta: int = 8


class SEHProxy():
    def __init__(self, params: Parameters):
        self.device = params.device
        self.model = params.model
        self.log_dir = params.log_dir
        self.beta = params.beta
        self.model.to(self.device)

    def reward_transform(self, rewards):
        return rewards

    def __call__(self, smiles: List[str]) -> np.array:
        mols = [Chem.MolFromSmiles(s) for s in smiles]
        graphs = [mol2graph(i) for i in mols]
        smiles_to_save = [Chem.MolToSmiles(mol) for mol in mols]
        is_valid = torch.tensor([i is not None for i in graphs]).bool()
        batch = gd.Batch.from_data_list([i for i in graphs if i is not None])
        batch.to(self.device)
        self.model.to(self.device)
        raw_scores = self.model(batch).reshape((-1,)).data.cpu()
        raw_scores[raw_scores.isnan()] = 0
        scores_to_save = list(raw_scores)
        #with open(os.path.join(self.log_dir, "visited.txt"), 'a') as file:
        #    # Write each molecule and its score to the file
        #    for molecule, score in zip(smiles_to_save, scores_to_save):
        #        file.write(f"{molecule}, {score}\n")

        transformed_scores = self.reward_transform(raw_scores).clip(1e-4, 100).reshape((-1,))
        #print(f"Proxy Mean: {raw_scores.mean()}, Proxy Max: {raw_scores.max()}, Mean Reward: {transformed_scores.mean()}, Max Reward: {transformed_scores.max()}")
        return np.array(transformed_scores, dtype=float)