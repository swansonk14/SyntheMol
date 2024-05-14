import random
import subprocess
import textwrap
from distutils.dir_util import copy_tree
from pathlib import Path
from queue import PriorityQueue
from tempfile import TemporaryDirectory
from typing import Dict, List, Optional, Union
import os
import numpy as np
from meeko import MoleculePreparation, PDBQTWriterLegacy
from rdkit import Chem
from rdkit.Chem import AllChem
from wurlitzer import pipes
from pydantic import Field
from pydantic.dataclasses import dataclass


@dataclass
class Parameters:
    RECEPTOR_ROOT_PATH = "/workspace/Tyers/Reaction-GFN/data/molecule/targets"
    RECEPTOR_CENTERS = {
        "Mpro": [-20.458, 18.109, -26.914],
        "TBLR1": [-1.014, 42.097, 39.750],
        "ClpP": [-38.127, 45.671, -20.898],
        "LRRK2_WD40": [-16.386, -15.911, 7.779],
        "sEH": [-13.4, 26.3, -13.3],
    }
    RECEPTOR_BOX_SIZES = {
        "Mpro": [20, 20, 20],
        "TBLR1": [25, 25, 25],
        "ClpP": [17, 17, 17],
        "LRRK2_WD40": [25, 25, 25],
        "sEH": [20.013, 16.3, 18.5],
    }
    @property

    def RECEPTOR_PATHS(self):
        return {k: f"{self.RECEPTOR_ROOT_PATH}/{k}.pdbqt" for k in self.RECEPTOR_CENTERS.keys()}

    qv_dir: Union[Path, str] = "/workspace/Tyers/Vina-GPU-2.1/QuickVina2-GPU-2.1"
    receptor_path: Optional[Union[Path, str]] = None
    receptor_name: Optional[str] = "ClpP"
    log_dir: str = "./logs_docking"
    center: Optional[List[float]] = None
    size: Optional[List[float]] = None
    norm: float = 1.0
    failed_score: float = 0.0
    conformer_attempts: int = 20
    docking_attempts: int = 10
    beta: int = 4


class GPUDockingScoreProxy():
    def __init__(self, params: Parameters):
        if params.receptor_path is None and params.receptor_name is None:
            raise ValueError(
                "Expected either receptor_path or receptor_name to be specified."
            )
        if params.receptor_path is not None and params.receptor_name is not None:
            raise ValueError(
                "Expected only one of receptor_path and receptor_name to be specified."
            )

        if params.receptor_name is not None:
            assert params.center is None

            self.receptor_path = params.RECEPTOR_PATHS[params.receptor_name]
            self.center = params.RECEPTOR_CENTERS[params.receptor_name]
        else:
            self.receptor_path = receptor_path
            self.center = center

        self.qv_dir = params.qv_dir
        self.size = params.RECEPTOR_BOX_SIZES[params.receptor_name] if params.size is None else params.size
        self.norm = params.norm
        self.failed_score = params.failed_score
        self.conformer_attempts = params.conformer_attempts
        self.docking_attempts = params.docking_attempts
        self.log_dir = params.log_dir
        self.beta = params.beta

    def compute_scores(self, smiles: List[str]) -> List[float]:
        scores = self.dock_batch_qv2gpu(smiles)
        scores = [-x / self.norm for x in scores]
        return scores

    def _docking_attempt(self, smiles, n):
        """
        Uses customized QuickVina2-GPU (Tang et al.) implementation to
        calculate docking score against target of choice.

        1. Unique temp directory (for scores, config, and pdbqt inputs) is created.
        2. Molecules are protonated and relaxed, and converted to pdbqt.
        3. QV2GPU runs on the molecule pdbqts.
        4. Scores are read from the output file.
        5. Temp directory is removed.

        Note: Failure at any point in the pipeline (reading molecule, pdbqt conversion,
            score calculation) returns self.failed_score for that molecule.
        """
        ligand_directory = TemporaryDirectory(suffix="_ligand")
        config_directory = TemporaryDirectory(suffix="_config")
        config_path = str(Path(config_directory.name) / "config.txt")
        scores_path = str(Path(ligand_directory.name) / "scores.txt")

        qv2cfg = textwrap.dedent(
            f"""
            receptor = {self.receptor_path}
            ligand_directory = {ligand_directory.name}
            opencl_binary_path = {self.qv_dir}
            center_x = {self.center[0]}
            center_y = {self.center[1]}
            center_z = {self.center[2]}
            size_x = {self.size[0]}
            size_y = {self.size[1]}
            size_z = {self.size[2]}
            thread = 8000
            num_modes = 15
        """
        )

        with open(config_path, "w") as file:
            file.write(qv2cfg)

        initial_pdbqts = []
        docked_pdbqts = []
        indices = []
        count = 0

        for idx, smi in enumerate(smiles):
            attempt = 0
            pdbqt_string = None

            while pdbqt_string is None and (
                attempt == 0 or attempt < self.conformer_attempts
            ):
                attempt += 1

                try:
                    mol = Chem.MolFromSmiles(smi)
                    mol = Chem.AddHs(mol)
                    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
                    AllChem.UFFOptimizeMolecule(mol)
                    preparator = MoleculePreparation()
                    mol_setups = preparator.prepare(mol)
                    setup = mol_setups[0]
                    pdbqt_string, _, _ = PDBQTWriterLegacy.write_string(setup)
                except Exception as e:
                    print(f"Failed embedding attempt #{attempt} with error: '{e}'.")

            if pdbqt_string is None:
                continue

            initial_pdbqts.append(pdbqt_string)
            indices.append(idx)

            output_file_path = str(Path(ligand_directory.name) / f"{str(count)}.pdbqt")
            with open(output_file_path, "w") as file:
                file.write(pdbqt_string)

            count += 1

        command = [
            "./QuickVina2-GPU-2-1",
            "--config",
            config_path,
        ]
        result = subprocess.run(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            cwd=self.qv_dir,
        )

        if not Path(scores_path).exists():
            print(f"Failed score loading attempt #{n}: {result.stderr}")

            debug_path = Path(__file__).parents[3] / "debug" / ("%08x" % random.getrandbits(32))

            print(f"Saving debug files into {debug_path}...")

            with pipes():
                copy_tree(str(ligand_directory.name), str(debug_path / "ligand"))
                copy_tree(str(config_directory.name), str(debug_path / "config"))

            scores = None
        else:
            scores = [self.failed_score for _ in smiles]
            with open(scores_path, "r") as file:
                for i, line in enumerate(file):
                    try:
                        score = float(line.strip())

                    except ValueError:
                        print(f"Failed line reading attempt: '{line.strip()}'.")
                        score = self.failed_score

                    scores[indices[i]] = score

        ligand_directory.cleanup()
        config_directory.cleanup()

        return scores

    def reward_transform(self, rewards):
        return rewards


    def __call__(self, smiles: List[str]) -> np.array:
        for attempt in range(1, self.docking_attempts + 1):
            scores = self._docking_attempt(smiles, attempt)
            if scores is not None:
                with open(os.path.join(self.log_dir, "visited.txt"), 'a') as file:
                    # Write each molecule and its score to the file
                    for molecule, score in zip(smiles, scores):
                        file.write(f"{molecule}, {score}\n")
                
                raw_scores = -1 * np.array(scores, dtype=float)
                transformed_scores = self.reward_transform(raw_scores)

                print(f"Proxy Mean: {raw_scores.mean()}, Proxy Max: {raw_scores.max()}, Mean Reward: {transformed_scores.mean()}, Max Reward: {transformed_scores.max()}")
                return transformed_scores
        
        scores = [self.failed_score] * len(smiles)
        print(f"Proxy Mean: {np.array(scores).mean()}, Proxy Max: {np.array(scores).max()}")
        return np.array(scores, dtype=float)
