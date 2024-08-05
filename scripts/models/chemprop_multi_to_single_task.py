"""Converts a Chemprop multi-task model to a series of single-task models."""
from copy import deepcopy
from pathlib import Path

import torch.nn as nn
from chemprop.utils import load_args, load_checkpoint, load_scalers, save_checkpoint
from tqdm import tqdm


def chemprop_multi_to_single_task(model_path: Path, save_dir: Path,) -> None:
    """Converts a Chemprop multi-task model to a series of single-task models.

    :param model_path: The path to the Chemprop multi-task model or directory of multi-task models.
    :param save_dir: The directory to save the single-task models to.
    """
    # Get model paths
    if model_path.is_dir():
        model_paths = sorted(path for path in model_path.glob("**/*.pt"))
    else:
        model_paths = [model_path]

    # Convert each model
    for model_index, path in enumerate(tqdm(model_paths)):
        # Load components
        args = load_args(path)
        model = load_checkpoint(path)
        (
            scaler,
            features_scaler,
            atom_descriptor_scaler,
            bond_descriptor_scaler,
            atom_bond_scaler,
        ) = load_scalers(path)

        # Save single-task models
        for task_index, task_name in enumerate(args.task_names):
            # Create task save directory
            task_save_dir = save_dir / task_name
            task_save_dir.mkdir(parents=True, exist_ok=True)

            # Copy components
            new_args = deepcopy(args)
            new_model = deepcopy(model)
            new_scaler = deepcopy(scaler)

            # Set task name
            new_args.task_names = [task_name]

            # Extract single-task model weights
            new_model.readout[-1] = nn.Linear(model.readout[-1].in_features, 1)

            new_model.readout[-1].weight.data = (
                model.readout[-1].weight[task_index : task_index + 1].data
            )
            new_model.readout[-1].bias.data = (
                model.readout[-1].bias[task_index : task_index + 1].data
            )

            # Extract single-task scaler stats
            if new_scaler is not None:
                new_scaler.means = scaler.means[task_index : task_index + 1]
                new_scaler.stds = scaler.stds[task_index : task_index + 1]

            # Save single-task model
            save_checkpoint(
                path=task_save_dir / f"model_{model_index}.pt",
                model=new_model,
                scaler=new_scaler,
                features_scaler=features_scaler,
                atom_descriptor_scaler=atom_descriptor_scaler,
                bond_descriptor_scaler=bond_descriptor_scaler,
                atom_bond_scaler=atom_bond_scaler,
                args=new_args,
            )


if __name__ == "__main__":
    from tap import tapify

    tapify(chemprop_multi_to_single_task)
