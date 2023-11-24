"""Utility functions for generating molecules."""
import operator
import re
from functools import partial
from pathlib import Path
from typing import Callable

import pandas as pd

from synthemol.constants import ROLLOUT_COL, SCORE_COL, SMILES_COL
from synthemol.generate.node import Node


OPERATORS = {
    "<": operator.lt,
    "<=": operator.le,
    "==": operator.eq,
    "!=": operator.ne,
    ">": operator.gt,
    ">=": operator.ge,
}


def save_generated_molecules(
    nodes: list[Node],
    chemical_space_to_building_block_id_to_smiles: dict[str, dict[str, str]],
    model_names: tuple[str, ...],
    save_path: Path,
) -> None:
    """Save generated molecules to a CSV file.

    :param nodes: A list of Nodes containing molecules. Only nodes with a single molecule are saved.
    :param chemical_space_to_building_block_id_to_smiles: A dictionary mapping building block IDs to SMILES.
    :param model_names: A list of model names corresponding to the individual scores of the Nodes.
    :param save_path: A path to a CSV file where the molecules will be saved.
    """
    construction_dicts = []
    individual_score_dicts = []
    max_reaction_num = 0
    reaction_num_to_max_reactant_num = {}

    # Process construction logs and individual scores for each node
    for node in nodes:
        # Convert construction logs from lists to dictionaries
        construction_dict = {"num_reactions": len(node.construction_log)}
        max_reaction_num = max(max_reaction_num, len(node.construction_log))

        for reaction_index, reaction_log in enumerate(node.construction_log):
            reaction_num = reaction_index + 1
            construction_dict[
                f"reaction_{reaction_num}_chemical_space"
            ] = reaction_log.chemical_space
            construction_dict[f"reaction_{reaction_num}_id"] = reaction_log.reaction_id

            reaction_num_to_max_reactant_num[reaction_num] = max(
                reaction_num_to_max_reactant_num.get(reaction_num, 0),
                len(reaction_log.reactant_ids),
            )

            for reactant_index, reactant_id in enumerate(reaction_log.reactant_ids):
                reactant_num = reactant_index + 1
                construction_dict[
                    f"building_block_{reaction_num}_{reactant_num}_id"
                ] = reactant_id
                construction_dict[
                    f"building_block_{reaction_num}_{reactant_num}_smiles"
                ] = chemical_space_to_building_block_id_to_smiles[
                    reaction_log.chemical_space
                ].get(
                    reactant_id, ""
                )

        construction_dicts.append(construction_dict)

        # Get individual scores
        individual_score_dicts.append(
            {
                model_name: score
                for model_name, score in zip(model_names, node.individual_scores)
            }
        )

    # Specify column order for CSV file
    columns = [
        SMILES_COL,
        "node_id",
        "num_expansions",
        ROLLOUT_COL,
        SCORE_COL,
        *model_names,
        "Q_value",
        "num_reactions",
    ]

    for reaction_num in range(1, max_reaction_num + 1):
        columns.append(f"reaction_{reaction_num}_chemical_space")
        columns.append(f"reaction_{reaction_num}_id")

        for reactant_num in range(
            1, reaction_num_to_max_reactant_num[reaction_num] + 1
        ):
            columns.append(f"building_block_{reaction_num}_{reactant_num}_id")
            columns.append(f"building_block_{reaction_num}_{reactant_num}_smiles")

    # Save data
    save_path.parent.mkdir(parents=True, exist_ok=True)
    data = pd.DataFrame(
        data=[
            {
                SMILES_COL: ".".join(node.molecules),
                "node_id": node.node_id,
                "num_expansions": node.num_visits,
                ROLLOUT_COL: node.rollout_num,
                SCORE_COL: node.property_score,
                "Q_value": node.exploit_score,
                **construction_dict,
                **individual_score_dict,
            }
            for node, construction_dict, individual_score_dict in zip(
                nodes, construction_dicts, individual_score_dicts
            )
        ],
        columns=columns,
    )
    data.to_csv(save_path, index=False)


def compare_to_threshold(score: float, comparator: str, threshold: float,) -> bool:
    """Compares a score to a threshold using a comparator.

    :param score: A score.
    :param comparator: A comparator (e.g., >, >=, ==, !=, <, <=).
    :param threshold: A threshold.
    :return: Whether the score satisfies the comparator and threshold.
    """
    return OPERATORS[comparator](score, threshold)


def parse_success_threshold(success_threshold: str) -> Callable[[float], bool]:
    """Parses a success threshold string into a function that determines whether a molecule is successful.

    :param success_threshold: A string of the form "> 0.5" with a comparator and a threshold value.
    :return: A function that determines whether a molecule is successful according to the threshold.
    """
    # Create regex to match the success threshold pattern
    success_threshold_pattern = (
        r"^(?P<comparator>[<>!=]{1,2})\s*(?P<threshold>-*\d+(?:\.\d+)?)$"
    )

    # Parse success threshold
    match = re.match(success_threshold_pattern, success_threshold)

    if match is None:
        raise ValueError(f"Invalid success threshold: {success_threshold}")

    comparator = match.group("comparator")
    threshold = float(match.group("threshold"))

    # Create function that determines whether a molecule is successful according to the threshold
    compare_to_threshold_fn = partial(
        compare_to_threshold, comparator=comparator, threshold=threshold
    )

    return compare_to_threshold_fn
