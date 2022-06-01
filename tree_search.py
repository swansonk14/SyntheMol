"""Contains classes and functions for performing a tree search to generate molecules combinatorially."""
import math
from functools import partial
from itertools import product
from typing import Any, Callable, List, Optional

import numpy as np
from rdkit import Chem

from real_reactions import convert_to_mol, Reaction, REAL_REACTIONS


# TODO: construction logs
class TreeNode:
    """A node in a tree search representing a step in the molecule construction process."""

    def __init__(self,
                 c_puct: float,
                 scoring_fn: Callable[[str], float],
                 fragments: Optional[tuple[str]] = None,
                 construction_log: Optional[dict[str, Any]] = None) -> None:
        """Initializes the TreeNode object.

        :param c_puct: The hyperparameter that encourages exploration.
        :param scoring_fn: A function that takes as input a SMILES representing a molecule and returns a score.
        :param fragments: A tuple of SMILES containing the fragments for the next reaction.
                         The first element is the currently constructed molecule while the remaining elements
                         are the fragments that are about to be added.
        :param construction_log: A list of dictionaries containing information about each reaction.
        """
        self.c_puct = c_puct
        self.scoring_fn = scoring_fn
        self.fragments = fragments or tuple()
        self.construction_log = construction_log or []
        # TODO: maybe change sum (really mean) to max since we don't care about finding the best leaf node, just the best node along the way?
        self.W = 0.0  # The sum of the leaf node values for leaf nodes that descend from this node.
        self.N = 0  # The number of leaf nodes that have been visited from this node.
        self.P = self.compute_score(fragments=fragments, scoring_fn=scoring_fn)  # The property score of this node.
        self.children: List[TreeNode] = []

    @classmethod
    def compute_score(cls, fragments: tuple[str], scoring_fn: Callable[[str], float]) -> float:
        """Computes the score of the fragments.

        :param fragments: A tuple of SMILES containing the fragments for the next reaction.
                         The first element is the currently constructed molecule while the remaining elements
                         are the fragments that are about to be added.
        :param scoring_fn: A function that takes as input a SMILES representing a molecule and returns a score.
        :return: The score of the fragments.
        """
        # TODO: change this!!! to weight the fragments differently
        return sum(scoring_fn(fragment) for fragment in fragments)

    def Q(self) -> float:
        """Value that encourages exploitation of nodes with high reward."""
        return self.W / self.N if self.N > 0 else 0.0

    def U(self, n: int) -> float:
        """Value that encourages exploration of nodes with few visits."""
        return self.c_puct * self.P * math.sqrt(n) / (1 + self.N)

    @property
    def num_fragments(self) -> int:
        """Returns the number of fragments for the upcoming reaction.

        :return: The number of fragments for the upcoming reaction.
        """
        return len(self.fragments)

    @property
    def num_reactions(self) -> int:
        """Returns the number of reactions used so far to generate the current molecule.

        :return: The number of reactions used so far to generate the current molecule.
        """
        return len(self.construction_log)


class TreeSearcher:
    """A class that runs a tree search to generate high scoring molecules."""

    def __init__(self,
                 reagent_to_fragments: dict[str, list[str]],
                 max_reactions: int,
                 scoring_fn: Callable[[str], float],
                 n_rollout: int,
                 c_puct: float,
                 num_expand_nodes: int) -> None:
        """Creates the TreeSearcher object.

        :param reagent_to_fragments: A dictionary mapping a reagent SMARTS to all fragment SMILES that match that reagent.
        :param max_reactions: The maximum number of reactions to use to construct a molecule.
        :param scoring_fn: A function that takes as input a SMILES representing a molecule and returns a score.
        :param n_rollout: The number of times to run the tree search.
        :param c_puct: The hyperparameter that encourages exploration.
        :param num_expand_nodes: The number of tree nodes to expand when extending the child nodes in the search tree.
        """
        self.reagent_to_fragments = reagent_to_fragments
        self.all_fragment_tuples = sorted(
            (fragment,) for fragments in self.reagent_to_fragments.values() for fragment in fragments
        )
        self.max_reactions = max_reactions
        self.scoring_fn = scoring_fn
        self.n_rollout = n_rollout
        self.c_puct = c_puct
        self.num_expand_nodes = num_expand_nodes
        self.rng = np.random.default_rng(seed=0)

        self.TreeNodeClass = partial(TreeNode, c_puct=c_puct, scoring_fn=scoring_fn)
        self.root = self.TreeNodeClass()
        self.state_map: dict[tuple[str, ...], TreeNode] = {tuple(): self.root}

    # TODO: documentation
    @classmethod
    def get_reagent_matches_per_mol(cls, reaction: Reaction, mols: list[Chem.Mol]) -> list[list[int]]:
        return [
            [
                reagent_index
                for reagent_index, reagent in enumerate(reaction.reagents)
                if reagent.has_substruct_match(mol)
            ]
            for mol in mols
        ]

    # TODO: documentation
    def get_next_fragments(self, fragments: tuple[str]) -> list[str]:
        mols = [convert_to_mol(fragment, add_hs=True) for fragment in fragments]

        # For each reaction these fragments can participate in, get all the unfilled reagents
        unfilled_reagents = set()
        for reaction in REAL_REACTIONS:
            reagent_indices = set(range(reaction.num_reagents))

            # Skip reaction if there's no room to add more reagents
            if len(mols) >= reaction.num_reagents:
                continue

            # For each mol, get a list of indices of reagents it matches
            reagent_matches_per_mol = self.get_reagent_matches_per_mol(reaction=reaction, mols=mols)

            # Loop through products of reagent indices that the mols match to
            # and for each product, if it matches to all separate reagents,
            # then include the missing reagents in the set of unfilled reagents
            for matched_reagent_indices in product(*reagent_matches_per_mol):
                matched_reagent_indices = set(matched_reagent_indices)

                if len(matched_reagent_indices) == len(mols):
                    unfilled_reagents |= {reaction.fragments[index] for index in
                                          reagent_indices - matched_reagent_indices}

        if len(unfilled_reagents) == 0:
            return []

        # Get all the fragments that match the other reagents
        available_fragments = list(dict.fromkeys(
            fragment
            for reagent in sorted(unfilled_reagents, key=lambda reagent: reagent.smarts)
            for fragment in self.reagent_to_fragments[reagent.smarts]
        ))

        return available_fragments

    # TODO: documentation
    def get_next_fragment_tuples(self, fragments: tuple[str]) -> list[tuple[str, ...]]:
        if len(fragments) == 0:
            return self.all_fragment_tuples

        next_fragments = self.get_next_fragments(fragments=fragments)
        next_fragment_tuples = [(*fragments, next_fragment) for next_fragment in next_fragments]

        return next_fragment_tuples

    def get_reactions_for_fragments(self, fragments: tuple[str]) -> list[tuple[Reaction, dict[str, int]]]:
        mols = [convert_to_mol(fragment, add_hs=True) for fragment in fragments]
        matching_reactions = []

        for reaction in REAL_REACTIONS:
            if len(mols) != reaction.num_reagents:
                continue

            # For each mol, get a list of indices of reagents it matches
            reagent_matches_per_mol = self.get_reagent_matches_per_mol(reaction=reaction, mols=mols)

            # Include every assignment of fragments to reagents that fills all the reagents
            for matched_reagent_indices in product(*reagent_matches_per_mol):
                if len(set(matched_reagent_indices)) == reaction.num_reagents:
                    fragment_to_reagent_index = dict(zip(fragments, matched_reagent_indices))
                    matching_reactions.append((reaction, fragment_to_reagent_index))

        return matching_reactions

    def run_all_reactions_for_fragments(self, fragments: tuple[str]) -> list[str]:
        matching_reactions = self.get_reactions_for_fragments(fragments=fragments)

        products = []
        for reaction, fragment_to_reagent_index in matching_reactions:
            # Put fragments in the right order for the reaction
            sorted_fragments = sorted(fragments, key=lambda frag: fragment_to_reagent_index[frag])

            # Run reaction
            reaction_products = reaction.run_reactants(sorted_fragments)

            if len(reaction_products) == 0:
                raise ValueError('Reaction failed to produce products.')

            assert all(len(product) == 1 for product in reaction_products)

            # Convert product mols to SMILES, remove Hs, and deduplicate (preserving order for random reproducibility)
            products += list(dict.fromkeys(Chem.MolToSmiles(Chem.RemoveHs(product[0])) for product in reaction_products))

            # TODO: reaction log
            # # Create reaction log
            # reaction_log = {
            #     'reaction_id': reaction.id,
            #     'reagent_ids': [fragment_to_index.get(fragment, -1) for fragment in fragments],
            #     'reagent_smiles': [fragment for fragment in fragments]
            # }

        return products

    def rollout(self, tree_node: TreeNode) -> float:
        """Performs a Monte Carlo Tree Search rollout.

        :param tree_node: An TreeNode representing the root of the MCTS search.
        :return: The value (reward) of the rollout.
        """
        # Stop the search if we've reached the maximum number of reactions
        if tree_node.num_reactions >= self.max_reactions:
            return tree_node.P

        # Expand if this node has never been visited
        if len(tree_node.children) == 0:
            # Run all valid reactions on the current fragments to generate new molecules
            new_fragment_tuples = [
                (product,) for product in self.run_all_reactions_for_fragments(fragments=tree_node.fragments)
            ]

            # Add all possible next fragments to the current fragments
            new_fragment_tuples += self.get_next_fragment_tuples(fragments=tree_node.fragments)

            if len(new_fragment_tuples) == 0:
                raise NotImplementedError  # TODO: handle this case of no valid next fragments

            # Limit the number of nodes to expand
            if len(new_fragment_tuples) > self.num_expand_nodes:
                new_fragment_tuples = self.rng.choice(new_fragment_tuples, size=self.num_expand_nodes, replace=False)

            # Add the expanded nodes as children
            child_fragment_tuples = set()
            for new_fragment_tuple in new_fragment_tuples:

                # Check whether this node was already added as a child
                if new_fragment_tuple not in child_fragment_tuples:

                    # Check the state map and merge with an existing node if available
                    if new_fragment_tuple not in self.state_map:
                        self.state_map[new_fragment_tuple] = self.TreeNodeClass(fragments=new_fragment_tuple)

                    new_node = self.state_map[new_fragment_tuple]

                    # Add the new node as a child of the tree node
                    tree_node.children.append(new_node)
                    child_fragment_tuples.add(new_fragment_tuple)

        # Select the best child node and unroll it
        sum_count = sum(child.N for child in tree_node.children)
        selected_node = max(tree_node.children, key=lambda x: x.Q() + x.U(n=sum_count))
        v = self.rollout(tree_node=selected_node)
        selected_node.W += v
        selected_node.N += 1

        return v

    def search(self) -> List[TreeNode]:
        """Runs the tree search.

        :return: A list of TreeNode objects representing text masks sorted from highest to lowest reward.
        """
        for _ in range(self.n_rollout):
            self.rollout(tree_node=self.root)

        # Get all the nodes representing fully constructed molecules
        nodes = [node for _, node in self.state_map.items() if node.num_fragments == 1]

        # Sort by highest reward and break ties by preferring less unmasked results
        nodes = sorted(nodes, key=lambda node: node.P, reverse=True)

        return nodes


class TreeSearchRunner:
    """A class that creates instances of TreeSearcher to search for high scoring molecules."""

    def __init__(self,
                 reagent_to_fragments: dict[str, list[str]],
                 max_reactions: int,
                 scoring_fn: Callable[[str], float],
                 n_rollout: int,
                 c_puct: float,
                 num_expand_nodes: int) -> None:
        """Creates the TreeSearchRunner object.

        :param reagent_to_fragments: A dictionary mapping a reagent SMARTS to all fragment SMILES that match that reagent.
        :param max_reactions: The maximum number of reactions to use to construct a molecule.
        :param scoring_fn: A function that takes as input a SMILES representing a molecule and returns a score.
        :param n_rollout: The number of times to run the tree search.
        :param c_puct: The hyperparameter that encourages exploration.
        :param num_expand_nodes: The number of tree nodes to expand when extending the child nodes in the search tree.
        """
        self.reagent_to_fragments = reagent_to_fragments
        self.max_reactions = max_reactions
        self.scoring_fn = scoring_fn
        self.n_rollout = n_rollout
        self.c_puct = c_puct
        self.num_expand_nodes = num_expand_nodes

    def search(self) -> List[TreeNode]:
        """Runs the tree search.

        :return: A list of TreeNode objects representing text masks sorted from highest to lowest reward.
        """
        return TreeSearcher(
            reagent_to_fragments=self.reagent_to_fragments,
            max_reactions=self.max_reactions,
            scoring_fn=self.scoring_fn,
            n_rollout=self.n_rollout,
            c_puct=self.c_puct,
            num_expand_nodes=self.num_expand_nodes
        ).search()
