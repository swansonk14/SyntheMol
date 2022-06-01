"""Contains classes and functions for performing a tree search to generate molecules combinatorially."""
import itertools
import math
from functools import partial
from typing import Any, Callable, Optional

import numpy as np
from rdkit import Chem

from real_reactions import convert_to_mol, Reaction, REAL_REACTIONS


class TreeNode:
    """A node in a tree search representing a step in the molecule construction process."""

    def __init__(self,
                 c_puct: float,
                 scoring_fn: Callable[[str], float],
                 fragments: Optional[tuple[str]] = None,
                 construction_log: Optional[tuple[dict[str, Any]]] = None) -> None:
        """Initializes the TreeNode object.

        :param c_puct: The hyperparameter that encourages exploration.
        :param scoring_fn: A function that takes as input a SMILES representing a molecule and returns a score.
        :param fragments: A tuple of SMILES containing the fragments for the next reaction.
                         The first element is the currently constructed molecule while the remaining elements
                         are the fragments that are about to be added.
        :param construction_log: A tuple of dictionaries containing information about each reaction.
        """
        self.c_puct = c_puct
        self.scoring_fn = scoring_fn
        self.fragments = fragments if fragments is not None else tuple()
        self.construction_log = construction_log if construction_log is not None else tuple()
        # TODO: maybe change sum (really mean) to max since we don't care about finding the best leaf node, just the best node along the way?
        self.W = 0.0  # The sum of the leaf node values for leaf nodes that descend from this node.
        self.N = 0  # The number of leaf nodes that have been visited from this node.
        self._P = None  # The property score of this node. (Computed lazily.)
        self.children: set[TreeNode] = set()

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

    @property
    def P(self) -> float:
        """The property score of this node."""
        if self._P is None:
            self._P = self.compute_score(fragments=self.fragments, scoring_fn=self.scoring_fn)

        return self._P

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

    def __hash__(self) -> int:
        return hash(self.fragments)

    def __eq__(self, other: Any) -> bool:
        if not isinstance(other, TreeNode):
            return False

        return self.fragments == other.fragments


class TreeSearcher:
    """A class that runs a tree search to generate high scoring molecules."""

    def __init__(self,
                 fragment_to_index: dict[str, int],
                 reagent_to_fragments: dict[str, list[str]],
                 max_reactions: int,
                 scoring_fn: Callable[[str], float],
                 n_rollout: int,
                 c_puct: float,
                 num_expand_nodes: int) -> None:
        """Creates the TreeSearcher object.

        :param fragment_to_index: A dictionary mapping fragment SMILES to their index in the fragment file.
        :param reagent_to_fragments: A dictionary mapping a reagent SMARTS to all fragment SMILES that match that reagent.
        :param max_reactions: The maximum number of reactions to use to construct a molecule.
        :param scoring_fn: A function that takes as input a SMILES representing a molecule and returns a score.
        :param n_rollout: The number of times to run the tree search.
        :param c_puct: The hyperparameter that encourages exploration.
        :param num_expand_nodes: The number of tree nodes to expand when extending the child nodes in the search tree.
        """
        self.fragment_to_index = fragment_to_index
        self.reagent_to_fragments = reagent_to_fragments
        self.all_fragments = sorted(
            fragment for fragments in self.reagent_to_fragments.values() for fragment in fragments
        )
        self.max_reactions = max_reactions
        self.scoring_fn = scoring_fn
        self.n_rollout = n_rollout
        self.c_puct = c_puct
        self.num_expand_nodes = num_expand_nodes
        self.rng = np.random.default_rng(seed=0)

        self.TreeNodeClass = partial(TreeNode, c_puct=c_puct, scoring_fn=scoring_fn)
        self.root = self.TreeNodeClass()
        self.state_map: dict[TreeNode, TreeNode] = {self.root: self.root}

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
            for matched_reagent_indices in itertools.product(*reagent_matches_per_mol):
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
    def get_all_next_fragment_nodes(self, node: TreeNode) -> list[TreeNode]:
        if len(node.fragments) == 0:
            next_fragments = self.all_fragments
        else:
            next_fragments = self.get_next_fragments(fragments=node.fragments)

        new_nodes = [
            self.TreeNodeClass(
                fragments=(*node.fragments, next_fragment),
                construction_log=node.construction_log
            )
            for next_fragment in next_fragments
        ]

        return new_nodes

    def get_reactions_for_fragments(self, fragments: tuple[str]) -> list[tuple[Reaction, dict[str, int]]]:
        mols = [convert_to_mol(fragment, add_hs=True) for fragment in fragments]
        matching_reactions = []

        for reaction in REAL_REACTIONS:
            if len(mols) != reaction.num_reagents:
                continue

            # For each mol, get a list of indices of reagents it matches
            reagent_matches_per_mol = self.get_reagent_matches_per_mol(reaction=reaction, mols=mols)

            # Include every assignment of fragments to reagents that fills all the reagents
            for matched_reagent_indices in itertools.product(*reagent_matches_per_mol):
                if len(set(matched_reagent_indices)) == reaction.num_reagents:
                    fragment_to_reagent_index = dict(zip(fragments, matched_reagent_indices))
                    matching_reactions.append((reaction, fragment_to_reagent_index))

        return matching_reactions

    def run_all_reactions(self, node: TreeNode) -> list[TreeNode]:
        matching_reactions = self.get_reactions_for_fragments(fragments=node.fragments)

        product_nodes = []
        for reaction, fragment_to_reagent_index in matching_reactions:
            # Put fragments in the right order for the reaction
            fragments = sorted(node.fragments, key=lambda frag: fragment_to_reagent_index[frag])

            # Run reaction
            products = reaction.run_reactants(fragments)

            if len(products) == 0:
                raise ValueError('Reaction failed to produce products.')

            assert all(len(product) == 1 for product in products)

            # Convert product mols to SMILES, remove Hs, and deduplicate (preserving order for random reproducibility)
            products = list(dict.fromkeys(Chem.MolToSmiles(Chem.RemoveHs(product[0])) for product in products))

            # Create reaction log
            reaction_log = {
                'reaction_id': reaction.id,
                'reagent_ids': [self.fragment_to_index.get(fragment, -1) for fragment in fragments],
                'reagent_smiles': [fragment for fragment in fragments]
            }

            product_nodes += [
                self.TreeNodeClass(
                    fragments=(product,),
                    construction_log=node.construction_log + (reaction_log,)
                )
                for product in products
            ]

        return product_nodes

    def expand_node(self, node: TreeNode) -> list[TreeNode]:
        # Run all valid reactions on the current fragments to generate new molecules
        new_nodes = self.run_all_reactions(node=node)

        # Add all possible next fragments to the current fragment
        new_nodes += self.get_all_next_fragment_nodes(node=node)

        return new_nodes

    def rollout(self, node: TreeNode) -> float:
        """Performs a Monte Carlo Tree Search rollout.

        :param node: An TreeNode representing the root of the MCTS search.
        :return: The value (reward) of the rollout.
        """
        # Stop the search if we've reached the maximum number of reactions
        if node.num_reactions >= self.max_reactions:
            return node.P

        # Expand if this node has never been visited
        if len(node.children) == 0:
            # Expand the node both by running reactions with the current fragments and adding new fragments
            new_nodes = self.expand_node(node=node)

            if len(new_nodes) == 0:
                raise NotImplementedError  # TODO: handle this case of no valid next fragments

            # Limit the number of nodes to expand
            if len(new_nodes) > self.num_expand_nodes:
                new_nodes = self.rng.choice(new_nodes, size=self.num_expand_nodes, replace=False)

            # Add the expanded nodes as children
            for new_node in new_nodes:

                # Check whether this node was already added as a child
                if new_node not in node.children:

                    # Check the state map and merge with an existing node if available
                    new_node = self.state_map.setdefault(new_node, new_node)

                    # Add the new node as a child of the tree node
                    node.children.add(new_node)

        # Select the best child node and unroll it
        sum_count = sum(child.N for child in node.children)
        selected_node = max(node.children, key=lambda x: x.Q() + x.U(n=sum_count))
        v = self.rollout(node=selected_node)
        selected_node.W += v
        selected_node.N += 1

        return v

    def search(self) -> list[TreeNode]:
        """Runs the tree search.

        :return: A list of TreeNode objects representing text masks sorted from highest to lowest reward.
        """
        for _ in range(self.n_rollout):
            self.rollout(node=self.root)

        # Get all the nodes representing fully constructed molecules
        nodes = [node for _, node in self.state_map.items() if node.num_fragments == 1]

        # Sort by highest reward and break ties by preferring less unmasked results
        nodes = sorted(nodes, key=lambda node: node.P, reverse=True)

        return nodes


class TreeSearchRunner:
    """A class that creates instances of TreeSearcher to search for high scoring molecules."""

    def __init__(self,
                 fragment_to_index: dict[str, int],
                 reagent_to_fragments: dict[str, list[str]],
                 max_reactions: int,
                 scoring_fn: Callable[[str], float],
                 n_rollout: int,
                 c_puct: float,
                 num_expand_nodes: int) -> None:
        """Creates the TreeSearchRunner object.

        :param fragment_to_index: A dictionary mapping fragment SMILES to their index in the fragment file.
        :param reagent_to_fragments: A dictionary mapping a reagent SMARTS to all fragment SMILES that match that reagent.
        :param max_reactions: The maximum number of reactions to use to construct a molecule.
        :param scoring_fn: A function that takes as input a SMILES representing a molecule and returns a score.
        :param n_rollout: The number of times to run the tree search.
        :param c_puct: The hyperparameter that encourages exploration.
        :param num_expand_nodes: The number of tree nodes to expand when extending the child nodes in the search tree.
        """
        self.fragment_to_index = fragment_to_index
        self.reagent_to_fragments = reagent_to_fragments
        self.max_reactions = max_reactions
        self.scoring_fn = scoring_fn
        self.n_rollout = n_rollout
        self.c_puct = c_puct
        self.num_expand_nodes = num_expand_nodes

    def search(self) -> list[TreeNode]:
        """Runs the tree search.

        :return: A list of TreeNode objects representing text masks sorted from highest to lowest reward.
        """
        return TreeSearcher(
            fragment_to_index=self.fragment_to_index,
            reagent_to_fragments=self.reagent_to_fragments,
            max_reactions=self.max_reactions,
            scoring_fn=self.scoring_fn,
            n_rollout=self.n_rollout,
            c_puct=self.c_puct,
            num_expand_nodes=self.num_expand_nodes
        ).search()
