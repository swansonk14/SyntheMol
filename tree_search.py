"""Contains classes and functions for performing a tree search to generate molecules combinatorially."""
import math
from functools import partial
from itertools import product
from random import Random
from typing import Any, Callable, List, Optional

from rdkit import Chem

from real_reactions import convert_to_mol, Reaction, REAL_REACTIONS


# TODO: all documentation
def get_reagent_matches_per_mol(reaction: Reaction, mols: list[Chem.Mol]) -> list[list[int]]:
    return [
        [
            reagent_index
            for reagent_index, reagent in enumerate(reaction.reagents)
            if reagent.has_substruct_match(mol)
        ]
        for mol in mols
    ]


def get_next_fragments(fragments: list[str], reagent_to_fragments: dict[str, list[str]]) -> list[str]:
    mols = [convert_to_mol(fragment, add_hs=True) for fragment in fragments]

    # For each reaction these fragments can participate in, get all the unfilled reagents
    unfilled_reagents = set()
    for reaction in REAL_REACTIONS:
        reagent_indices = set(range(reaction.num_reagents))

        # Skip reaction if there's no room to add more reagents
        if len(mols) >= reaction.num_reagents:
            continue

        # For each mol, get a list of indices of reagents it matches
        reagent_matches_per_mol = get_reagent_matches_per_mol(reaction=reaction, mols=mols)

        # Loop through products of reagent indices that the mols match to
        # and for each product, if it matches to all separate reagents,
        # then include the missing reagents in the set of unfilled reagents
        for matched_reagent_indices in product(*reagent_matches_per_mol):
            matched_reagent_indices = set(matched_reagent_indices)

            if len(matched_reagent_indices) == len(mols):
                unfilled_reagents |= {reaction.reagents[index] for index in reagent_indices - matched_reagent_indices}

    if len(unfilled_reagents) == 0:
        return []

    # Get all the fragments that match the other reagents
    available_fragments = list(dict.fromkeys(
        fragment
        for reagent in sorted(unfilled_reagents, key=lambda reagent: reagent.smarts)
        for fragment in reagent_to_fragments[reagent.smarts]
    ))

    return available_fragments


# TODO: construction logs
class TreeNode:
    """A node in a tree search representing a step in the molecule construction process."""

    def __init__(self,
                 c_puct: float,
                 scoring_fn: Callable[[str], float],
                 reagents: Optional[tuple[str]] = None,
                 construction_log: Optional[dict[str, Any]] = None) -> None:
        """Initializes the TreeNode object.

        :param c_puct: The hyperparameter that encourages exploration.
        :param scoring_fn: A function that takes as input a SMILES representing a molecule and returns a score.
        :param reagents: A list of SMILES containing the reagents for the next reaction.
                         The first element is the currently constructed molecule while the remaining elements
                         are the reagents that are about to be added.
        :param construction_log: A list of dictionaries containing information about each reaction.
        """
        self.c_puct = c_puct
        self.scoring_fn = scoring_fn
        self.reagents = reagents or tuple()
        self.construction_log = construction_log or []
        self._W = 0.0  # The sum of the node value.
        self._N = 0  # The number of times of arrival at this node.
        self._P = None  # The property score (reward) of this node.
        self.children: List[TreeNode] = []

    def compute_score(self) -> float:
        """Computes the score of this node in terms of the scores of the reagents.

        :return: The score of this node.
        """
        # TODO: change this!!! to weight the reagents differently
        return sum(self.scoring_fn(reagent) for reagent in self.reagents)

    # TODO: maybe change sum (really mean) to max since we don't care about finding the best leaf node, just the best node along the way?
    @property
    def W(self) -> float:
        """The sum of the leaf node values for leaf nodes that descend from this node."""
        return self._W

    @W.setter
    def W(self, W: float) -> None:
        self._W = W

    @property
    def N(self) -> int:
        """The number of leaf nodes that have been visited from this node."""
        return self._N

    @N.setter
    def N(self, N: int) -> None:
        self._N = N

    @property
    def P(self) -> float:
        """The score of this node."""
        if self._P is None:
            self._P = self.compute_score()

        return self._P

    def Q(self) -> float:
        """Value that encourages exploitation of nodes with high reward."""
        return self.W / self.N if self.N > 0 else 0.0

    def U(self, n: int) -> float:
        """Value that encourages exploration of nodes with few visits."""
        return self.c_puct * self.P * math.sqrt(n) / (1 + self.N)

    @property
    def num_reagents(self) -> int:
        """Returns the number of reagents for the upcoming reaction.

        :return: The number of reagents for the upcoming reaction.
        """
        return len(self.reagents)

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
        self.all_fragments = sorted(
            fragment for fragments in self.reagent_to_fragments.values() for fragment in fragments
        )
        self.max_reactions = max_reactions
        self.scoring_fn = scoring_fn
        self.n_rollout = n_rollout
        self.c_puct = c_puct
        self.num_expand_nodes = num_expand_nodes
        self.random = Random(0)

        self.TreeNodeClass = partial(TreeNode, c_puct=c_puct, scoring_fn=scoring_fn)
        self.root = self.TreeNodeClass()
        self.state_map = {tuple(): self.root}

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
            # TODO: fill this in
            # Select the first fragment
            if tree_node.num_reagents == 0:
                pass
            # Select the second fragment
            elif tree_node.num_reagents == 1:
                pass
            # Either complete the reaction with two reagents or select and complete the reaction with three reagents
            elif tree_node.num_reagents == 2:
                pass
            else:
                raise ValueError(f'Cannot handle reactions with "{tree_node.num_reagents}" reagents.')

            # Limit the number of nodes to expand
            self.random.shuffle(new_nodes)
            new_nodes = new_nodes[:self.num_expand_nodes]

            # Add the expanded nodes as children
            child_states = set()
            for new_node in new_nodes:
                # Check whether this node was already added as a child
                if new_node.state not in child_states:
                    # Check the state map and merge with an existing node if available
                    new_node = self.state_map.setdefault(TreeNode.compute_state(new_node.reagents), new_node)

                    # Add the new node as a child of the tree node
                    tree_node.children.append(new_node)
                    child_states.add(new_node.state)

            # For each child in the tree, compute its reward
            for child_node in tree_node.children:
                if child_node.P == 0:
                    child_node.P = state_score(reagents=child_node.reagents, scoring_fn=self.scoring_fn)

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

        # Get all the nodes representing fully constructed molecules (i.e., no reagents about to be added)
        nodes = [node for _, node in self.state_map.items() if node.reagents is not None]

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
            reagent_to_fragments=reagent_to_fragments,
            max_reactions=self.max_reactions,
            scoring_fn=self.scoring_fn,
            n_rollout=self.n_rollout,
            c_puct=self.c_puct,
            num_expand_nodes=self.num_expand_nodes
        ).search()
