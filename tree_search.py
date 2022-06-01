"""Contains classes and functions for performing a tree search to generate molecules combinatorially."""
import math
from functools import partial
from random import Random
from typing import Any, Callable, List, Optional


# TODO: construction log
class TreeNode:
    """A node in a Monte Carlo Tree Search representing a step in the molecule construction process."""

    def __init__(self,
                 c_puct: float,
                 molecule: Optional[str] = None,
                 reagents: Optional[list[str]] = None,
                 W: float = 0.0,
                 N: int = 0,
                 P: float = 0.0) -> None:
        """Initializes the TreeNode object.

        :param c_puct: The hyperparameter that encourages exploration.
        :param molecule: The SMILES representing the currently constructed molecule.
        :param reagents: A list of SMILES containing the reagents that will be added to the current molecule.
        :param W: The sum of the node value.
        :param N: The number of times of arrival at this node.
        :param P: The property score (reward) of this node.
        """
        self.c_puct = c_puct
        self.molecule = molecule
        self.reagents = reagents
        self.W = W
        self.N = N
        self.P = P
        self.children: List[TreeNode] = []

    def Q(self) -> float:
        """Value that encourages exploitation of nodes with high reward."""
        return self.W / self.N if self.N > 0 else 0.0

    def U(self, n: int) -> float:
        """Value that encourages exploration of nodes with few visits."""
        return self.c_puct * self.P * math.sqrt(n) / (1 + self.N)

    def __str__(self) -> str:
        """Gets a string representation of the TreeNode.

        :return: The SMILES for the molecule followedd by the SMILES for each reagent, all space separated.
        """
        return ' '.join(
            ([self.molecule] if self.molecule is not None else []) +
            (self.reagents if self.reagents is not None else [])
        )

    def __hash__(self) -> int:
        """Computes the hash of the TreeNode, which is the hash of the string representation.

        :return: The hash of the TreeNode.
        """
        return hash(str(self))

    def __eq__(self, other: Any) -> bool:
        """Determines whether this TreeNode is equal to another object.

        :param other: Any object.
        :return: True if this TreeNode is equal to the other object (same molecule and reagents), False otherwise.
        """
        if not isinstance(other, TreeNode):
            return False

        return self.molecule == other.molecule and self.reagents == other.reagents


class TreeSearcher:
    """A class that runs a tree search to generate high scoring molecules."""

    def __init__(self,
                 scoring_fn: Callable[[str], float],
                 n_rollout: int,
                 c_puct: float,
                 num_expand_nodes: int) -> None:
        """Creates the TreeSearcher object.

        :param scoring_fn: A function that takes as input a SMILES representing a molecule and returns a score.
        :param n_rollout: The number of times to run the tree search.
        :param c_puct: The hyperparameter that encourages exploration.
        :param num_expand_nodes: The number of tree nodes to expand when extending the child nodes in the search tree.
        """
        self.scoring_fn = scoring_fn
        self.n_rollout = n_rollout
        self.c_puct = c_puct
        self.num_expand_nodes = num_expand_nodes
        self.random = Random(0)

        self.TreeNodeClass = partial(TreeNode, c_puct=c_puct)
        self.root = self.TreeNodeClass()
        self.state_map = {self.root: self.root}

    def rollout(self, tree_node: TreeNode) -> float:
        """Performs a Monte Carlo Tree Search rollout.

        :param tree_node: An TreeNode representing the root of the MCTS search.
        :return: The value (reward) of the rollout.
        """
        # Get a list of currently existing phrases (list of lists of indices)
        phrases = get_phrases(tree_node.mask)

        # Filter to only include phrases that are above the minimum phrase length
        shrinkable_phrases = [phrase for phrase in phrases if len(phrase) > self.min_phrase_length]

        # Stop the search and call this a leaf node if the node can't be expanded further
        if tree_node.percent_unmasked <= self.min_percent_unmasked or len(shrinkable_phrases) == 0:
            return tree_node.P

        # Expand if this node has never been visited
        if len(tree_node.children) == 0:
            # Maintain a set of all the masks added as children of the tree
            tree_children_masks = set()

            # Include left and right ends of each phrase as options to mask
            indices_to_mask = [i for phrase in shrinkable_phrases for i in [phrase[0], phrase[-1]]]

            # If not yet at max number of phrases, allow breaking a phrase
            # (as long as the resulting phrases are long enough)
            if len(phrases) < self.max_phrases:
                indices_to_mask += [i for phrase in shrinkable_phrases for i in
                                    phrase[self.min_phrase_length:-self.min_phrase_length]]

            # Limit the number of MCTS nodes to expand
            self.random.shuffle(indices_to_mask)
            indices_to_mask = indices_to_mask[:self.num_expand_nodes]

            # For each index, prune (mask) it and create a new TreeNode (if it doesn't already exist)
            for index_to_mask in indices_to_mask:
                # Create new mask with additional index masked
                new_mask = list(tree_node.mask)
                new_mask[index_to_mask] = 0
                new_mask = tuple(new_mask)

                # Check the state map and merge with an existing mask if available
                new_node = self.state_map.setdefault(new_mask, self.TreeNodeClass(mask=new_mask))

                # Add the mask to the children of the tree
                if new_mask not in tree_children_masks:
                    tree_node.children.append(new_node)
                    tree_children_masks.add(new_mask)

            # For each child in the tree, compute its reward
            for child in tree_node.children:
                if child.P == 0:
                    child.P = text_score(words=child.words, mask=child.mask, scoring_fn=self.scoring_fn)

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
                 scoring_fn: Callable[[str], float],
                 n_rollout: int,
                 c_puct: float,
                 num_expand_nodes: int) -> None:
        """Creates the TreeSearchRunner object.

        :param scoring_fn: A function that takes as input a SMILES representing a molecule and returns a score.
        :param n_rollout: The number of times to run the tree search.
        :param c_puct: The hyperparameter that encourages exploration.
        :param num_expand_nodes: The number of tree nodes to expand when extending the child nodes in the search tree.
        """
        self.scoring_fn = scoring_fn
        self.n_rollout = n_rollout
        self.c_puct = c_puct
        self.num_expand_nodes = num_expand_nodes

    def search(self) -> List[TreeNode]:
        """Runs the tree search.

        :return: A list of TreeNode objects representing text masks sorted from highest to lowest reward.
        """
        return TreeSearcher(
            scoring_fn=self.scoring_fn,
            n_rollout=self.n_rollout,
            c_puct=self.c_puct,
            num_expand_nodes=self.num_expand_nodes
        ).search()
