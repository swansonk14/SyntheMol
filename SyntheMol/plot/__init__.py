"""SyntheMol.plot module."""
from SyntheMol.plot.plot_auc import compute_curve, plot_auc
from SyntheMol.plot.plot_building_block_scores import plot_building_block_scores
from SyntheMol.plot.plot_building_block_vs_molecule_scores import plot_building_block_vs_molecule_scores
from SyntheMol.plot.plot_generated_molecule_analysis import (
    plot_building_block_usage,
    plot_generated_molecule_analysis,
    plot_similarity,
    plot_reactions_numbers,
    plot_reaction_usage,
    plot_scores
)
from SyntheMol.plot.plot_mcts_over_time import plot_mcts_over_time
from SyntheMol.plot.plot_model_generalization import plot_model_generalization
from SyntheMol.plot.plot_molecule_analysis import plot_molecule_analysis
from SyntheMol.plot.plot_real_counts import plot_real_counts
from SyntheMol.plot.plot_regression_values import plot_regression_values
from SyntheMol.plot.plot_toxicity import plot_toxicity
