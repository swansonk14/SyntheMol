"""Generate molecules combinatorially using a Monte Carlo tree search guided by a molecular property predictor."""
from tap import tapify

from SyntheMol.generate.generate import generate

tapify(generate)
