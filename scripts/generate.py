"""Generate molecules combinatorially using a Monte Carlo tree search guided by a molecular property predictor."""
from tap import tapify

from SyntheMol.generate import generate


if __name__ == '__main__':
    tapify(generate)
