"""Trains a machine learning property prediction model."""
from tap import tapify

from SyntheMol.models import train_model


if __name__ == '__main__':
    tapify(train_model)
