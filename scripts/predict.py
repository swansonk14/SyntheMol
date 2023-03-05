"""Make predictions with a model or ensemble of models and save them to a file."""
from tap import tapify

from SyntheMol.models import predict_model


if __name__ == '__main__':
    tapify(predict_model)
