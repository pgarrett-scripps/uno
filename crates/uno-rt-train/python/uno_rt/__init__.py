"""Optional Uno CNN training. Search and evaluation do not need this package."""
__version__ = "0.1.0"


def train(data, out, *, units, method, **kwargs):
    """Train/fine-tune a CNN; return its report. See ``uno rt train --help``."""
    from .training import train as implementation
    return implementation(data, out, units=units, method=method, **kwargs)


def load_model(path):
    """Return (PyTorch CNN, normalization/architecture info, metadata)."""
    from .artifact import load_model as implementation
    return implementation(path)
