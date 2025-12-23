import numpy as np
from typing import List, Dict, Any
from src.ml_models.nextgen_dimer_input import encode_dimer_features

# Multi-modal fusion: simple concatenation (can be replaced with attention/fusion layers)
def fuse_features(feature_list: List[np.ndarray]) -> np.ndarray:
    return np.concatenate(feature_list, axis=-1)

# Prepare batch input for model
def prepare_batch(entries: List[Dict[str, Any]]) -> np.ndarray:
    features = [encode_dimer_features(e) for e in entries]
    return np.stack(features)

# Example usage:
# batch = prepare_batch(list_of_dimer_entries)
