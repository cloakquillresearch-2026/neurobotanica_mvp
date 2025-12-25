import json
import numpy as np
import sys
from pathlib import Path

# Add project root to path
PROJECT_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from src.analysis.nextgen_dimer_ingest import load_nextgen_dimer_evidence
from src.ml_models.nextgen_dimer_input import encode_dimer_features
from src.ml_models.nextgen_dimer_fusion import prepare_batch
from src.ml_models.nextgen_dimer_heads import NextGenDimerModel

# Minimal SMILES lookup for demo (expand as needed)
SMILES_LOOKUP = {
    'THC': 'CCCCCc1cc(c2c(c1)OC(C3C2C=C(CC3)C)(C)C)O',
    'CBD': 'CC(C)(C1=CC(=C(C=C1)O)C2=CC(=C(C=C2)O)C)C',
    'CBG': 'CCCCCc1cc(c2c(c1)OC(C3C2C=C(CC3)C)(C)C)O',
    'CBDV': 'CC(C)(C1=CC(=C(C=C1)O)C2=CC(=C(C=C2)O)C)C',
    'myrcene': 'CC(=C)CCC=C(C)C',
    'linalool': 'CC(C)=CCC=C(C)CO',
}

if __name__ == "__main__":
    # Load data
    dimers = load_nextgen_dimer_evidence("data/processed/nextgen_dimers.json")
    entries = [d.__dict__ for d in dimers]
    X = prepare_batch(entries)
    # Dummy targets for demonstration
    y_synergy = np.array([1 if e['interaction_effect'] == 'synergistic' else 0 for e in entries])
    y_effect = np.array([float(e.get('synergy_score', 1.0)) for e in entries])
    # Train/test split (all for train in this demo)
    model = NextGenDimerModel()
    model.fit(X, y_synergy, y_effect)
    synergy_pred, effect_pred = model.predict(X)
    print("Synergy predictions:", synergy_pred)
    print("Effect size predictions:", effect_pred)
