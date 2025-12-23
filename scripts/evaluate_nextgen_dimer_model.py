import numpy as np
from sklearn.metrics import accuracy_score, roc_auc_score, mean_squared_error, r2_score
from src.analysis.nextgen_dimer_ingest import load_nextgen_dimer_evidence
from src.ml_models.nextgen_dimer_input import encode_dimer_features
from src.ml_models.nextgen_dimer_fusion import prepare_batch
from src.ml_models.nextgen_dimer_heads import NextGenDimerModel

if __name__ == "__main__":
    # Load data
    dimers = load_nextgen_dimer_evidence("data/processed/nextgen_dimers.json")
    entries = [d.__dict__ for d in dimers]
    X = prepare_batch(entries)
    y_synergy = np.array([1 if e['interaction_effect'] == 'synergistic' else 0 for e in entries])
    y_effect = np.array([float(e.get('synergy_score', 1.0)) for e in entries])

    # Train/test split (all for train in this demo)
    model = NextGenDimerModel()
    model.fit(X, y_synergy, y_effect)
    synergy_pred, effect_pred = model.predict(X)

    # Evaluation metrics
    synergy_pred_label = (synergy_pred > 0.5).astype(int)
    acc = accuracy_score(y_synergy, synergy_pred_label)
    auc = roc_auc_score(y_synergy, synergy_pred)
    mse = mean_squared_error(y_effect, effect_pred)
    r2 = r2_score(y_effect, effect_pred)

    print(f"Synergy classification: Accuracy={acc:.2f}, AUC={auc:.2f}")
    print(f"Effect size regression: MSE={mse:.4f}, R2={r2:.2f}")
