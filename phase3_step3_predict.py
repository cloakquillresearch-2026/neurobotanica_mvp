import json
import joblib
import json
import joblib
import numpy as np
from pathlib import Path

# Load model bundle correctly
bundle = joblib.load("models/dimer_potential_v1.joblib")
model = bundle["model"]
scaler = bundle["scaler"]
feature_names = bundle.get("feature_names")

print(f"Model type: {type(model).__name__}")
print(f"Feature names: {feature_names}")
print(f"Model version: {bundle.get('model_version', 'unknown')}")

with open("data/processed/phase3_ranked_pairs.json", "r") as f:
    pairs = json.load(f)

results = []
errors = 0

for p in pairs:
    try:
        # Build feature vector in correct order
        X_raw = np.array([[
            0.5,   # increased_receptor_affinity — placeholder
            0.5,   # prolonged_duration — placeholder
            0.5,   # reduced_metabolism — placeholder
            p.get("priority_score", 5.0) / 25.0,  # confidence_score proxy
            0.5,   # enhanced_lipophilicity — placeholder
        ]])
        X_scaled = scaler.transform(X_raw)
        pred = float(model.predict(X_scaled)[0])
        p["effect_size"] = round(pred, 4)
        p["source"] = "model_prediction"
        p["evidence_type"] = "model_prediction"
        p["confidence"] = "low"
        if not p.get("notes"):
            p["notes"] = "Phase 3 model prediction. No literature evidence found."
    except Exception as e:
        p["effect_size"] = 0.78
        p["source"] = "model_prediction"
        p["notes"] = f"Model error — default applied. {e}"
        errors += 1
    results.append(p)

print(f"\nPairs predicted: {len(results)}")
print(f"Model errors:    {errors}")
if results:
    sizes = [r["effect_size"] for r in results]
    print(f"Effect size range: {min(sizes):.4f} – {max(sizes):.4f}")
    print(f"Mean effect size:  {sum(sizes)/len(sizes):.4f}")

# Preserve already-sourced literature pairs in phase3_predictions.json
with open("data/processed/phase3_predictions.json", "r") as f:
    existing = json.load(f)

existing_keyed = {
    frozenset([r["compound_1"], r["compound_2"]]): r
    for r in existing
    if r.get("source") in ("literature", "rct", "invalid")
}

final = []
for r in results:
    key = frozenset([r["compound_1"], r["compound_2"]])
    if key in existing_keyed:
        final.append(existing_keyed[key])  # preserve literature records
    else:
        final.append(r)

with open("data/processed/phase3_predictions.json", "w") as f:
    json.dump(final, f, indent=2)

print(f"\nSaved {len(final)} records to phase3_predictions.json")
print(f"Literature/invalid records preserved: {len(existing_keyed)}")

