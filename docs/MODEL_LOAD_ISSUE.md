Issue: `dimer_potential_v1.joblib` intentionally stores a bundle dict containing the estimator plus scaler and metadata.

Status: RESOLVED (2026-02-27)

Details:
- The bundle keys include `model`, `scaler`, `feature_names`, `metrics`, `model_name`, `model_version`, and `saved_at`.
- Loading should extract both the estimator and the scaler from the bundle rather than assuming a bare estimator object.

Correct load pattern (apply everywhere):

```python
bundle = joblib.load("models/dimer_potential_v1.joblib")
model = bundle["model"]
scaler = bundle["scaler"]
```

Files updated:
- `phase3_step3_predict.py` — now loads the bundle and applies `scaler.transform` before `model.predict`.

Rationale:
- Preserving the bundle ensures `scaler` and `feature_names` travel with the estimator and avoids creating duplicate artifacts.

Logged by: automation (session run) — 2026-02-27
