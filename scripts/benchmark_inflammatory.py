"""Simple benchmark for TS-PS-001 predictive call latency and cache effectiveness."""
import os
import time
import json
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from backend.services.inflammatory_synergy_engine import get_inflammatory_synergy_engine

# Opt-in bypass for dev benchmarking
os.environ['TS_PS_001_DEV_BYPASS'] = 'true'

engine = get_inflammatory_synergy_engine()

payload = {
    "biomarkers": {"tnf_alpha": 12.5, "il6": 6.0, "crp": 3.1, "il1b": 0.8},
    "condition_profile": {"conditions": [{"name": "inflammation", "severity": 5}]},
    "available_kingdoms": ["cannabis", "fungal", "plant"]
}

iterations = int(os.getenv('TS_PS_001_BENCH_ITERS', '1000'))

# Warm-up
_ = engine.predict_inflammatory_synergy(payload['biomarkers'], payload['condition_profile'], payload['available_kingdoms'])

start = time.perf_counter()
for _ in range(iterations):
    engine.predict_inflammatory_synergy(payload['biomarkers'], payload['condition_profile'], payload['available_kingdoms'])

elapsed = time.perf_counter() - start
out_lines = []
line1 = f"{iterations} runs completed in {elapsed:.4f}s — avg {elapsed/iterations*1000:.4f} ms/run"
print(line1)
out_lines.append(line1)

# Change input slightly to miss cache
payload2 = {
    "biomarkers": {"tnf_alpha": 13.5, "il6": 6.0, "crp": 3.1, "il1b": 0.8},
    "condition_profile": payload['condition_profile'],
    "available_kingdoms": payload['available_kingdoms']
}
start = time.perf_counter()
for _ in range(iterations):
    engine.predict_inflammatory_synergy(payload2['biomarkers'], payload2['condition_profile'], payload2['available_kingdoms'])
elapsed2 = time.perf_counter() - start
line2 = f"{iterations} cache-miss runs completed in {elapsed2:.4f}s — avg {elapsed2/iterations*1000:.4f} ms/run"
print(line2)
out_lines.append(line2)

# Write a small artifact for CI
with open('benchmark_out.txt', 'w') as fh:
    fh.write('\n'.join(out_lines) + '\n')