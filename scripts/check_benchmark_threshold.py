"""Check benchmark_out.txt and fail if average latency (ms) exceeds threshold.

Expected file format: two lines like:
"100 runs completed in 0.1234s — avg 0.0012 ms/run"
"100 cache-miss runs completed in 1.2345s — avg 0.0123 ms/run"
"""
import os
import re
import sys

THRESHOLD_MS = float(os.getenv('BENCH_THRESHOLD_MS', '5.0'))

if not os.path.exists('benchmark_out.txt'):
    print('No benchmark_out.txt found; skipping threshold check')
    sys.exit(0)

with open('benchmark_out.txt') as fh:
    lines = [l.strip() for l in fh if l.strip()]

avg_ms = None
for line in lines:
    m = re.search(r'avg\s+([0-9.]+)\s*ms/run', line)
    if m:
        avg_ms = float(m.group(1))
        print(f'Found avg ms: {avg_ms} (threshold {THRESHOLD_MS} ms)')
        if avg_ms > THRESHOLD_MS:
            print('Benchmark threshold exceeded — failing')
            sys.exit(2)

if avg_ms is None:
    print('Could not parse benchmark output; failing')
    sys.exit(1)

print('Benchmark within threshold')
sys.exit(0)
