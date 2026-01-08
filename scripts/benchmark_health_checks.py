"""Benchmark health monitor checks (force refresh) to find hotspots."""
import asyncio
import time
import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from backend.services import health_monitor

async def run_benchmark(iters: int = 100):
    # Ensure one immediate run to warm caches
    await health_monitor.get_cached_health(force_refresh=True)

    start = time.perf_counter()
    for _ in range(iters):
        await health_monitor.get_cached_health(force_refresh=True)
    elapsed = time.perf_counter() - start
    print(f"{iters} health refresh runs completed in {elapsed:.4f}s — avg {elapsed/iters*1000:.4f} ms/run")

if __name__ == '__main__':
    iters = int(os.getenv('HEALTH_BENCH_ITERS', '100'))
    asyncio.run(run_benchmark(iters))