"""Check that critical Python imports succeed in the active environment.

This script is intended for CI and local dev to ensure dependencies are installable
and importable, preventing situations with broken/missing packages.

Usage: python scripts/check_imports.py
Exits with non-zero status code when any import fails.
"""
import sys

modules = [
    "fastapi",
    "annotated_doc",
    "anyio",
    "annotated_types",
    "redis",
]

failed = []
for m in modules:
    try:
        __import__(m)
        print(f"OK: imported {m}")
    except Exception as e:
        print(f"FAIL: could not import {m}: {e}")
        failed.append((m, str(e)))

# Also run pip check to surface dependency conflicts
print("Running 'pip check'...")
import subprocess
proc = subprocess.run([sys.executable, "-m", "pip", "check"], capture_output=True, text=True)
if proc.returncode != 0:
    print("pip check failed:\n", proc.stdout, proc.stderr)
    # treat pip check failures as fatal
    sys.exit(2)

if failed:
    print("\nImport validation failed for modules:")
    for m, e in failed:
        print(f" - {m}: {e}")
    sys.exit(1)

print("All imports OK")
sys.exit(0)
