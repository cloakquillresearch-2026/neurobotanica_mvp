#!/usr/bin/env python3
"""Scan the repository for forbidden XML-like artifacts and exit non-zero if any are found."""
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FORBIDDEN = [r"</content>", r"<parameter\s+name=", r"<filePath=", r"</parameter>"]
EXTS = {'.md', '.py', '.rst', '.txt', '.yaml', '.yml', '.json', '.toml', '.cfg'}

found = []
for dirpath, dirnames, filenames in os.walk(ROOT):
    # Skip common binary and env dirs and our own scripts
    if any(part in ('node_modules', '.venv', '.git', '__pycache__', '.pytest_cache', 'scripts') for part in dirpath.split(os.sep)):
        continue
    for fn in filenames:
        ext = os.path.splitext(fn)[1]
        if ext.lower() not in EXTS:
            continue
        path = os.path.join(dirpath, fn)
        # avoid scanning this script (it contains the FORBIDDEN literal strings)
        if os.path.abspath(path) == os.path.abspath(__file__):
            continue
        try:
            with open(path, 'r', encoding='utf-8') as f:
                for i, line in enumerate(f, start=1):
                    for pat in FORBIDDEN:
                        if re.search(pat, line):
                            found.append((path, i, line.strip()))
        except Exception:
            # skip binary or unreadable
            continue

if found:
    print('Forbidden XML-like artifacts detected:')
    for p, i, ln in found:
        print(f'{p}:{i}: {ln}')
    sys.exit(1)

print('No forbidden artifacts found.')
sys.exit(0)
