import json
from collections import defaultdict

filepath = "data/training/neurobotanica_enriched_with_dimers.json"

with open(filepath, "r") as f:
    data = json.load(f)

print(f"Total records loaded: {len(data)}")

isomer_pairs = [
    {"THC", "dronabinol"},
    {"Delta8-THC", "THC"},
    {"Delta8-THC", "dronabinol"},
]

def extract_compounds(compound_str):
    if not compound_str:
        return []
    for sep in ["+", "/", " and ", " & "]:
        if sep in compound_str:
            return [c.strip() for c in compound_str.split(sep)]
    return [compound_str.strip()]

self_pairs = []
isomer_hits = []
seen_pairs = defaultdict(list)

for i, e in enumerate(data):
    compound = e.get("compound", "")
    parts = extract_compounds(compound)

    if len(parts) >= 2 and len(set(parts)) == 1:
        self_pairs.append((i, compound))

    if len(parts) >= 2:
        pair_set = set(parts[:2])
        for iso in isomer_pairs:
            if iso.issubset(pair_set) or pair_set == iso:
                isomer_hits.append((i, compound, e.get("source")))

    key = frozenset(parts[:2]) if len(parts) >= 2 else compound
    seen_pairs[key].append(i)

duplicates = {k: v for k, v in seen_pairs.items() if len(v) > 1 and isinstance(k, frozenset)}

print(f"\nSelf-pairs found: {len(self_pairs)}")
for i, c in self_pairs:
    print(f"  [{i}] {c}")

print(f"\nStructural isomer pairs found: {len(isomer_hits)}")
for i, c, src in isomer_hits:
    print(f"  [{i}] {c} | source: {src}")

print(f"\nDuplicate compound pairs: {len(duplicates)}")
for k, indices in list(duplicates.items())[:20]:
    print(f"  {sorted(list(k))} → {len(indices)} records")
import json
from collections import defaultdict

filepath = "data/training/neurobotanica_enriched_with_dimers.json"

with open(filepath, "r") as f:
    data = json.load(f)

print(f"Total records loaded: {len(data)}")

# Known structural isomer pairs
isomer_pairs = [
    {"THC", "dronabinol"},
    {"Delta8-THC", "THC"},
    {"Delta8-THC", "dronabinol"},
]

# Normalize compound names for matching
def extract_compounds(compound_str):
    """Try to extract individual compounds from combo strings like 'THC+CBD'"""
    if not compound_str:
        return []
    for sep in ["+", "/", " and ", " & "]:
        if sep in compound_str:
            return [c.strip() for c in compound_str.split(sep)]
    return [compound_str.strip()]

self_pairs = []
isomer_hits = []
seen_pairs = defaultdict(list)

for i, e in enumerate(data):
    compound = e.get("compound", "")
    parts = extract_compounds(compound)

    # Self-pairs
    if len(parts) >= 2 and len(set(parts)) == 1:
        self_pairs.append((i, compound))

    # Structural isomers
    if len(parts) >= 2:
        pair_set = set(parts[:2])
        for iso in isomer_pairs:
            if iso.issubset(pair_set) or pair_set == iso:
                isomer_hits.append((i, compound, e.get("source")))

    # Duplicate detection
    key = frozenset(parts[:2]) if len(parts) >= 2 else compound
    seen_pairs[key].append(i)

duplicates = {k: v for k, v in seen_pairs.items() if len(v) > 1 and len(k) > 1}

print(f"\nSelf-pairs found: {len(self_pairs)}")
for i, c in self_pairs:
    print(f"  [{i}] {c}")

print(f"\nStructural isomer pairs found: {len(isomer_hits)}")
for i, c, src in isomer_hits:
    print(f"  [{i}] {c} | source: {src}")

print(f"\nDuplicate compound pairs (both directions): {len(duplicates)}")
for k, indices in list(duplicates.items())[:20]:
    print(f"  {list(k)} → {len(indices)} records at indices {indices[:5]}")
