#!/usr/bin/env python3
"""Debug exact condition vs target mismatches."""
import json
from pathlib import Path

PROJECT_ROOT = Path(__file__).parent.parent
data_file = PROJECT_ROOT / 'data' / 'training' / 'neurobotanica_fully_enriched.json'

with open(data_file, 'r') as f:
    dataset = json.load(f)

compounds = dataset['compounds']

# Get first compound with clinical studies
for compound in compounds:
    if compound.get('clinical_studies'):
        print(f"Compound: {compound['compound_name']}")
        print(f"Therapeutic targets: {compound.get('therapeutic_targets', [])}")
        print()
        print("Clinical study conditions (first 10):")
        for i, study in enumerate(compound['clinical_studies'][:10]):
            print(f"  {i+1}. {study.get('condition')}")
        print()
        break
