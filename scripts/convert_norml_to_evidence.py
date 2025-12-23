#!/usr/bin/env python3
"""
Convert NORML extraction JSON into per-study evidence JSON files for review.
"""
import json
from pathlib import Path

ROOT = Path(__file__).parent.parent
NORML = ROOT / 'norml_complete_200plus_studies.json'
OUT = ROOT / 'data' / 'clinical_evidence'
OUT.mkdir(parents=True, exist_ok=True)

print(f"Loading NORML data from {NORML}")
with open(NORML, 'r', encoding='utf-8') as f:
    data = json.load(f)

# support both top-level list and NORML object with `studies` key
studies = data if isinstance(data, list) else data.get('studies', [])

print(f"Found {len(studies)} studies in NORML data")

count = 0
for item in studies:
    pmid = item.get('study_id') or item.get('pmid') or item.get('id') or f"norml_{count}"
    entry = {
        'compound': item.get('intervention', {}).get('cannabinoid', 'UNKNOWN') if 'intervention' in item else item.get('compound', 'UNKNOWN'),
        'condition': item.get('condition', item.get('indication', 'UNKNOWN')),
        'effect_size': item.get('outcomes', {}).get('effect_size', 'none') if 'outcomes' in item else item.get('effect_size', 'none'),
        'study_type': item.get('study_type', 'clinical_study'),
        'source': f"NORML:{pmid}",
        'participants': item.get('sample_size'),
        'year': item.get('year'),
        'notes': item.get('study_title', item.get('title', '')),
        'confidence': item.get('confidence', 'medium'),
        'abstract': item.get('outcomes', {}).get('results', '') if 'outcomes' in item else item.get('abstract', '')
    }
    fname = f"norml_{pmid}.json"
    with open(OUT / fname, 'w', encoding='utf-8') as of:
        json.dump(entry, of, indent=2)
    count += 1
    if count % 50 == 0:
        print(f"Processed {count} studies...")

print(f"Converted {count} NORML studies into {OUT}")
