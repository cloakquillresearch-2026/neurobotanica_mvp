#!/usr/bin/env python3
"""
Import validated evidence JSON files into the main NeuroBotanica training dataset.
- Loads all *.json from data/clinical_evidence (excluding guide/query/report files)
- Merges with existing data/training/neurobotanica_enriched_dataset.json (if present)
- Writes merged output to data/training/neurobotanica_enriched_dataset.json
"""
import json
from pathlib import Path

EVIDENCE_DIR = Path(__file__).parent.parent / 'data' / 'clinical_evidence'
TRAINING_PATH = Path(__file__).parent.parent / 'data' / 'training' / 'neurobotanica_enriched_dataset.json'

# Load all evidence files
entries = []
for f in EVIDENCE_DIR.glob('*.json'):
    if f.name.startswith('EVIDENCE_') or f.name.endswith('report.json') or f.name.endswith('queries.txt'):
        continue
    try:
        entry = json.loads(f.read_text(encoding='utf-8'))
        entries.append(entry)
    except Exception as e:
        print(f"Skipping {f.name}: {e}")

print(f"Loaded {len(entries)} evidence entries from {EVIDENCE_DIR}")

# Load existing training data if present
if TRAINING_PATH.exists():
    with open(TRAINING_PATH, 'r', encoding='utf-8') as f:
        try:
            existing = json.load(f)
            print(f"Loaded {len(existing)} existing training entries.")
        except Exception as e:
            print(f"Error loading existing training data: {e}")
            existing = []
else:
    existing = []

# Only keep dict entries with a 'source' field
existing_dicts = [e for e in existing if isinstance(e, dict) and 'source' in e]
existing_sources = {e['source'] for e in existing_dicts}
new_entries = [e for e in entries if isinstance(e, dict) and e.get('source') not in existing_sources]
merged = existing_dicts + new_entries
print(f"Merged dataset: {len(merged)} entries (added {len(new_entries)})")

with open(TRAINING_PATH, 'w', encoding='utf-8') as f:
    json.dump(merged, f, indent=2)
print(f"Wrote merged training dataset to {TRAINING_PATH}")
