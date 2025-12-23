#!/usr/bin/env python3
"""
All-in-one NORML enrichment and validation script.
- Loads norml_complete_200plus_studies.json
- Writes per-study JSON files to data/clinical_evidence/
- Validates all outputs and writes a summary report
- Prints summary to stdout and logs
"""
import json
from pathlib import Path

ROOT = Path(__file__).parent.parent
NORML = ROOT / 'norml_complete_200plus_studies.json'
OUT = ROOT / 'data' / 'clinical_evidence'
OUT.mkdir(parents=True, exist_ok=True)
LOG = ROOT / 'logs' / 'norml_enrich_and_validate.log'
LOG.parent.mkdir(parents=True, exist_ok=True)

print(f"Loading NORML data from {NORML}")
with open(NORML, 'r', encoding='utf-8') as f:
    data = json.load(f)

studies = data if isinstance(data, list) else data.get('studies', [])
print(f"Found {len(studies)} studies in NORML data")

required_keys = ['compound','condition','effect_size','study_type','source','year','notes','abstract']
valid = 0
invalid = 0
errors = []
count = 0

for item in studies:
    pmid = item.get('study_id') or item.get('pmid') or item.get('id') or f"norml_{count}"
    # Robustly extract compound from intervention (dict, str, None, or other)
    try:
        intervention = item.get('intervention')
        if isinstance(intervention, dict):
            compound = intervention.get('cannabinoid', 'UNKNOWN')
        elif isinstance(intervention, str):
            compound = intervention
        elif intervention is None:
            compound = item.get('compound', 'UNKNOWN')
        else:
            compound = str(intervention)
        entry = {
            'compound': compound,
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
        # validation
        missing = [k for k in required_keys if k not in entry]
        if missing:
            invalid += 1
            errors.append({'file': fname, 'error': 'missing_keys', 'missing': missing})
        elif not isinstance(entry.get('compound'), str) or not isinstance(entry.get('condition'), str):
            invalid += 1
            errors.append({'file': fname, 'error': 'invalid_types_compound_condition'})
        elif entry.get('year') is not None and not isinstance(entry.get('year'), int):
            invalid += 1
            errors.append({'file': fname, 'error': 'invalid_year_type'})
        else:
            valid += 1
    except Exception as e:
        invalid += 1
        errors.append({'file': f'norml_{pmid}.json', 'error': f'processing_error: {e}'})
    count += 1
    if count % 50 == 0:
        print(f"Processed {count} studies...")

summary = {
    'checked_files': count,
    'valid': valid,
    'invalid': invalid,
    'errors': errors[:10]  # show only first 10 errors for brevity
}

report_path = OUT / 'validation_report.json'
with open(report_path, 'w', encoding='utf-8') as rf:
    json.dump(summary, rf, indent=2)

with open(LOG, 'w', encoding='utf-8') as lf:
    lf.write(json.dumps(summary, indent=2))

print(f"\nChecked {count} files: {valid} valid, {invalid} invalid")
print(f"Validation report written to {report_path}")
if errors:
    print(f"First error: {errors[0]}")
else:
    print('No errors found')
