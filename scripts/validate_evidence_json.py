#!/usr/bin/env python3
"""
Validate generated evidence JSON files in data/clinical_evidence
Writes a summary report to data/clinical_evidence/validation_report.json
"""
import json
from pathlib import Path

EVIDENCE_DIR = Path(__file__).parent.parent / 'data' / 'clinical_evidence'
REPORT_PATH = EVIDENCE_DIR / 'validation_report.json'

required_keys = ['compound','condition','effect_size','study_type','source','year','notes','abstract']

files = list(EVIDENCE_DIR.glob('*.json'))
summary = {'checked_files': len(files), 'valid': 0, 'invalid': 0, 'errors': []}

for p in files:
    try:
        data = json.loads(p.read_text(encoding='utf-8'))
    except Exception as e:
        summary['invalid'] += 1
        summary['errors'].append({'file': p.name, 'error': f'json_load_error: {e}'})
        continue
    missing = [k for k in required_keys if k not in data]
    if missing:
        summary['invalid'] += 1
        summary['errors'].append({'file': p.name, 'error': 'missing_keys', 'missing': missing})
        continue
    # basic type checks
    if not isinstance(data.get('compound'), str) or not isinstance(data.get('condition'), str):
        summary['invalid'] += 1
        summary['errors'].append({'file': p.name, 'error': 'invalid_types_compound_condition'})
        continue
    # year can be int or None
    if data.get('year') is not None and not isinstance(data.get('year'), int):
        summary['invalid'] += 1
        summary['errors'].append({'file': p.name, 'error': 'invalid_year_type'})
        continue
    summary['valid'] += 1

REPORT_PATH.write_text(json.dumps(summary, indent=2), encoding='utf-8')
print(f"Checked {summary['checked_files']} files: {summary['valid']} valid, {summary['invalid']} invalid")
if summary['errors']:
    print('Errors written to', REPORT_PATH)
else:
    print('No errors found')
