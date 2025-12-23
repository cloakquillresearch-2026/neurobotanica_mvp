#!/usr/bin/env python3
"""
Fix Dataset: Add Therapeutic Targets from Clinical Studies
Merges study conditions into compound therapeutic_targets.
"""
import json
from pathlib import Path
from collections import Counter

PROJECT_ROOT = Path(__file__).parent.parent
input_file = PROJECT_ROOT / 'data' / 'training' / 'neurobotanica_fully_enriched.json'
output_file = PROJECT_ROOT / 'data' / 'training' / 'neurobotanica_fully_enriched_fixed.json'

with open(input_file, 'r') as f:
    dataset = json.load(f)

print("Fixing therapeutic targets based on clinical studies...")
print()

compounds_updated = 0
total_targets_added = 0

for compound in dataset['compounds']:
    clinical_studies = compound.get('clinical_studies', [])
    
    if not clinical_studies:
        continue
    
    # Get existing targets
    existing_targets = set(compound.get('therapeutic_targets', []))
    original_count = len(existing_targets)
    
    # Add conditions from clinical studies
    study_conditions = [s.get('condition') for s in clinical_studies if s.get('condition')]
    condition_counts = Counter(study_conditions)
    
    # Add conditions that appear in at least 2 studies (filter noise)
    for condition, count in condition_counts.items():
        if count >= 2:  # Threshold: at least 2 studies
            existing_targets.add(condition)
    
    # Update compound
    compound['therapeutic_targets'] = sorted(list(existing_targets))
    
    new_count = len(existing_targets)
    if new_count > original_count:
        compounds_updated += 1
        total_targets_added += (new_count - original_count)
        print(f"{compound['compound_name']}: {original_count} → {new_count} targets")

print()
print(f"Updated {compounds_updated} compounds")
print(f"Added {total_targets_added} new therapeutic targets from clinical evidence")
print()

# Save
with open(output_file, 'w') as f:
    json.dump(dataset, f, indent=2)

print(f"Saved: {output_file}")
