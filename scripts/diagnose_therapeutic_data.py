#!/usr/bin/env python3
"""
Diagnose TherapeuticPredictionModel Data Issues
"""
import json
import pandas as pd
import numpy as np
from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(PROJECT_ROOT / 'backend'))

# Load dataset
data_file = PROJECT_ROOT / 'data' / 'training' / 'neurobotanica_fully_enriched.json'
with open(data_file, 'r') as f:
    dataset = json.load(f)

print("=" * 70)
print("TherapeuticPredictionModel Data Diagnostics")
print("=" * 70)
print()

# Check clinical studies integration
compounds = dataset['compounds']
compounds_with_studies = [c for c in compounds if c.get('clinical_studies')]

print(f"Total compounds: {len(compounds)}")
print(f"Compounds with clinical_studies field: {len(compounds_with_studies)}")
print()

if compounds_with_studies:
    print("Sample compound with clinical studies:")
    sample = compounds_with_studies[0]
    print(f"  Name: {sample.get('compound_name')}")
    print(f"  Clinical studies: {len(sample.get('clinical_studies', []))}")
    print(f"  First study: {sample.get('clinical_studies', [{}])[0]}")
    print()

# Check therapeutic targets distribution
all_targets = []
for compound in compounds:
    targets = compound.get('therapeutic_targets', [])
    all_targets.extend(targets)

print(f"Total therapeutic target entries: {len(all_targets)}")
print(f"Unique targets: {len(set(all_targets))}")
print()

# Check efficacy score distribution
print("Checking efficacy estimation logic...")

def _estimate_efficacy(compound, target, studies):
    """Same logic as training script - UPDATED with fuzzy matching."""
    # Normalize target name for matching
    target_normalized = target.upper().replace(' ', '_').replace('-', '_')
    
    # Find relevant studies with fuzzy matching
    relevant_studies = []
    for s in studies:
        condition = s.get('condition', '').upper().replace(' ', '_').replace('-', '_')
        # Check if target matches condition or is a substring
        if (target_normalized in condition or 
            condition in target_normalized or
            target_normalized.replace('_', '') in condition.replace('_', '')):
            relevant_studies.append(s)
    
    if not relevant_studies:
        return 0.5
    
    effect_map = {
        'large': 0.85,
        'medium': 0.70,
        'small': 0.55,
        'minimal': 0.40,
        None: 0.50
    }
    
    efficacies = [effect_map.get(s.get('effect_size'), 0.5) for s in relevant_studies]
    return np.mean(efficacies)

efficacy_scores = []
targets_with_evidence = 0

for compound in compounds:
    if 'calculation_note' in compound.get('rdkit_descriptors', {}):
        continue  # Skip estimated descriptors
    
    clinical_studies = compound.get('clinical_studies', [])
    therapeutic_targets = compound.get('therapeutic_targets', [])
    
    for target in therapeutic_targets:
        efficacy = _estimate_efficacy(compound, target, clinical_studies)
        efficacy_scores.append(efficacy)
        
        if efficacy != 0.5:
            targets_with_evidence += 1

print(f"Total efficacy scores: {len(efficacy_scores)}")
print(f"Targets with direct evidence (≠0.5): {targets_with_evidence}")
print(f"Targets with default score (=0.5): {len(efficacy_scores) - targets_with_evidence}")
print()

if efficacy_scores:
    print(f"Efficacy score distribution:")
    print(f"  Mean: {np.mean(efficacy_scores):.3f}")
    print(f"  Std: {np.std(efficacy_scores):.3f}")
    print(f"  Min: {np.min(efficacy_scores):.3f}")
    print(f"  Max: {np.max(efficacy_scores):.3f}")
    print(f"  Unique values: {sorted(set(efficacy_scores))}")
    print()

# Check feature variance
print("Checking feature variance...")
features_list = []

for compound in compounds:
    descriptors = compound.get('rdkit_descriptors', {})
    if 'calculation_note' in descriptors:
        continue
    
    therapeutic_targets = compound.get('therapeutic_targets', [])
    if not therapeutic_targets:
        continue
    
    for target in therapeutic_targets:
        features = {
            'molecular_weight': descriptors.get('molecular_weight', 0),
            'logP': descriptors.get('logP', 0),
            'complexity': descriptors.get('complexity', 0),
            'target_encoded': hash(target) % 100
        }
        features_list.append(features)

if features_list:
    df = pd.DataFrame(features_list)
    print(f"Feature variance:")
    for col in df.columns:
        print(f"  {col}: {df[col].var():.6f} (std={df[col].std():.4f})")
    print()

print("=" * 70)
print("DIAGNOSIS:")

if targets_with_evidence == 0:
    print("❌ CRITICAL: No therapeutic targets have clinical evidence!")
    print("   All efficacy scores are defaulting to 0.5 (no signal)")
    print()
    print("ROOT CAUSE:")
    print("   The 'clinical_studies' field is not matching 'therapeutic_targets'")
    print("   because condition names don't align.")
    print()
    print("SOLUTION:")
    print("   Need to map NORML conditions to therapeutic_targets:")
    print("   - 'CHRONIC_PAIN' → 'chronic_pain'")
    print("   - 'ANXIETY' → 'anxiety'")
    print("   - 'EPILEPSY' → 'epilepsy'")
    print("   etc.")
elif targets_with_evidence < len(efficacy_scores) * 0.1:
    print("⚠️  WARNING: <10% of targets have clinical evidence")
    print(f"   Only {targets_with_evidence}/{len(efficacy_scores)} targets matched to studies")
else:
    print("✅ Data quality appears OK")
    print(f"   {targets_with_evidence}/{len(efficacy_scores)} targets have evidence")

print("=" * 70)
