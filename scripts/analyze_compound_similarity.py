#!/usr/bin/env python3
"""
Improved Efficacy Estimation via Compound Similarity Transfer
Instead of defaulting to 0.5, infer efficacy from similar compounds with evidence.
"""
import json
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Dict, List
from sklearn.metrics.pairwise import cosine_similarity

PROJECT_ROOT = Path(__file__).parent.parent

# Load data
data_file = PROJECT_ROOT / 'data' / 'training' / 'neurobotanica_fully_enriched_fixed.json'
with open(data_file, 'r') as f:
    dataset = json.load(f)

print("=" * 70)
print("Improved Efficacy Estimation via Similarity Transfer")
print("=" * 70)
print()

# Extract compound features and efficacies
compounds_with_evidence = []
compounds_without_evidence = []

for compound in dataset['compounds']:
    descriptors = compound.get('rdkit_descriptors', {})
    if 'calculation_note' in descriptors:
        continue
    
    clinical_studies = compound.get('clinical_studies', [])
    has_studies = len(clinical_studies) > 0
    
    feature_vector = [
        descriptors.get('molecular_weight', 0),
        descriptors.get('logP', 0),
        descriptors.get('h_bond_donors', 0),
        descriptors.get('h_bond_acceptors', 0),
        descriptors.get('tpsa', 0),
        descriptors.get('complexity', 0),
        descriptors.get('qed', 0.5)
    ]
    
    compound_data = {
        'name': compound['compound_name'],
        'features': feature_vector,
        'targets': compound.get('therapeutic_targets', []),
        'studies': clinical_studies
    }
    
    if has_studies:
        compounds_with_evidence.append(compound_data)
    else:
        compounds_without_evidence.append(compound_data)

print(f"Compounds with clinical evidence: {len(compounds_with_evidence)}")
print(f"Compounds without evidence: {len(compounds_without_evidence)}")
print()

# Build similarity matrix
if len(compounds_with_evidence) > 0:
    evidence_features = np.array([c['features'] for c in compounds_with_evidence])
    no_evidence_features = np.array([c['features'] for c in compounds_without_evidence])
    
    # Compute cosine similarity
    similarities = cosine_similarity(no_evidence_features, evidence_features)
    
    print("Top 3 most similar evidence compounds for each no-evidence compound:")
    print()
    
    for i, compound in enumerate(compounds_without_evidence):
        top_3_idx = np.argsort(similarities[i])[-3:][::-1]
        
        print(f"{compound['name']}:")
        for idx in top_3_idx:
            similar_compound = compounds_with_evidence[idx]
            similarity_score = similarities[i][idx]
            num_studies = len(similar_compound['studies'])
            print(f"  → {similar_compound['name']}: {similarity_score:.3f} similarity ({num_studies} studies)")
        print()
    
    print("=" * 70)
    print("Recommendation:")
    print("  Transfer efficacy estimates from top-3 similar compounds")
    print("  weighted by similarity score and study quality.")
    print("=" * 70)
else:
    print("⚠️  No compounds with clinical evidence found!")
    print("    Cannot perform similarity transfer.")
