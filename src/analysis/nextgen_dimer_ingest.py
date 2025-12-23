import json
from typing import List, Dict
from src.analysis.nextgen_dimer_schema import NextGenDimerEvidence

def load_nextgen_dimer_evidence(json_path: str) -> List[NextGenDimerEvidence]:
    with open(json_path, 'r', encoding='utf-8') as f:
        data = json.load(f)
    dimers = []
    for entry in data:
        try:
            dimers.append(NextGenDimerEvidence(**entry))
        except Exception as e:
            print(f"Skipping invalid entry: {e}")
    return dimers

# Example usage for ingestion:
# nextgen_dimers = load_nextgen_dimer_evidence('data/processed/nextgen_dimers.json')
