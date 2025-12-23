import json
from typing import List
from .dimer_schema import DimerEvidence
from .dimer_validation import validate_dimer_entry

def load_dimer_evidence(json_path: str) -> List[DimerEvidence]:
    with open(json_path, 'r') as f:
        data = json.load(f)
    dimers = []
    for entry in data:
        if validate_dimer_entry(entry):
            dimers.append(DimerEvidence(**entry))
    return dimers
