from typing import Dict

def validate_dimer_entry(entry: Dict) -> bool:
    required = ['compound_1', 'compound_2', 'dimer_type', 'interaction_effect', 'evidence_type']
    for field in required:
        if field not in entry or not entry[field]:
            return False
    return True
