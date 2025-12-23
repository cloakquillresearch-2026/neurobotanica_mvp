from typing import Dict

def validate_dimer_entry(entry: Dict) -> bool:
    required = ['compound_1', 'compound_2', 'dimer_type', 'interaction_effect', 'evidence_type']
    # Accept new evidence_type values
    valid_evidence_types = {'clinical', 'preclinical', 'in silico', 'prediction', 'omics', 'pathway', 'ai-prediction'}
    for field in required:
        if field not in entry or not entry[field]:
            return False
    if entry['evidence_type'] not in valid_evidence_types:
        return False
    # Optionally: check synergy_score is float if present
    if 'synergy_score' in entry and entry['synergy_score'] is not None:
        try:
            float(entry['synergy_score'])
        except (ValueError, TypeError):
            return False
    # Optionally: ai_predicted must be bool if present
    if 'ai_predicted' in entry and entry['ai_predicted'] is not None:
        if not isinstance(entry['ai_predicted'], bool):
            return False
    # Optionally: ai_features must be dict if present
    if 'ai_features' in entry and entry['ai_features'] is not None:
        if not isinstance(entry['ai_features'], dict):
            return False
    return True
