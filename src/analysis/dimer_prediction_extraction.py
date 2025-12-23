import json
from typing import List
from src.analysis.dimer_schema import DimerEvidence
from src.analysis.dimer_validation import validate_dimer_entry

def convert_prediction_to_dimer_evidence(pred: dict) -> DimerEvidence:
    return DimerEvidence(
        compound_1=pred.get('monomer_a'),
        compound_2=pred.get('monomer_b'),
        dimer_type=pred.get('dimer_type'),
        interaction_effect='synergistic' if pred.get('scores', {}).get('synergy_score', 0) > 0.7 else 'unknown',
        evidence_type='prediction',
        effect_size=str(pred.get('scores', {}).get('synergy_score')),
        study_type=None,
        source='dimeric_predictions',
        confidence=str(pred.get('scores', {}).get('formation_probability')),
        notes=pred.get('dimer_name')
    )

def extract_dimeric_predictions(json_path: str) -> List[DimerEvidence]:
    with open(json_path, 'r') as f:
        data = json.load(f)
    predictions = data.get('predictions', [])
    dimers = []
    for pred in predictions:
        dimer = convert_prediction_to_dimer_evidence(pred)
        if validate_dimer_entry(dimer.__dict__):
            dimers.append(dimer)
    return dimers
