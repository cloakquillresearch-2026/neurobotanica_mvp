import json
from typing import List, Dict
from src.analysis.dimer_schema import DimerEvidence
from src.analysis.dimer_validation import validate_dimer_entry

# Helper: Map AI-predicted fields to schema
AI_FIELDS = [
    'enhanced_lipophilicity', 'increased_receptor_affinity', 'prolonged_duration',
    'reduced_metabolism', 'therapeutic_potential_score', 'confidence_score'
]

def extract_dimeric_predictions(json_path: str) -> List[DimerEvidence]:
    """Extract and validate AI-predicted dimer evidence from neurobotanica_dimeric_predictions.json."""
    with open(json_path, 'r') as f:
        data = json.load(f)
    dimers = []
    for entry in data.get('dimeric_cannabinoids', []):
        # Map AI fields
        ai_features = {k: entry[k] for k in AI_FIELDS if k in entry}
        dimer_entry = {
            'compound_1': entry['compound_1'],
            'compound_2': entry['compound_2'],
            'dimer_type': entry.get('dimer_type', 'unknown'),
            'interaction_effect': 'unknown',
            'evidence_type': 'ai-prediction',
            'ai_predicted': True,
            'ai_features': ai_features,
            'confidence': str(entry.get('confidence_score')) if 'confidence_score' in entry else None
        }
        if validate_dimer_entry(dimer_entry):
            dimers.append(DimerEvidence(**dimer_entry))
    return dimers

def extract_synergy_evidence(json_path: str) -> List[DimerEvidence]:
    """Extract and validate synergy evidence from neurobotanica_enriched_synergies.json."""
    with open(json_path, 'r') as f:
        data = json.load(f)
    dimers = []
    for study in data:
        if 'terpene_synergies' in study:
            for synergy in study['terpene_synergies']:
                dimer_entry = {
                    'compound_1': study.get('compound', 'unknown'),
                    'compound_2': synergy.get('partner_compound', 'unknown'),
                    'dimer_type': f"cannabinoid-{synergy.get('partner_type', 'unknown')}",
                    'interaction_effect': 'synergistic',
                    'evidence_type': 'clinical',
                    'effect_size': str(synergy.get('enhancement_factor')) if 'enhancement_factor' in synergy else None,
                    'source': study.get('source'),
                    'notes': synergy.get('mechanism'),
                    'confidence': synergy.get('evidence_quality'),
                    'synergy_score': float(synergy['enhancement_factor']) if 'enhancement_factor' in synergy else None
                }
                if validate_dimer_entry(dimer_entry):
                    dimers.append(DimerEvidence(**dimer_entry))
    return dimers

if __name__ == "__main__":
    ai_dimers = extract_dimeric_predictions("data/training/neurobotanica_dimeric_predictions.json")
    synergy_dimers = extract_synergy_evidence("data/training/neurobotanica_enriched_synergies.json")
    print(f"AI-predicted dimers: {len(ai_dimers)}")
    print(f"Synergy dimers: {len(synergy_dimers)}")
    # Optionally: Save or further process
