#!/usr/bin/env python3
"""
Generate Synthetic Patient Response Data
NeuroBotanica Data Enrichment Step 3/4

Creates anonymized patient response profiles from clinical study outcomes.
HIPAA-compliant synthetic data generation.
"""
import json
import random
from pathlib import Path
from typing import Dict, List
from datetime import datetime, timedelta


# Patient demographic distributions
AGE_RANGES = ["18-25", "26-35", "36-45", "46-55", "56-65", "66+"]
WEIGHT_CATEGORIES = ["underweight", "normal", "overweight", "obese"]
EXPERIENCE_LEVELS = ["naive", "occasional", "regular", "daily"]
THC_TOLERANCE = ["low", "medium", "high"]
DELIVERY_METHODS = ["inhalation", "sublingual", "edible", "topical", "transdermal"]

# Response ratings (1-10 scale)
RESPONSE_CATEGORIES = {
    'excellent': (8, 10),
    'good': (6, 7),
    'moderate': (4, 5),
    'minimal': (2, 3),
    'none': (0, 1)
}

# Effect size to response mapping
EFFECT_SIZE_TO_RESPONSE = {
    'large': 'excellent',
    'medium': 'good',
    'small': 'moderate',
    None: 'minimal'
}


def generate_patient_profiles(study: Dict, num_patients: int = 5) -> List[Dict]:
    """Generate synthetic patient profiles from study outcomes."""
    profiles = []
    
    # Extract study details
    condition = study.get('condition', 'UNKNOWN')
    
    # Handle outcomes field (can be dict or list)
    outcomes = study.get('outcomes', {})
    if not isinstance(outcomes, dict):
        outcomes = {}
    
    effect_size = outcomes.get('effect_size_category')
    key_findings = outcomes.get('key_findings', [])
    
    # Determine response distribution based on effect size
    response_category = EFFECT_SIZE_TO_RESPONSE.get(effect_size, 'minimal')
    
    for i in range(num_patients):
        # Handle intervention field (can be dict or string)
        intervention = study.get('intervention', {})
        if not isinstance(intervention, dict):
            intervention = {}
        
        # Generate demographic profile
        profile = {
            'patient_id': f"SYNTH_{study.get('study_id', 'UNKNOWN')}_{i+1}",
            'demographics': {
                'age_range': random.choice(AGE_RANGES),
                'weight_category': random.choice(WEIGHT_CATEGORIES),
                'experience_level': random.choice(EXPERIENCE_LEVELS),
                'thc_tolerance': random.choice(THC_TOLERANCE)
            },
            'condition': condition,
            'treatment': {
                'study_source': study.get('study_id'),
                'cannabinoid_profile': intervention.get('cannabinoid_profile', 'mixed'),
                'delivery_method': _extract_delivery_method(study),
                'duration_days': _estimate_duration(study)
            },
            'response': {
                'overall_efficacy_score': _generate_response_score(response_category),
                'symptom_improvement': _generate_symptom_scores(response_category),
                'adverse_effects': _generate_adverse_effects(),
                'satisfaction_score': _generate_response_score(response_category),
                'adherence_score': random.randint(7, 10)
            },
            'outcome_notes': _select_finding(key_findings) if key_findings else None,
            'data_quality': {
                'synthetic': True,
                'source_study_quality': _assess_study_quality(study),
                'generated_date': datetime.now().isoformat()
            }
        }
        
        profiles.append(profile)
    
    return profiles


def _extract_delivery_method(study: Dict) -> str:
    """Extract delivery method from study intervention."""
    intervention = study.get('intervention', {})
    if not isinstance(intervention, dict):
        intervention = {}
    
    delivery = intervention.get('delivery_method', '')
    
    if 'inhal' in delivery.lower() or 'smok' in delivery.lower():
        return 'inhalation'
    elif 'oral' in delivery.lower() or 'capsul' in delivery.lower():
        return 'edible'
    elif 'sublingual' in delivery.lower():
        return 'sublingual'
    elif 'topical' in delivery.lower():
        return 'topical'
    elif 'transdermal' in delivery.lower():
        return 'transdermal'
    else:
        return random.choice(DELIVERY_METHODS)


def _estimate_duration(study: Dict) -> int:
    """Estimate treatment duration in days."""
    intervention = study.get('intervention', {})
    if not isinstance(intervention, dict):
        intervention = {}
    
    duration_text = intervention.get('treatment_duration', '').lower()
    
    if 'acute' in duration_text or 'single' in duration_text:
        return 1
    elif 'week' in duration_text:
        weeks = 4  # default
        if '12' in duration_text:
            weeks = 12
        elif '8' in duration_text:
            weeks = 8
        return weeks * 7
    elif 'month' in duration_text:
        months = 3  # default
        if '6' in duration_text:
            months = 6
        return months * 30
    else:
        return random.randint(14, 84)  # 2-12 weeks


def _generate_response_score(category: str) -> int:
    """Generate response score within category range."""
    min_score, max_score = RESPONSE_CATEGORIES[category]
    return random.randint(min_score, max_score)


def _generate_symptom_scores(response_category: str) -> Dict[str, int]:
    """Generate symptom-specific improvement scores."""
    base_score = _generate_response_score(response_category)
    
    # Add some variance to specific symptoms
    return {
        'primary_symptom': base_score,
        'secondary_symptoms': max(0, base_score + random.randint(-2, 1)),
        'quality_of_life': max(0, base_score + random.randint(-1, 2)),
        'functional_improvement': max(0, base_score + random.randint(-2, 1))
    }


def _generate_adverse_effects() -> List[str]:
    """Generate realistic adverse effect profile."""
    common_effects = [
        'dry_mouth', 'dizziness', 'drowsiness', 'increased_appetite',
        'mild_anxiety', 'headache', 'fatigue'
    ]
    
    # Most patients have 0-2 mild side effects
    num_effects = random.choices([0, 1, 2, 3], weights=[30, 40, 20, 10])[0]
    
    if num_effects == 0:
        return []
    
    return random.sample(common_effects, num_effects)


def _select_finding(findings: List[str]) -> str:
    """Select a representative finding."""
    if not findings:
        return None
    return random.choice(findings)


def _assess_study_quality(study: Dict) -> str:
    """Assess study quality tier."""
    quality = study.get('study_quality', {})
    study_type = study.get('study_type', '')
    
    if study_type == 'RCT' and quality.get('randomized') and quality.get('placebo_controlled'):
        return 'high'
    elif study_type == 'RCT':
        return 'good'
    elif study_type in ['COHORT', 'OBSERVATIONAL']:
        return 'moderate'
    else:
        return 'low'


def integrate_patient_data(dataset: Dict, clinical_studies: Dict) -> Dict:
    """Generate and integrate synthetic patient response data."""
    enriched = dataset.copy()
    
    all_patient_profiles = []
    stats = {
        'total_patients': 0,
        'by_condition': {},
        'by_response_category': {cat: 0 for cat in RESPONSE_CATEGORIES.keys()}
    }
    
    # Generate patients from clinical studies
    for study in clinical_studies.get('studies', []):
        # Generate 3-5 patients per high-quality study, 1-2 for lower quality
        quality = _assess_study_quality(study)
        num_patients = {
            'high': random.randint(4, 6),
            'good': random.randint(3, 5),
            'moderate': random.randint(2, 3),
            'low': random.randint(1, 2)
        }.get(quality, 2)
        
        profiles = generate_patient_profiles(study, num_patients)
        all_patient_profiles.extend(profiles)
        
        # Update stats
        condition = study.get('condition', 'UNKNOWN')
        stats['by_condition'][condition] = stats['by_condition'].get(condition, 0) + num_patients
        stats['total_patients'] += num_patients
    
    enriched['patient_response_data'] = {
        'profiles': all_patient_profiles,
        'metadata': {
            'total_patients': stats['total_patients'],
            'data_type': 'synthetic',
            'hipaa_compliant': True,
            'source': 'NORML clinical studies',
            'generation_date': datetime.now().isoformat(),
            'statistics': stats
        }
    }
    
    enriched['enrichment_metadata'] = enriched.get('enrichment_metadata', {})
    enriched['enrichment_metadata']['patient_data_integration'] = {
        'date': '2025-12-23',
        'total_patients': stats['total_patients'],
        'conditions_represented': len(stats['by_condition'])
    }
    
    return enriched


def main():
    """Main integration workflow."""
    print("=" * 70)
    print("Synthetic Patient Response Data Generation")
    print("=" * 70)
    print()
    
    # Input/output paths
    SCRIPT_DIR = Path(__file__).parent
    PROJECT_ROOT = SCRIPT_DIR.parent
    
    input_file = PROJECT_ROOT / 'data' / 'training' / 'neurobotanica_enriched_synergies.json'
    clinical_file = PROJECT_ROOT / 'norml_complete_200plus_studies.json'
    output_file = PROJECT_ROOT / 'data' / 'training' / 'neurobotanica_enriched_patients.json'
    
    # Load enriched dataset from step 2
    print(f"Loading enriched dataset: {input_file}")
    with open(input_file, 'r', encoding='utf-8') as f:
        dataset = json.load(f)
    print(f"✅ Loaded {len(dataset['compounds'])} compounds")
    print()
    
    # Load clinical studies
    print(f"Loading clinical studies: {clinical_file}")
    with open(clinical_file, 'r', encoding='utf-8') as f:
        clinical_data = json.load(f)
    print(f"✅ Loaded {len(clinical_data.get('studies', []))} studies")
    print()
    
    # Generate patient data
    print("Generating synthetic patient response profiles...")
    enriched = integrate_patient_data(dataset, clinical_data)
    
    metadata = enriched['patient_response_data']['metadata']
    print(f"   Total synthetic patients: {metadata['total_patients']}")
    print(f"   Conditions represented: {metadata['statistics']['total_patients']}")
    print()
    
    # Save
    print(f"Saving: {output_file}")
    with open(output_file, 'w', encoding='utf-8') as f:
        json.dump(enriched, f, indent=2)
    
    print("=" * 70)
    print("✅ Patient Data Integration Complete!")
    print(f"✅ Output: {output_file}")
    print("✅ HIPAA-compliant synthetic data")
    print("=" * 70)


if __name__ == '__main__':
    main()
