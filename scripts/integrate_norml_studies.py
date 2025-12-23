"""
Integrate NORML Clinical Studies into NeuroBotanica Training Dataset
=====================================================================

Combines 320+ clinical studies from NORML extraction with the 63-compound
dataset to create enriched training data for therapeutic prediction models.

Input:
- data/training/neurobotanica_complete_dataset_63compounds.json
- data/norml_extraction/*.json (22 condition-specific files)

Output:
- data/training/neurobotanica_enriched_dataset.json
"""

import json
from pathlib import Path
from typing import Dict, List
from collections import defaultdict
import sys

def load_norml_studies(norml_dir: Path) -> Dict:
    """Load all NORML clinical studies."""
    print(f"Loading NORML studies from {norml_dir}...")
    
    all_studies = []
    condition_counts = {}
    
    # Get all JSON files
    json_files = list(norml_dir.glob("*_studies.json"))
    
    for json_file in json_files:
        with open(json_file, 'r', encoding='utf-8') as f:
            data = json.load(f)
            condition = data.get('condition', json_file.stem.replace('_studies', ''))
            studies = data.get('studies', [])
            
            all_studies.extend(studies)
            condition_counts[condition] = len(studies)
            
    print(f"✅ Loaded {len(all_studies)} studies across {len(condition_counts)} conditions")
    for condition, count in sorted(condition_counts.items()):
        print(f"   {condition}: {count} studies")
    
    return {
        'total_studies': len(all_studies),
        'studies': all_studies,
        'condition_counts': condition_counts
    }


def map_studies_to_compounds(compounds: List[Dict], studies: List[Dict]) -> Dict:
    """Map clinical studies to relevant cannabinoid compounds."""
    
    # Build cannabinoid-to-study mappings
    compound_studies = defaultdict(list)
    
    # Keywords for matching cannabinoids mentioned in studies
    cannabinoid_keywords = {
        'THC': ['thc', 'delta-9', 'tetrahydrocannabinol', 'dronabinol', 'marinol'],
        'CBD': ['cbd', 'cannabidiol', 'epidiolex'],
        'CBG': ['cbg', 'cannabigerol'],
        'CBN': ['cbn', 'cannabinol'],
        'CBC': ['cbc', 'cannabichromene'],
        'THCA': ['thca', 'tetrahydrocannabinolic acid'],
        'CBDA': ['cbda', 'cannabidiolic acid'],
        'THCV': ['thcv', 'tetrahydrocannabivarin'],
        'CBDV': ['cbdv', 'cannabidivarin']
    }
    
    for study in studies:
        study_text = json.dumps(study).lower()
        matched_cannabinoids = []
        
        # Check which cannabinoids are mentioned
        for compound_name, keywords in cannabinoid_keywords.items():
            if any(kw in study_text for kw in keywords):
                matched_cannabinoids.append(compound_name)
        
        # If no specific cannabinoids found, assume whole-plant cannabis (THC+CBD)
        if not matched_cannabinoids:
            if 'cannabis' in study_text or 'marijuana' in study_text:
                matched_cannabinoids = ['THC', 'CBD']
        
        # Add study to matched compounds
        for compound_name in matched_cannabinoids:
            # Handle outcomes field (can be dict or missing)
            outcomes = study.get('outcomes', {})
            if not isinstance(outcomes, dict):
                outcomes = {}
            
            compound_studies[compound_name].append({
                'study_id': study.get('study_id'),
                'condition': study.get('condition'),
                'study_type': study.get('study_type'),
                'year': study.get('publication_year') or study.get('year'),
                'key_findings': outcomes.get('key_findings', []),
                'effect_size': outcomes.get('effect_size_category')
            })
    
    return dict(compound_studies)


def enrich_dataset(base_dataset: Dict, norml_data: Dict) -> Dict:
    """Enrich base dataset with NORML clinical studies."""
    
    print("\nEnriching dataset with clinical studies...")
    
    # Map studies to compounds
    compound_studies = map_studies_to_compounds(
        base_dataset['compounds'],
        norml_data['studies']
    )
    
    # Add clinical_studies field to each compound
    enriched_compounds = []
    for compound in base_dataset['compounds']:
        compound_copy = compound.copy()
        compound_name = compound['compound_name']
        
        # Add clinical studies if available
        if compound_name in compound_studies:
            compound_copy['clinical_studies'] = compound_studies[compound_name]
            print(f"   {compound_name}: {len(compound_studies[compound_name])} studies")
        else:
            compound_copy['clinical_studies'] = []
        
        enriched_compounds.append(compound_copy)
    
    # Create enriched dataset
    enriched = {
        'metadata': {
            **base_dataset['metadata'],
            'enrichment_version': '1.0',
            'enrichment_date': '2025-12-23',
            'norml_studies_integrated': norml_data['total_studies'],
            'conditions_covered': list(norml_data['condition_counts'].keys())
        },
        'compounds': enriched_compounds,
        'norml_summary': {
            'total_studies': norml_data['total_studies'],
            'condition_counts': norml_data['condition_counts']
        }
    }
    
    return enriched


def main():
    """Main integration workflow."""
    print("=" * 70)
    print("NORML Clinical Studies Integration")
    print("=" * 70)
    
    project_root = Path(__file__).parent.parent
    
    # Load base dataset
    base_path = project_root / 'data/training/neurobotanica_complete_dataset_63compounds.json'
    print(f"\nLoading base dataset: {base_path}")
    with open(base_path, 'r', encoding='utf-8') as f:
        base_dataset = json.load(f)
    print(f"✅ Loaded {len(base_dataset['compounds'])} compounds")
    
    # Load NORML studies
    norml_dir = project_root / 'data/norml_extraction'
    norml_data = load_norml_studies(norml_dir)
    
    # Enrich dataset
    enriched_dataset = enrich_dataset(base_dataset, norml_data)
    
    # Save enriched dataset
    output_path = project_root / 'data/training/neurobotanica_enriched_dataset.json'
    print(f"\nSaving enriched dataset: {output_path}")
    with open(output_path, 'w', encoding='utf-8') as f:
        json.dump(enriched_dataset, f, indent=2)
    
    print("\n" + "=" * 70)
    print("✅ Integration Complete!")
    print(f"✅ Output: {output_path}")
    print(f"✅ Total compounds: {len(enriched_dataset['compounds'])}")
    print(f"✅ Total studies integrated: {norml_data['total_studies']}")
    print("=" * 70)


if __name__ == "__main__":
    main()
