import json
import os

def validate_study_extraction(filename):
    """Quick validation of extracted studies"""
    with open(filename, encoding='utf-8') as f:
        data = json.load(f)
    
    studies = data['studies']
    
    print(f"\n{'='*60}")
    print(f"VALIDATION: {data['condition']}")
    print(f"{'='*60}")
    print(f"Total studies: {len(studies)}")
    print(f"Phase: {data['phase']}")
    print(f"Extraction date: {data['extraction_date']}")
    
    # Check required fields
    required_fields = ['study_id', 'study_type', 'condition', 'study_title', 'citation']
    
    complete_count = 0
    for i, study in enumerate(studies, 1):
        missing = [f for f in required_fields if f not in study]
        if missing:
            print(f"⚠️  Study {i}: Missing fields: {missing}")
        else:
            print(f"✓ Study {i}: {study['study_id']} - Complete")
            complete_count += 1
    
    print(f"{'='*60}")
    print(f"Complete studies: {complete_count}/{len(studies)}")
    print(f"Completion rate: {(complete_count/len(studies)*100):.1f}%")
    print(f"{'='*60}\n")

def merge_phase_data(output_filename='norml_complete_200plus_studies.json'):
    """Merge Phase 1, Phase 2, and Phase 3 data"""
    
    # Load existing Phase 1 legacy data (if exists)
    phase1_legacy_file = 'data/norml_extraction/norml_training_dataset_FINAL.json'
    
    # Load Phase 1 individual condition files (NEW)
    phase1_files = [
        'data/norml_extraction/ptsd_studies.json',
        'data/norml_extraction/epilepsy_studies.json',
        'data/norml_extraction/insomnia_studies.json',
        'data/norml_extraction/alzheimers_studies.json'
    ]
    
    # Load Phase 2 data files
    phase2_files = [
        'data/norml_extraction/chronic_pain_studies.json',
        'data/norml_extraction/anxiety_studies.json',
        'data/norml_extraction/depression_studies.json',
        'data/norml_extraction/arthritis_studies.json'
    ]
    
    # Load Phase 3 data files (expansion)
    phase3_files = [
        'data/norml_extraction/multiple_sclerosis_studies.json',
        'data/norml_extraction/nausea_chemotherapy_studies.json',
        'data/norml_extraction/ibd_crohns_studies.json',
        'data/norml_extraction/parkinsons_studies.json',
        'data/norml_extraction/glaucoma_studies.json',
        'data/norml_extraction/cancer_palliative_studies.json',
        'data/norml_extraction/appetite_cachexia_studies.json',
        'data/norml_extraction/tourette_syndrome_studies.json'
    ]
    
    all_studies = []
    conditions_loaded = []
    
    # Load Phase 1 legacy file if exists (skip if individual files present)
    if os.path.exists(phase1_legacy_file) and not os.path.exists(phase1_files[0]):
        with open(phase1_legacy_file, encoding='utf-8') as f:
            phase1 = json.load(f)
            all_studies.extend(phase1.get('studies', []))
            print(f"✓ Phase 1 (legacy): {len(phase1.get('studies', []))} studies loaded")
    
    # Load Phase 1 individual condition files
    for file_path in phase1_files:
        if os.path.exists(file_path):
            with open(file_path, encoding='utf-8') as f:
                data = json.load(f)
                all_studies.extend(data.get('studies', []))
                conditions_loaded.append(data['condition'])
                print(f"✓ {data['condition']}: {len(data.get('studies', []))} studies loaded")
    
    # Load Phase 2 files
    for file_path in phase2_files:
        if os.path.exists(file_path):
            with open(file_path, encoding='utf-8') as f:
                data = json.load(f)
                all_studies.extend(data.get('studies', []))
                conditions_loaded.append(data['condition'])
                print(f"✓ {data['condition']}: {len(data.get('studies', []))} studies loaded")
    
    # Load Phase 3 files
    for file_path in phase3_files:
        if os.path.exists(file_path):
            with open(file_path, encoding='utf-8') as f:
                data = json.load(f)
                all_studies.extend(data.get('studies', []))
                conditions_loaded.append(data['condition'])
                print(f"✓ {data['condition']}: {len(data.get('studies', []))} studies loaded")
    
    # Create merged dataset
    merged_data = {
        "extraction_date": "2025-12-19",
        "phases_complete": ["Phase 1", "Phase 2", "Phase 3"],
        "total_conditions": len(conditions_loaded),
        "total_studies": len(all_studies),
        "conditions": conditions_loaded,
        "studies": all_studies
    }
    
    # Save merged dataset
    with open(output_filename, 'w', encoding='utf-8') as f:
        json.dump(merged_data, f, indent=2)
    
    print(f"\n{'='*60}")
    print(f"✓ MERGED DATASET CREATED: {output_filename}")
    print(f"Total studies: {merged_data['total_studies']}")
    print(f"Conditions: {len(merged_data['conditions'])}")
    print(f"{'='*60}\n")

if __name__ == "__main__":
    print("\n🔬 NORML Study Extraction Validator")
    print("=" * 60)
    
    # Validate chronic pain extraction
    chronic_pain_file = 'data/norml_extraction/chronic_pain_studies.json'
    if os.path.exists(chronic_pain_file):
        validate_study_extraction(chronic_pain_file)
    else:
        print(f"⚠️  File not found: {chronic_pain_file}")
    
    # Try to merge if other files exist
    print("\n📊 Attempting to merge all available data...")
    merge_phase_data()
