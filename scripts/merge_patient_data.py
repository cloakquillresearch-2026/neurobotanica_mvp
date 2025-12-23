import json
from pathlib import Path

def merge_patient_data():
    """Merge patient_response_data into the fully enriched dataset."""
    project_root = Path(__file__).parent.parent
    enriched_path = project_root / 'data' / 'training' / 'neurobotanica_fully_enriched_fixed.json'
    patient_path = project_root / 'data' / 'training' / 'neurobotanica_enriched_patients.json'
    if not enriched_path.exists() or not patient_path.exists():
        print("Required files not found.")
        return
    with open(enriched_path, 'r', encoding='utf-8') as f:
        enriched = json.load(f)
    with open(patient_path, 'r', encoding='utf-8') as f:
        patient = json.load(f)
    # Merge patient_response_data at the top level
    enriched['patient_response_data'] = patient['patient_response_data']
    with open(enriched_path, 'w', encoding='utf-8') as f:
        json.dump(enriched, f, indent=2)
    print(f"Merged patient_response_data into {enriched_path}")

if __name__ == "__main__":
    merge_patient_data()
