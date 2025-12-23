import json
from pathlib import Path

def update_fully_enriched_dataset():
    """Update the fully enriched dataset to include merged dimeric data."""
    project_root = Path(__file__).parent.parent
    merged_path = project_root / 'data' / 'training' / 'neurobotanica_enriched_with_dimers.json'
    target_path = project_root / 'data' / 'training' / 'neurobotanica_fully_enriched_fixed.json'
    if merged_path.exists():
        with open(merged_path, 'r', encoding='utf-8') as f:
            merged = json.load(f)
        # Wrap in expected structure if needed
        if isinstance(merged, list):
            merged = {'compounds': merged}
        with open(target_path, 'w', encoding='utf-8') as f:
            json.dump(merged, f, indent=2)
        print(f"Updated {target_path} with merged dimeric data.")
    else:
        print(f"Merged dimeric file not found: {merged_path}")

if __name__ == "__main__":
    update_fully_enriched_dataset()
