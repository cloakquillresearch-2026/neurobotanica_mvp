import json
from src.analysis.dimer_schema import DimerEvidence
from scripts.extract_and_validate_dimers import extract_dimeric_predictions, extract_synergy_evidence

def save_dimers(dimers, out_path):
    # Convert dataclass objects to dicts for JSON serialization
    with open(out_path, 'w') as f:
        json.dump([d.__dict__ for d in dimers], f, indent=2)

def main():
    ai_dimers = extract_dimeric_predictions("data/training/neurobotanica_dimeric_predictions.json")
    synergy_dimers = extract_synergy_evidence("data/training/neurobotanica_enriched_synergies.json")
    all_dimers = ai_dimers + synergy_dimers
    save_dimers(all_dimers, "data/processed/validated_dimers_nextgen.json")
    print(f"Integrated and saved {len(all_dimers)} validated dimer/entourage entries.")

if __name__ == "__main__":
    main()
