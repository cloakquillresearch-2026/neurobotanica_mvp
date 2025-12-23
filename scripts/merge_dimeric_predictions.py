import json
from src.analysis.dimer_prediction_extraction import extract_dimeric_predictions
from src.analysis.dimer_integration import merge_dimeric_evidence

# Example usage: merge dimeric predictions into existing training data
if __name__ == "__main__":
    # Load dimeric predictions
    dimeric_dimers = extract_dimeric_predictions(
        "data/training/neurobotanica_dimeric_predictions.json"
    )
    # Load existing training data
    with open("data/training/neurobotanica_enriched_dataset.json", "r") as f:
        training_data = json.load(f)
    # Merge
    merged = merge_dimeric_evidence(training_data, [d.__dict__ for d in dimeric_dimers])
    # Save merged dataset
    with open("data/training/neurobotanica_enriched_with_dimers.json", "w") as f:
        json.dump(merged, f, indent=2)
    print(f"Merged {len(dimeric_dimers)} dimeric entries. Total: {len(merged)}.")
