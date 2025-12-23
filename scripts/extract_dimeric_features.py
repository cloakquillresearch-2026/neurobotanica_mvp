import json
from src.analysis.dimer_features import dimer_to_features
from src.analysis.dimer_prediction_extraction import extract_dimeric_predictions

if __name__ == "__main__":
    # Example: Convert dimeric predictions to features for model input
    dimers = extract_dimeric_predictions("data/training/neurobotanica_dimeric_predictions.json")
    features = [dimer_to_features(d) for d in dimers]
    with open("data/training/dimeric_features.json", "w") as f:
        json.dump(features, f, indent=2)
    print(f"Extracted features for {len(features)} dimeric entries.")
