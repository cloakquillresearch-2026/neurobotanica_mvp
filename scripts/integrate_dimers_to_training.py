import json
from pathlib import Path

# Integrate validated dimer/entourage evidence into the main training dataset
VALIDATED_DIMERS_PATH = Path('data/processed/validated_dimers_nextgen.json')
TRAINING_PATH = Path('data/training/neurobotanica_enriched_dataset.json')

# Load validated dimers
with open(VALIDATED_DIMERS_PATH, 'r', encoding='utf-8') as f:
    validated_dimers = json.load(f)

# Load existing training data
if TRAINING_PATH.exists():
    with open(TRAINING_PATH, 'r', encoding='utf-8') as f:
        training_data = json.load(f)
else:
    training_data = []

# Merge: add validated dimers as a new field if not present
if isinstance(training_data, dict):
    training_data['validated_dimers_nextgen'] = validated_dimers
else:
    # If list, wrap in dict
    training_data = {'compounds': training_data, 'validated_dimers_nextgen': validated_dimers}

with open(TRAINING_PATH, 'w', encoding='utf-8') as f:
    json.dump(training_data, f, indent=2)

print(f"Integrated {len(validated_dimers)} validated dimers into {TRAINING_PATH}")
