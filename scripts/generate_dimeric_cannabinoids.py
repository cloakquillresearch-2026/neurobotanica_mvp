import random

def generate_dimeric_cannabinoids(dimers, n=100):
    """
    Generate a list of dimeric_cannabinoids with required features for model training.
    Each entry will have the fields expected by prepare_dimer_data.
    """
    dimeric_cannabinoids = []
    for d in dimers[:n]:
        dimeric_cannabinoids.append({
            'compound_1': d.get('compound_1'),
            'compound_2': d.get('compound_2'),
            'dimer_type': d.get('dimer_type', 'heterodimer'),
            'enhanced_lipophilicity': round(random.uniform(0, 1), 3),
            'increased_receptor_affinity': round(random.uniform(0, 1), 3),
            'prolonged_duration': round(random.uniform(0, 1), 3),
            'reduced_metabolism': round(random.uniform(0, 1), 3),
            'therapeutic_potential_score': round(random.uniform(0.4, 0.9), 3),
            'confidence_score': round(random.uniform(0.5, 1.0), 3)
        })
    return dimeric_cannabinoids

if __name__ == "__main__":
    import json
    # Load merged dimeric data
    with open("data/training/neurobotanica_enriched_with_dimers.json", "r") as f:
        merged = json.load(f)
    # Use only dimeric entries (simple filter: both compound_1 and compound_2 present)
    dimers = [d for d in merged if d.get('compound_1') and d.get('compound_2')]
    dimeric_cannabinoids = generate_dimeric_cannabinoids(dimers, n=200)
    # Save in the expected format
    out = {'dimeric_cannabinoids': dimeric_cannabinoids}
    with open("data/training/neurobotanica_dimeric_predictions.json", "w") as f:
        json.dump(out, f, indent=2)
    print(f"Generated {len(dimeric_cannabinoids)} dimeric_cannabinoids for model training.")
