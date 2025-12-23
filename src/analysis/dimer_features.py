from .dimer_schema import DimerEvidence

def dimer_to_features(dimer: DimerEvidence) -> dict:
    features = {
        'compound_1': dimer.compound_1,
        'compound_2': dimer.compound_2,
        'dimer_type': dimer.dimer_type,
        'interaction_effect': dimer.interaction_effect,
        # Add more feature engineering as needed
    }
    return features
