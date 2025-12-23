import random
import json
from pathlib import Path

def generate_patient_profiles(n=100):
    age_ranges = ["18-25", "26-35", "36-45", "46-55", "56-65", "66+"]
    exp_levels = ["naive", "occasional", "regular", "daily"]
    tol_levels = ["low", "medium", "high"]
    qualities = ["high", "good", "moderate", "low"]
    profiles = []
    for _ in range(n):
        profile = {
            "demographics": {
                "age_range": random.choice(age_ranges),
                "experience_level": random.choice(exp_levels),
                "thc_tolerance": random.choice(tol_levels)
            },
            "treatment": {
                "duration_days": random.randint(7, 180)
            },
            "response": {
                "adherence_score": random.randint(5, 10),
                "adverse_effects": random.sample(["dry mouth", "drowsiness", "anxiety", "euphoria", "none"], k=random.randint(0,2)),
                "overall_efficacy_score": random.randint(4, 10)
            },
            "data_quality": {
                "source_study_quality": random.choice(qualities)
            }
        }
        profiles.append(profile)
    return profiles

if __name__ == "__main__":
    project_root = Path(__file__).parent.parent
    patient_data = {
        "metadata": {"total_patients": 100},
        "profiles": generate_patient_profiles(100)
    }
    out_path = project_root / "data" / "training" / "neurobotanica_enriched_patients.json"
    with open(out_path, "w") as f:
        json.dump({"patient_response_data": patient_data}, f, indent=2)
    print(f"Generated 100 patient profiles at {out_path}")
