import json
from src.analysis.nextgen_dimer_ingest import load_nextgen_dimer_evidence

if __name__ == "__main__":
    # Load and validate next-gen dimer/entourage evidence
    dimers = load_nextgen_dimer_evidence("data/processed/nextgen_dimers.json")
    print(f"Loaded {len(dimers)} next-gen dimer/entourage entries.")
    # Placeholder: pass to model pipeline or save for further processing
    # (Model pipeline implementation to follow)
