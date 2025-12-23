from dataclasses import dataclass
from typing import Optional

@dataclass
class DimerEvidence:
    compound_1: str
    compound_2: str
    dimer_type: str  # e.g., 'cannabinoid-cannabinoid', 'cannabinoid-terpene'
    interaction_effect: str  # 'synergistic', 'antagonistic', 'additive', 'unknown'
    evidence_type: str  # 'clinical', 'preclinical', 'in silico', 'prediction', 'omics', 'pathway', 'ai-prediction'
    effect_size: Optional[str] = None
    study_type: Optional[str] = None
    source: Optional[str] = None
    confidence: Optional[str] = None
    notes: Optional[str] = None
    # Next-gen/omics/pathway/AI fields
    omics_signature: Optional[str] = None  # e.g., transcriptomic/proteomic/metabolomic evidence
    pathway: Optional[str] = None  # e.g., KEGG/Reactome/MetaCyc pathway name or ID
    synergy_score: Optional[float] = None  # Quantitative synergy metric (e.g., Bliss, Loewe)
    ai_predicted: Optional[bool] = None  # True if entry is AI-predicted
    ai_features: Optional[dict] = None  # Dict of AI-predicted features (e.g., {'therapeutic_potential_score': 0.85})
    patient_stratification: Optional[str] = None  # e.g., 'genotype', 'phenotype', 'age', 'sex', etc.
