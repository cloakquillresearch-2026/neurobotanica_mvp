from dataclasses import dataclass, field
from typing import Optional, List, Dict

@dataclass
class NextGenDimerEvidence:
    compound_1: str
    compound_2: str
    dimer_type: str  # e.g., 'cannabinoid-cannabinoid', 'cannabinoid-terpene', 'cannabinoid-flavonoid', etc.
    interaction_effect: str  # 'synergistic', 'antagonistic', 'additive', 'unknown'
    evidence_type: str  # 'clinical', 'preclinical', 'omics', 'pathway', 'ai-prediction', 'real-world', 'patient-stratified', etc.
    effect_size: Optional[str] = None
    study_type: Optional[str] = None
    source: Optional[str] = None
    confidence: Optional[str] = None
    notes: Optional[str] = None
    omics_signature: Optional[str] = None
    pathway: Optional[str] = None
    synergy_score: Optional[float] = None
    ai_predicted: Optional[bool] = None
    ai_features: Optional[Dict] = field(default_factory=dict)
    patient_stratification: Optional[str] = None
    population_group: Optional[str] = None  # e.g., 'elderly', 'female', 'genotype: CYP2C9*3', etc.
    time_series: Optional[List[Dict]] = field(default_factory=list)  # For longitudinal/real-world evidence
    adverse_effects: Optional[List[str]] = field(default_factory=list)
    benefit_risk_ratio: Optional[float] = None
    regulatory_status: Optional[str] = None  # e.g., 'FDA-approved', 'experimental', etc.
    external_links: Optional[List[str]] = field(default_factory=list)  # PubMed, clinicaltrials.gov, etc.
