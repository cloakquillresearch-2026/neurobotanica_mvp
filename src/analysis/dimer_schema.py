from dataclasses import dataclass
from typing import Optional

@dataclass
class DimerEvidence:
    compound_1: str
    compound_2: str
    dimer_type: str  # e.g., 'cannabinoid-cannabinoid', 'cannabinoid-terpene'
    interaction_effect: str  # 'synergistic', 'antagonistic', 'additive', 'unknown'
    evidence_type: str  # 'clinical', 'preclinical', 'in silico', 'prediction'
    effect_size: Optional[str] = None
    study_type: Optional[str] = None
    source: Optional[str] = None
    confidence: Optional[str] = None
    notes: Optional[str] = None
