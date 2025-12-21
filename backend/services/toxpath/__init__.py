"""
ToxPath - Toxicity Risk Assessment Engine

Trade Secret Module for NeuroBotanica.
Exposes capabilities via API while keeping scoring logic/weights/heuristics proprietary.

Features:
- Structural alert detection
- Risk tier classification
- Route-aware risk framing
- Testing plan generation
- Toxicity memo generation
- Cannabis safety assessment (Schedule III compliance)
- Drug-drug interaction prediction
- Abuse liability assessment
"""

from .assessor import (
    ToxPathAssessor,
    ToxPathRequest,
    ToxPathResponse,
    RiskTier,
    RiskSummary,
    StructuralAlert,
    AlertSeverity,
    TestingStep,
    ExposureProfile,
    ImpurityEntry,
    AdministrationRoute
)

from .memo_generator import (
    ToxPathMemoGenerator,
    MemoFormat,
    MemoConfig
)

from .cannabis_safety import (
    CannabisSafetyAssessor,
    AbuseClass,
    RiskLevel,
    InteractionType,
    InteractionMechanism,
    VulnerablePopulation,
    AdverseEvent,
    DrugInteraction,
    AbuseLiabilityProfile,
    SafetyProfile,
    SafetyRecommendation,
)

__all__ = [
    # Main assessor
    "ToxPathAssessor",
    
    # Input models
    "ToxPathRequest",
    "ExposureProfile",
    "ImpurityEntry",
    "AdministrationRoute",
    
    # Output models
    "ToxPathResponse",
    "RiskTier",
    "RiskSummary",
    "StructuralAlert",
    "AlertSeverity",
    "TestingStep",
    
    # Memo generation
    "ToxPathMemoGenerator",
    "MemoFormat",
    "MemoConfig",
    
    # Cannabis safety (Schedule III)
    "CannabisSafetyAssessor",
    "AbuseClass",
    "RiskLevel",
    "InteractionType",
    "InteractionMechanism",
    "VulnerablePopulation",
    "AdverseEvent",
    "DrugInteraction",
    "AbuseLiabilityProfile",
    "SafetyProfile",
    "SafetyRecommendation",
]
