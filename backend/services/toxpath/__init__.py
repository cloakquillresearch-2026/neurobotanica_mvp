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
    "MemoConfig"
]
