"""
RegPath - Regulatory Pathway Optimization Engine

Trade Secret Module for NeuroBotanica.
Produces manufacturer-facing Regulatory Strategy Memos with:
- Pathway recommendations (IND/NDA/505(b)(2)/ANDA)
- Readiness checklists
- Timeline estimates

Exposes capabilities via API while keeping decision matrices/weights/heuristics proprietary.
"""

from .strategist import (
    RegPathStrategist,
    RegPathRequest,
    RegPathResponse,
    ProductProfile,
    ProductType,
    NoveltyLevel,
    RegulatoryPathway,
    PathwayRecommendation,
    ChecklistItem,
    ChecklistCategory,
    TimelineMilestone,
    MilestonePhase,
    EvidenceInputs,
    StrategyConstraints
)

from .memo_generator import (
    RegPathMemoGenerator,
    MemoFormat,
    MemoConfig
)

__all__ = [
    # Main strategist
    "RegPathStrategist",
    
    # Input models
    "RegPathRequest",
    "ProductProfile",
    "ProductType",
    "NoveltyLevel",
    "EvidenceInputs",
    "StrategyConstraints",
    
    # Output models
    "RegPathResponse",
    "RegulatoryPathway",
    "PathwayRecommendation",
    "ChecklistItem",
    "ChecklistCategory",
    "TimelineMilestone",
    "MilestonePhase",
    
    # Memo generation
    "RegPathMemoGenerator",
    "MemoFormat",
    "MemoConfig"
]
