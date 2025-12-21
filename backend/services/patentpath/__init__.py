"""
PatentPath Lite - NeuroBotanica Week 6-8
AI-powered patent analysis for cannabis innovation.

Features:
1. Prior Art Search - USPTO PatentsView & EPO Open Patent Services
2. Novelty Scoring - Similarity-based novelty assessment
3. FTO Checker - Freedom to operate analysis
4. Claim Generator - Templated patent claim drafting
5. TK Attribution Checker - Traditional knowledge compliance
6. Cost Estimator - USPTO fees and attorney cost projection

This module provides lightweight patent research tools
for early-stage innovation assessment.

DISCLAIMER: This is decision support only. Not legal advice.
Consult qualified patent counsel before filing or commercial decisions.
"""
from .prior_art import PriorArtSearcher, PriorArtResult, SearchQuery, PatentSource
from .novelty import NoveltyScorer, NoveltyReport, NoveltyLevel
from .fto_checker import FTOChecker, FTOReport, RiskLevel, FTO_DISCLAIMER
from .claim_generator import ClaimGenerator, ClaimGenerationReport, ClaimType, CLAIM_DISCLAIMER
from .tk_checker import TKAttributionChecker, TKCheckResult, TKIndicatorStrength, SacredKnowledgeStatus
from .cost_estimator import CostEstimator, CostEstimate, ProtectionStrategy, EntitySize, ComplexityLevel

__all__ = [
    # Prior Art
    "PriorArtSearcher",
    "PriorArtResult",
    "SearchQuery",
    "PatentSource",
    # Novelty
    "NoveltyScorer",
    "NoveltyReport",
    "NoveltyLevel",
    # FTO
    "FTOChecker",
    "FTOReport",
    "RiskLevel",
    "FTO_DISCLAIMER",
    # Claims
    "ClaimGenerator",
    "ClaimGenerationReport",
    "ClaimType",
    "CLAIM_DISCLAIMER",
    # TK Attribution (Week 8)
    "TKAttributionChecker",
    "TKCheckResult",
    "TKIndicatorStrength",
    "SacredKnowledgeStatus",
    # Cost Estimation (Week 8)
    "CostEstimator",
    "CostEstimate",
    "ProtectionStrategy",
    "EntitySize",
    "ComplexityLevel"
]
