"""
Triangulation Scoring Framework
Multi-source validation for dimeric cannabinoid predictions

Supports NeuroBotanica patent claims:
- Cross-validation of computational predictions
- Confidence bounds calculation
- Evidence triangulation from multiple sources
- Uncertainty quantification

Reference: NeuroBotanica MVP Development Plan - Week 4 Task 4.3
"""
import logging
from typing import Optional, List, Dict, Any, Tuple
from dataclasses import dataclass, field
from datetime import datetime
from enum import Enum
import math

import numpy as np

logger = logging.getLogger(__name__)


class ValidationSource(str, Enum):
    """Sources for prediction validation."""
    COMPUTATIONAL = "computational"  # Formation probability from algorithms
    STRUCTURAL = "structural"  # Similarity to known dimers
    ML_THERAPEUTIC = "ml_therapeutic"  # ML-predicted therapeutic potential
    EXPERIMENTAL = "experimental"  # Experimental validation
    LITERATURE = "literature"  # Literature evidence
    DATABASE = "database"  # Database cross-reference


class ExperimentalStatus(str, Enum):
    """Status of experimental validation."""
    CONFIRMED = "confirmed"  # Experimentally confirmed
    TENTATIVE = "tentative"  # Preliminary evidence
    PREDICTED = "predicted"  # No experimental data
    CONTRADICTED = "contradicted"  # Experimental evidence against


@dataclass
class ValidationEvidence:
    """Evidence from a single validation source."""
    source: ValidationSource
    score: float  # 0.0-1.0
    confidence: float  # Confidence in this source
    weight: float  # Relative weight for triangulation
    metadata: Dict[str, Any] = field(default_factory=dict)
    
    def weighted_contribution(self) -> float:
        """Calculate weighted contribution to triangulation."""
        return self.score * self.confidence * self.weight


@dataclass
class TriangulationResult:
    """Result of triangulation scoring."""
    # Core scores
    triangulation_score: float  # 0.0-1.0 composite score
    uncertainty: float  # Standard deviation of inputs
    
    # Confidence interval
    confidence_interval_lower: float
    confidence_interval_upper: float
    confidence_level: float  # e.g., 0.95 for 95%
    
    # Source breakdown
    num_sources: int
    sources_used: List[str]
    source_contributions: Dict[str, float]
    
    # Validation status
    experimental_status: ExperimentalStatus
    validation_grade: str  # A, B, C, D
    
    # Agreement metrics
    source_agreement: float  # 0.0-1.0 how well sources agree
    consensus_strength: str  # "strong", "moderate", "weak", "conflicting"
    
    # Recommendations
    confidence_level_label: str  # "high", "medium", "low"
    recommended_actions: List[str]
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for storage/API."""
        return {
            "triangulation_score": self.triangulation_score,
            "uncertainty": self.uncertainty,
            "confidence_interval": {
                "lower": self.confidence_interval_lower,
                "upper": self.confidence_interval_upper,
                "level": self.confidence_level
            },
            "num_sources": self.num_sources,
            "sources_used": self.sources_used,
            "source_contributions": self.source_contributions,
            "experimental_status": self.experimental_status.value,
            "validation_grade": self.validation_grade,
            "source_agreement": self.source_agreement,
            "consensus_strength": self.consensus_strength,
            "confidence_level_label": self.confidence_level_label,
            "recommended_actions": self.recommended_actions
        }


class TriangulationScorer:
    """Triangulate dimeric predictions using multiple validation methods.
    
    Implements three primary validation sources:
    1. Computational prediction (formation probability)
    2. Structural similarity to known dimers
    3. ML-predicted therapeutic potential
    
    Plus optional additional sources:
    4. Experimental validation
    5. Literature evidence
    6. Database cross-reference
    
    The triangulation score represents the confidence that a predicted
    dimer will have the expected properties, based on agreement across
    multiple independent sources.
    """
    
    # Default source weights
    DEFAULT_WEIGHTS = {
        ValidationSource.COMPUTATIONAL: 0.35,
        ValidationSource.STRUCTURAL: 0.25,
        ValidationSource.ML_THERAPEUTIC: 0.25,
        ValidationSource.EXPERIMENTAL: 0.50,  # High weight when available
        ValidationSource.LITERATURE: 0.30,
        ValidationSource.DATABASE: 0.20,
    }
    
    # Grade thresholds
    GRADE_THRESHOLDS = {
        "A": 0.80,  # Strong evidence, high confidence
        "B": 0.65,  # Good evidence, moderate confidence
        "C": 0.45,  # Limited evidence
        "D": 0.0,   # Insufficient evidence
    }
    
    # Confidence level thresholds
    CONFIDENCE_THRESHOLDS = {
        "high": 0.70,
        "medium": 0.45,
        "low": 0.0
    }
    
    def __init__(
        self,
        custom_weights: Optional[Dict[ValidationSource, float]] = None,
        confidence_level: float = 0.95,
        require_minimum_sources: int = 2
    ):
        """Initialize triangulation scorer.
        
        Args:
            custom_weights: Custom weights for validation sources
            confidence_level: Confidence level for interval (default 0.95)
            require_minimum_sources: Minimum sources for valid triangulation
        """
        self.weights = custom_weights or self.DEFAULT_WEIGHTS.copy()
        self.confidence_level = confidence_level
        self.require_minimum_sources = require_minimum_sources
        
        # Z-score for confidence interval
        self.z_score = self._calculate_z_score(confidence_level)
        
        logger.info(f"TriangulationScorer initialized with {confidence_level*100}% confidence")
    
    def _calculate_z_score(self, confidence_level: float) -> float:
        """Calculate Z-score for confidence interval."""
        # Common values
        z_scores = {
            0.90: 1.645,
            0.95: 1.960,
            0.99: 2.576
        }
        return z_scores.get(confidence_level, 1.960)
    
    def calculate_triangulation_score(
        self,
        formation_probability: Optional[float] = None,
        structural_similarity: Optional[float] = None,
        ml_therapeutic_score: Optional[float] = None,
        experimental_validation: Optional[str] = None,
        literature_score: Optional[float] = None,
        database_score: Optional[float] = None,
        additional_metadata: Optional[Dict[str, Any]] = None
    ) -> TriangulationResult:
        """Calculate triangulated confidence score.
        
        Args:
            formation_probability: 0.0-1.0 from computational model
            structural_similarity: 0.0-1.0 Tanimoto to known dimers
            ml_therapeutic_score: 0.0-1.0 from ML model
            experimental_validation: "confirmed", "tentative", "predicted", "contradicted"
            literature_score: 0.0-1.0 from literature evidence
            database_score: 0.0-1.0 from database cross-reference
            additional_metadata: Additional information for tracking
            
        Returns:
            TriangulationResult with score, confidence bounds, and recommendations
        """
        # Collect evidence from all available sources
        evidence = []
        
        # Source 1: Computational (formation probability)
        if formation_probability is not None:
            evidence.append(ValidationEvidence(
                source=ValidationSource.COMPUTATIONAL,
                score=formation_probability,
                confidence=0.85,  # High confidence in computational methods
                weight=self.weights[ValidationSource.COMPUTATIONAL],
                metadata={"method": "formation_probability"}
            ))
        
        # Source 2: Structural similarity
        if structural_similarity is not None:
            evidence.append(ValidationEvidence(
                source=ValidationSource.STRUCTURAL,
                score=structural_similarity,
                confidence=0.90,  # Very reliable metric
                weight=self.weights[ValidationSource.STRUCTURAL],
                metadata={"method": "tanimoto_similarity"}
            ))
        
        # Source 3: ML therapeutic prediction
        if ml_therapeutic_score is not None:
            evidence.append(ValidationEvidence(
                source=ValidationSource.ML_THERAPEUTIC,
                score=ml_therapeutic_score,
                confidence=0.70,  # ML predictions have inherent uncertainty
                weight=self.weights[ValidationSource.ML_THERAPEUTIC],
                metadata={"method": "ml_prediction"}
            ))
        
        # Source 4: Literature evidence
        if literature_score is not None:
            evidence.append(ValidationEvidence(
                source=ValidationSource.LITERATURE,
                score=literature_score,
                confidence=0.95,  # Literature is reliable
                weight=self.weights[ValidationSource.LITERATURE],
                metadata={"method": "literature_review"}
            ))
        
        # Source 5: Database cross-reference
        if database_score is not None:
            evidence.append(ValidationEvidence(
                source=ValidationSource.DATABASE,
                score=database_score,
                confidence=0.90,
                weight=self.weights[ValidationSource.DATABASE],
                metadata={"method": "database_lookup"}
            ))
        
        # Handle insufficient evidence
        if len(evidence) < self.require_minimum_sources:
            return self._insufficient_evidence_result(evidence, experimental_validation)
        
        # Calculate triangulation score
        triangulation_score, uncertainty = self._calculate_weighted_score(evidence)
        
        # Adjust for experimental validation
        exp_status = self._parse_experimental_status(experimental_validation)
        triangulation_score, uncertainty = self._adjust_for_experimental(
            triangulation_score, uncertainty, exp_status
        )
        
        # Calculate confidence interval
        ci_lower, ci_upper = self._calculate_confidence_interval(
            triangulation_score, uncertainty, len(evidence)
        )
        
        # Calculate source agreement
        source_agreement = self._calculate_source_agreement(evidence)
        
        # Determine consensus strength
        consensus = self._determine_consensus_strength(source_agreement, len(evidence))
        
        # Calculate validation grade
        grade = self._calculate_grade(triangulation_score, len(evidence), exp_status)
        
        # Determine confidence level label
        confidence_label = self._get_confidence_label(triangulation_score)
        
        # Generate recommendations
        recommendations = self._generate_recommendations(
            triangulation_score, grade, exp_status, source_agreement
        )
        
        # Build source contributions
        contributions = {
            e.source.value: round(e.weighted_contribution(), 3)
            for e in evidence
        }
        
        return TriangulationResult(
            triangulation_score=round(triangulation_score, 4),
            uncertainty=round(uncertainty, 4),
            confidence_interval_lower=round(ci_lower, 4),
            confidence_interval_upper=round(ci_upper, 4),
            confidence_level=self.confidence_level,
            num_sources=len(evidence),
            sources_used=[e.source.value for e in evidence],
            source_contributions=contributions,
            experimental_status=exp_status,
            validation_grade=grade,
            source_agreement=round(source_agreement, 3),
            consensus_strength=consensus,
            confidence_level_label=confidence_label,
            recommended_actions=recommendations
        )
    
    def _calculate_weighted_score(
        self,
        evidence: List[ValidationEvidence]
    ) -> Tuple[float, float]:
        """Calculate weighted average score and uncertainty."""
        if not evidence:
            return 0.0, 1.0
        
        # Normalize weights
        total_weight = sum(e.weight * e.confidence for e in evidence)
        
        if total_weight == 0:
            return 0.5, 0.5
        
        # Calculate weighted average
        weighted_sum = sum(
            e.score * e.weight * e.confidence
            for e in evidence
        )
        score = weighted_sum / total_weight
        
        # Calculate uncertainty (weighted standard deviation)
        scores = np.array([e.score for e in evidence])
        weights = np.array([e.weight * e.confidence for e in evidence])
        weights = weights / weights.sum()
        
        if len(scores) > 1:
            # Weighted variance
            mean = np.average(scores, weights=weights)
            variance = np.average((scores - mean) ** 2, weights=weights)
            uncertainty = math.sqrt(variance)
        else:
            uncertainty = 0.2  # Default uncertainty for single source
        
        return score, uncertainty
    
    def _adjust_for_experimental(
        self,
        score: float,
        uncertainty: float,
        exp_status: ExperimentalStatus
    ) -> Tuple[float, float]:
        """Adjust score and uncertainty based on experimental validation."""
        if exp_status == ExperimentalStatus.CONFIRMED:
            # Strong experimental support
            score = min(1.0, score * 1.15)
            uncertainty *= 0.5
        elif exp_status == ExperimentalStatus.TENTATIVE:
            # Some experimental support
            score = min(1.0, score * 1.05)
            uncertainty *= 0.8
        elif exp_status == ExperimentalStatus.CONTRADICTED:
            # Experimental evidence against prediction
            score *= 0.5
            uncertainty *= 1.5
        # PREDICTED status doesn't modify
        
        return score, uncertainty
    
    def _calculate_confidence_interval(
        self,
        score: float,
        uncertainty: float,
        n_sources: int
    ) -> Tuple[float, float]:
        """Calculate confidence interval using t-distribution approach."""
        # Adjust for sample size
        se = uncertainty / math.sqrt(max(1, n_sources))
        
        ci_lower = max(0.0, score - self.z_score * se)
        ci_upper = min(1.0, score + self.z_score * se)
        
        return ci_lower, ci_upper
    
    def _calculate_source_agreement(
        self,
        evidence: List[ValidationEvidence]
    ) -> float:
        """Calculate how well sources agree with each other."""
        if len(evidence) < 2:
            return 1.0
        
        scores = [e.score for e in evidence]
        
        # Agreement based on coefficient of variation
        mean = np.mean(scores)
        if mean == 0:
            return 0.5
        
        cv = np.std(scores) / mean
        
        # Convert CV to agreement score (1.0 = perfect agreement)
        agreement = max(0, 1 - cv)
        
        return agreement
    
    def _determine_consensus_strength(
        self,
        agreement: float,
        n_sources: int
    ) -> str:
        """Determine consensus strength from agreement score."""
        if agreement >= 0.85 and n_sources >= 3:
            return "strong"
        elif agreement >= 0.65 and n_sources >= 2:
            return "moderate"
        elif agreement >= 0.45:
            return "weak"
        else:
            return "conflicting"
    
    def _calculate_grade(
        self,
        score: float,
        n_sources: int,
        exp_status: ExperimentalStatus
    ) -> str:
        """Calculate validation grade (A-D)."""
        # Adjust threshold based on experimental status
        bonus = 0.0
        if exp_status == ExperimentalStatus.CONFIRMED:
            bonus = 0.1
        elif exp_status == ExperimentalStatus.CONTRADICTED:
            bonus = -0.2
        
        adjusted_score = score + bonus
        
        # Also consider number of sources
        if n_sources < 2:
            # Downgrade if insufficient sources
            adjusted_score -= 0.15
        
        if adjusted_score >= self.GRADE_THRESHOLDS["A"]:
            return "A"
        elif adjusted_score >= self.GRADE_THRESHOLDS["B"]:
            return "B"
        elif adjusted_score >= self.GRADE_THRESHOLDS["C"]:
            return "C"
        else:
            return "D"
    
    def _get_confidence_label(self, score: float) -> str:
        """Get confidence level label from score."""
        if score >= self.CONFIDENCE_THRESHOLDS["high"]:
            return "high"
        elif score >= self.CONFIDENCE_THRESHOLDS["medium"]:
            return "medium"
        else:
            return "low"
    
    def _parse_experimental_status(
        self,
        status_str: Optional[str]
    ) -> ExperimentalStatus:
        """Parse experimental status string."""
        if status_str is None:
            return ExperimentalStatus.PREDICTED
        
        status_map = {
            "confirmed": ExperimentalStatus.CONFIRMED,
            "tentative": ExperimentalStatus.TENTATIVE,
            "predicted": ExperimentalStatus.PREDICTED,
            "contradicted": ExperimentalStatus.CONTRADICTED,
            "none": ExperimentalStatus.PREDICTED,
        }
        
        return status_map.get(status_str.lower(), ExperimentalStatus.PREDICTED)
    
    def _generate_recommendations(
        self,
        score: float,
        grade: str,
        exp_status: ExperimentalStatus,
        agreement: float
    ) -> List[str]:
        """Generate actionable recommendations."""
        recommendations = []
        
        # High confidence - ready for next steps
        if grade == "A":
            recommendations.append(
                "High confidence prediction - suitable for experimental validation"
            )
            if exp_status != ExperimentalStatus.CONFIRMED:
                recommendations.append(
                    "Consider prioritizing for synthesis and testing"
                )
        
        # Good confidence - promising
        elif grade == "B":
            recommendations.append(
                "Promising prediction - additional validation recommended"
            )
            if agreement < 0.7:
                recommendations.append(
                    "Sources show some disagreement - investigate discrepancies"
                )
        
        # Limited confidence - needs more data
        elif grade == "C":
            recommendations.append(
                "Limited evidence - gather additional validation sources"
            )
            recommendations.append(
                "Consider computational refinement before experimental work"
            )
        
        # Insufficient confidence
        else:
            recommendations.append(
                "Insufficient evidence for reliable prediction"
            )
            recommendations.append(
                "Improve computational models or gather literature evidence"
            )
        
        # Experimental status recommendations
        if exp_status == ExperimentalStatus.TENTATIVE:
            recommendations.append(
                "Preliminary experimental evidence exists - consider follow-up studies"
            )
        elif exp_status == ExperimentalStatus.CONTRADICTED:
            recommendations.append(
                "WARNING: Experimental evidence contradicts prediction - re-evaluate"
            )
        
        return recommendations
    
    def _insufficient_evidence_result(
        self,
        evidence: List[ValidationEvidence],
        experimental: Optional[str]
    ) -> TriangulationResult:
        """Return result for insufficient evidence."""
        return TriangulationResult(
            triangulation_score=0.0,
            uncertainty=1.0,
            confidence_interval_lower=0.0,
            confidence_interval_upper=1.0,
            confidence_level=self.confidence_level,
            num_sources=len(evidence),
            sources_used=[e.source.value for e in evidence],
            source_contributions={},
            experimental_status=self._parse_experimental_status(experimental),
            validation_grade="D",
            source_agreement=0.0,
            consensus_strength="insufficient",
            confidence_level_label="low",
            recommended_actions=[
                f"Insufficient evidence - need at least {self.require_minimum_sources} validation sources",
                "Gather additional computational, structural, or literature evidence"
            ]
        )
    
    def batch_triangulate(
        self,
        predictions: List[Dict[str, Any]]
    ) -> List[TriangulationResult]:
        """Triangulate multiple predictions in batch.
        
        Args:
            predictions: List of dicts with prediction scores
            
        Returns:
            List of TriangulationResult objects
        """
        results = []
        
        for pred in predictions:
            result = self.calculate_triangulation_score(
                formation_probability=pred.get("formation_probability"),
                structural_similarity=pred.get("structural_similarity"),
                ml_therapeutic_score=pred.get("ml_therapeutic_score"),
                experimental_validation=pred.get("experimental_validation"),
                literature_score=pred.get("literature_score"),
                database_score=pred.get("database_score")
            )
            results.append(result)
        
        return results
    
    def rank_by_confidence(
        self,
        results: List[Tuple[str, TriangulationResult]]
    ) -> List[Tuple[str, TriangulationResult]]:
        """Rank predictions by triangulation confidence.
        
        Args:
            results: List of (identifier, result) tuples
            
        Returns:
            Sorted list (highest confidence first)
        """
        return sorted(
            results,
            key=lambda x: (
                x[1].triangulation_score,
                -x[1].uncertainty,
                x[1].num_sources
            ),
            reverse=True
        )
