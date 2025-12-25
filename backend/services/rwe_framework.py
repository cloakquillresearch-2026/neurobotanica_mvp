"""
Real-World Evidence (RWE) Framework

Collects and analyzes real-world evidence from Nevada dispensaries
to validate therapeutic predictions and improve model accuracy.

Critical for:
- BioPath validation (community evidence)
- ClinPath optimization (outcome data)
- TherapeuticPredictionModel training data

Reference: Nevada Pilot Launch - March 2026
"""

from typing import Dict, List, Optional, Any, Tuple
from dataclasses import dataclass
from datetime import datetime, timedelta
from enum import Enum
import logging
import hashlib

logger = logging.getLogger(__name__)


class EvidenceType(str, Enum):
    """Types of real-world evidence collected."""
    SALE = "sale"  # Product sale record
    PATIENT_FEEDBACK = "patient_feedback"  # Patient-reported outcome
    ADVERSE_EVENT = "adverse_event"  # Reported side effect
    OUTCOME_REPORT = "outcome_report"  # Structured outcome data
    HEALER_VALIDATION = "healer_validation"  # Community healer assessment


class OutcomeCategory(str, Enum):
    """Patient outcome categories."""
    SIGNIFICANT_IMPROVEMENT = "significant_improvement"
    MODERATE_IMPROVEMENT = "moderate_improvement"
    SLIGHT_IMPROVEMENT = "slight_improvement"
    NO_CHANGE = "no_change"
    SLIGHT_WORSENING = "slight_worsening"
    SIGNIFICANT_WORSENING = "significant_worsening"


@dataclass
class PatientOutcome:
    """Anonymized patient outcome data."""
    patient_hash: str
    condition: str
    cannabinoid: str
    dose_mg: float
    frequency: str
    duration_days: int
    outcome: OutcomeCategory
    outcome_score: float  # 0-10
    side_effects: List[str]
    confidence: float  # Data quality score


@dataclass
class RWEAggregation:
    """Aggregated RWE statistics for a compound-condition pair."""
    compound: str
    condition: str
    total_patients: int
    improvement_rate: float
    average_outcome_score: float
    confidence_interval: Tuple[float, float]
    common_side_effects: List[Dict[str, Any]]
    recommended_dose_range: Tuple[float, float]
    evidence_quality: float


class RealWorldEvidenceFramework:
    """Framework for collecting and analyzing real-world evidence.
    
    Features:
    - Patient privacy protection (hashing, anonymization)
    - Outcome aggregation and analysis
    - Integration with BioPath validation
    - Dispensary performance tracking
    """
    
    # Minimum samples for reliable aggregation
    MIN_SAMPLES_FOR_AGGREGATION = 30
    
    # Outcome score thresholds
    IMPROVEMENT_THRESHOLD = 6.0  # Score >= 6 is improvement
    SIGNIFICANT_IMPROVEMENT_THRESHOLD = 8.0
    
    def __init__(self):
        """Initialize RWE framework."""
        self.evidence_buffer: List[PatientOutcome] = []
        self.aggregation_cache: Dict[str, RWEAggregation] = {}
    
    def anonymize_patient_id(self, patient_id: str, dispensary_id: str) -> str:
        """Create anonymous patient hash for tracking.
        
        Uses SHA-256 with dispensary salt to prevent cross-dispensary tracking
        while allowing within-dispensary outcome correlation.
        """
        salt = f"neurobotanica_rwe_{dispensary_id}"
        combined = f"{patient_id}:{salt}"
        return hashlib.sha256(combined.encode()).hexdigest()[:32]
    
    def record_patient_outcome(
        self,
        patient_id: str,
        dispensary_id: str,
        condition: str,
        cannabinoid: str,
        dose_mg: float,
        frequency: str,
        duration_days: int,
        outcome_score: float,
        side_effects: Optional[List[str]] = None,
        consent_verified: bool = True
    ) -> Optional[PatientOutcome]:
        """Record a patient outcome with privacy protection.
        
        Args:
            patient_id: Raw patient ID (will be hashed)
            dispensary_id: Source dispensary
            condition: Target condition being treated
            cannabinoid: Primary cannabinoid used
            dose_mg: Daily dose in mg
            frequency: Dosing frequency
            duration_days: Duration of treatment
            outcome_score: Patient-reported outcome (0-10)
            side_effects: List of reported side effects
            consent_verified: Whether patient consented to data use
            
        Returns:
            PatientOutcome if recorded, None if consent not verified
        """
        if not consent_verified:
            logger.warning("Cannot record outcome without consent verification")
            return None
        
        # Anonymize patient
        patient_hash = self.anonymize_patient_id(patient_id, dispensary_id)
        
        # Categorize outcome
        if outcome_score >= self.SIGNIFICANT_IMPROVEMENT_THRESHOLD:
            outcome_category = OutcomeCategory.SIGNIFICANT_IMPROVEMENT
        elif outcome_score >= self.IMPROVEMENT_THRESHOLD:
            outcome_category = OutcomeCategory.MODERATE_IMPROVEMENT
        elif outcome_score >= 5.0:
            outcome_category = OutcomeCategory.SLIGHT_IMPROVEMENT
        elif outcome_score >= 4.0:
            outcome_category = OutcomeCategory.NO_CHANGE
        elif outcome_score >= 2.0:
            outcome_category = OutcomeCategory.SLIGHT_WORSENING
        else:
            outcome_category = OutcomeCategory.SIGNIFICANT_WORSENING
        
        # Calculate confidence based on data completeness
        confidence = self._calculate_confidence(
            duration_days=duration_days,
            has_side_effects=side_effects is not None,
            outcome_extremity=abs(outcome_score - 5.0)
        )
        
        outcome = PatientOutcome(
            patient_hash=patient_hash,
            condition=condition.lower(),
            cannabinoid=cannabinoid.upper(),
            dose_mg=dose_mg,
            frequency=frequency,
            duration_days=duration_days,
            outcome=outcome_category,
            outcome_score=outcome_score,
            side_effects=side_effects or [],
            confidence=confidence
        )
        
        self.evidence_buffer.append(outcome)
        logger.info(f"Recorded RWE: {cannabinoid} for {condition}, score={outcome_score}")
        
        return outcome
    
    def _calculate_confidence(
        self,
        duration_days: int,
        has_side_effects: bool,
        outcome_extremity: float
    ) -> float:
        """Calculate confidence score for evidence entry.
        
        Higher confidence for:
        - Longer treatment duration (more reliable)
        - Complete side effect reporting
        - Moderate outcomes (extreme outcomes may be biased)
        """
        base_confidence = 0.5
        
        # Duration bonus (up to 0.2)
        if duration_days >= 30:
            base_confidence += 0.2
        elif duration_days >= 14:
            base_confidence += 0.15
        elif duration_days >= 7:
            base_confidence += 0.1
        
        # Side effect reporting bonus
        if has_side_effects:
            base_confidence += 0.1
        
        # Moderate outcome bonus (extreme outcomes get slight penalty)
        if outcome_extremity <= 2.0:
            base_confidence += 0.1
        elif outcome_extremity <= 3.5:
            base_confidence += 0.05
        
        return min(1.0, base_confidence)
    
    def aggregate_evidence(
        self,
        compound: str,
        condition: str,
        min_samples: Optional[int] = None
    ) -> Optional[RWEAggregation]:
        """Aggregate RWE for a compound-condition pair.
        
        Args:
            compound: Cannabinoid compound
            condition: Target condition
            min_samples: Minimum samples required (default: 30)
            
        Returns:
            RWEAggregation if sufficient data, None otherwise
        """
        min_samples = min_samples or self.MIN_SAMPLES_FOR_AGGREGATION
        
        # Filter relevant outcomes
        relevant = [
            o for o in self.evidence_buffer
            if o.cannabinoid == compound.upper() and o.condition == condition.lower()
        ]
        
        if len(relevant) < min_samples:
            logger.info(
                f"Insufficient samples for {compound}/{condition}: "
                f"{len(relevant)}/{min_samples}"
            )
            return None
        
        # Calculate statistics
        scores = [o.outcome_score for o in relevant]
        confidences = [o.confidence for o in relevant]
        doses = [o.dose_mg for o in relevant]
        
        avg_score = sum(scores) / len(scores)
        improvement_count = sum(1 for o in relevant if o.outcome_score >= self.IMPROVEMENT_THRESHOLD)
        improvement_rate = improvement_count / len(relevant)
        
        # Confidence interval (simple 95% CI approximation)
        import statistics
        if len(scores) > 1:
            std_dev = statistics.stdev(scores)
            margin = 1.96 * (std_dev / (len(scores) ** 0.5))
            ci = (max(0, avg_score - margin), min(10, avg_score + margin))
        else:
            ci = (avg_score, avg_score)
        
        # Aggregate side effects
        side_effect_counts: Dict[str, int] = {}
        for o in relevant:
            for se in o.side_effects:
                side_effect_counts[se] = side_effect_counts.get(se, 0) + 1
        
        common_side_effects = [
            {"effect": se, "frequency": count / len(relevant)}
            for se, count in sorted(side_effect_counts.items(), key=lambda x: -x[1])[:5]
        ]
        
        # Dose range (10th to 90th percentile)
        sorted_doses = sorted(doses)
        dose_10 = sorted_doses[int(len(sorted_doses) * 0.1)]
        dose_90 = sorted_doses[int(len(sorted_doses) * 0.9)]
        
        # Overall evidence quality
        evidence_quality = sum(confidences) / len(confidences)
        
        aggregation = RWEAggregation(
            compound=compound,
            condition=condition,
            total_patients=len(relevant),
            improvement_rate=improvement_rate,
            average_outcome_score=avg_score,
            confidence_interval=ci,
            common_side_effects=common_side_effects,
            recommended_dose_range=(dose_10, dose_90),
            evidence_quality=evidence_quality
        )
        
        # Cache result
        cache_key = f"{compound}:{condition}"
        self.aggregation_cache[cache_key] = aggregation
        
        return aggregation
    
    def get_biopath_evidence(
        self,
        compound: str,
        condition: str
    ) -> Dict[str, Any]:
        """Format RWE for BioPath validation input.
        
        Returns evidence in format expected by BiasAwareValidationEngine.
        """
        aggregation = self.aggregate_evidence(compound, condition, min_samples=10)
        
        if aggregation is None:
            return {
                "source_type": "real_world_evidence",
                "score": 0.0,
                "sample_size": 0,
                "available": False
            }
        
        return {
            "source_type": "real_world_evidence",
            "score": aggregation.improvement_rate,
            "sample_size": aggregation.total_patients,
            "study_quality": aggregation.evidence_quality,
            "confidence_interval": aggregation.confidence_interval,
            "available": True
        }
    
    def get_clinpath_outcomes(
        self,
        compound: str,
        condition: str
    ) -> Dict[str, Any]:
        """Format RWE for ClinPath trial optimization.
        
        Returns outcome data for approval probability prediction.
        """
        aggregation = self.aggregate_evidence(compound, condition, min_samples=10)
        
        if aggregation is None:
            return {
                "rwe_available": False,
                "improvement_rate": 0.0,
                "sample_size": 0
            }
        
        return {
            "rwe_available": True,
            "improvement_rate": aggregation.improvement_rate,
            "average_outcome": aggregation.average_outcome_score,
            "sample_size": aggregation.total_patients,
            "recommended_dose_range": aggregation.recommended_dose_range,
            "safety_profile": {
                "common_side_effects": aggregation.common_side_effects,
                "evidence_quality": aggregation.evidence_quality
            }
        }
    
    def get_dispensary_performance(self, dispensary_id: str) -> Dict[str, Any]:
        """Get performance metrics for a dispensary.
        
        Used for pilot program monitoring.
        """
        # This would query the database in production
        return {
            "dispensary_id": dispensary_id,
            "total_rwe_entries": 0,
            "unique_patients": 0,
            "conditions_covered": [],
            "average_confidence": 0.0,
            "data_quality_score": 0.0,
            "last_submission": None
        }


# Global instance
_rwe_framework: Optional[RealWorldEvidenceFramework] = None


def get_rwe_framework() -> RealWorldEvidenceFramework:
    """Get global RWE framework instance."""
    global _rwe_framework
    if _rwe_framework is None:
        _rwe_framework = RealWorldEvidenceFramework()
    return _rwe_framework
