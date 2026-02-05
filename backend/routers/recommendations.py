"""
ClinPath Recommendation Ranking Router
File: backend/routers/recommendations.py

Ranks available cultivars against a patient profile using ClinPath
therapeutic predictions. Returns top-N recommendations sorted by
predicted efficacy for the primary symptom.

🔌 PLUG-IN POINTS marked where trained ML models replace rule-based MVP logic.
"""

from fastapi import APIRouter, HTTPException
from pydantic import BaseModel, Field
from typing import Optional
from datetime import datetime
import uuid
import json

# Import prediction engine from clinpath router
from backend.routers.clinpath import (
    PatientProfile,
    FormulationInput,
    _calculate_therapeutic_scores,
    _apply_demographic_correction,
    _hash_patient_profile,
)

router = APIRouter(prefix="/api/clinpath", tags=["recommendations"])


# ============================================================================
# REQUEST/RESPONSE MODELS
# ============================================================================

class CultivarProfile(BaseModel):
    """A single cultivar available at the dispensary."""
    cultivar_id: str
    cultivar_name: str
    cannabinoid_profile: dict  # {"THC": 18.5, "CBD": 0.5, "CBG": 1.2, ...}
    terpene_profile: dict      # {"myrcene": 1.5, "limonene": 0.3, ...}
    polysaccharide_profile: Optional[dict] = None
    compound_ids: Optional[list[str]] = None


class RecommendationRequest(BaseModel):
    patient: PatientProfile
    primary_symptom: str          # 'anxiety', 'pain', 'insomnia', etc.
    secondary_symptoms: Optional[list[str]] = None
    contraindications: Optional[list[str]] = None
    available_cultivars: list[CultivarProfile]
    dispensary_id: Optional[str] = None
    budtender_id: Optional[str] = None
    top_n: int = 3                # Number of recommendations to return


class RankedRecommendation(BaseModel):
    rank: int
    cultivar_id: str
    cultivar_name: str
    predicted_efficacy: float     # 0-1 for primary symptom
    confidence_score: float
    evidence_quality: str

    # Therapeutic highlights for budtender talking points
    primary_score: float          # Score for the requested symptom
    side_effect_summary: dict     # Key risks to mention
    cognitive_summary: dict       # Memory/focus impacts
    talking_points: list[str]     # Plain-language points for the budtender

    # Full prediction reference
    prediction_id: str
    formulation_id: str
    cannabinoid_profile: dict
    terpene_profile: dict


class RecommendationResponse(BaseModel):
    recommendation_id: str
    primary_symptom: str
    patient_profile_hash: str
    recommendations: list[RankedRecommendation]
    total_cultivars_evaluated: int
    timestamp: str


# ============================================================================
# SYMPTOM-TO-SCORE MAPPING
# ============================================================================

SYMPTOM_SCORE_MAP = {
    "anxiety": "anxiolytic_score",
    "depression": "antidepressant_score",
    "insomnia": "sedative_score",
    "sleep": "sedative_score",
    "pain": "analgesic_score",
    "chronic_pain": "analgesic_score",
    "inflammation": "analgesic_score",
    "nausea": "analgesic_score",      # MVP approximation
    "appetite_loss": "sedative_score", # MVP approximation (THC-driven)
    "ptsd": "anxiolytic_score",        # MVP approximation
}


# ============================================================================
# CONTRAINDICATION CHECKS
# ============================================================================

def _check_contraindications(
    cultivar: CultivarProfile,
    contraindications: list[str],
) -> tuple[bool, str]:
    """
    Check if a cultivar is safe given patient contraindications.
    Returns (is_safe, reason).
    -------------------------------------------------------
    🔌 PLUG-IN POINT: Replace with DDI prediction engine
       from patent Claims 23-32.
    """
    cp = cultivar.cannabinoid_profile
    thc = cp.get("THC", 0)
    cbd = cp.get("CBD", 0)

    for contra in contraindications:
        contra_lower = contra.lower()

        if contra_lower == "pregnancy" and thc > 0:
            return False, "THC contraindicated during pregnancy"

        if contra_lower == "psychosis" and thc > 5:
            return False, f"THC {thc}% exceeds safe threshold for psychosis history"

        if contra_lower == "heart_condition" and thc > 10:
            return False, f"THC {thc}% may increase cardiovascular risk"

        if contra_lower in ("anticoagulant", "blood_thinner") and cbd > 20:
            return False, f"CBD {cbd}% may interact with anticoagulants (CYP2C9)"

        if contra_lower == "liver_disease" and (thc + cbd) > 15:
            return False, "High cannabinoid load contraindicated with liver disease"

    return True, ""


# ============================================================================
# TALKING POINT GENERATOR
# ============================================================================

def _generate_talking_points(
    cultivar: CultivarProfile,
    scores: dict,
    primary_symptom: str,
) -> list[str]:
    """
    Generate plain-language talking points for budtender education.
    """
    points = []
    cp = cultivar.cannabinoid_profile
    tp = cultivar.terpene_profile

    # Primary efficacy
    score_key = SYMPTOM_SCORE_MAP.get(primary_symptom, "analgesic_score")
    efficacy = scores.get(score_key, 0)
    if efficacy >= 0.7:
        points.append(f"Strong match for {primary_symptom} relief ({efficacy:.0%} predicted efficacy)")
    elif efficacy >= 0.4:
        points.append(f"Moderate match for {primary_symptom} relief ({efficacy:.0%} predicted efficacy)")
    else:
        points.append(f"Mild support for {primary_symptom} ({efficacy:.0%} predicted efficacy)")

    # THC/CBD ratio context
    thc = cp.get("THC", 0)
    cbd = cp.get("CBD", 0)
    if thc > 0 and cbd > 0:
        ratio = thc / cbd
        if ratio > 5:
            points.append(f"THC-dominant ({thc}% THC / {cbd}% CBD) — expect psychoactive effects")
        elif ratio > 1:
            points.append(f"Balanced-leaning THC ({thc}% THC / {cbd}% CBD) — moderate psychoactivity")
        else:
            points.append(f"CBD-dominant ({cbd}% CBD / {thc}% THC) — minimal psychoactivity")
    elif cbd > 0:
        points.append(f"CBD-only ({cbd}% CBD) — non-intoxicating")

    # Key terpene highlights
    if tp.get("myrcene", 0) > 0.5:
        points.append("High myrcene — may enhance sedation and body relaxation")
    if tp.get("limonene", 0) > 0.3:
        points.append("Notable limonene — associated with mood elevation")
    if tp.get("linalool", 0) > 0.2:
        points.append("Contains linalool — calming, may reduce anxiety")
    if tp.get("beta_caryophyllene", 0) > 0.2:
        points.append("Beta-caryophyllene present — anti-inflammatory via CB2 activation")
    if tp.get("pinene", 0) > 0.2:
        points.append("Pinene detected — may counteract THC memory impairment")

    # Side effect warning
    if scores.get("psychoactivity_risk", 0) > 0.6:
        points.append("⚠️ High psychoactivity risk — recommend starting with low dose")
    if scores.get("sedation_risk", 0) > 0.5:
        points.append("⚠️ Sedation likely — best for evening/nighttime use")

    return points


# ============================================================================
# RANKING ENGINE
# ============================================================================

def _rank_cultivars(
    request: RecommendationRequest,
) -> list[dict]:
    """
    Score and rank all available cultivars for the patient.
    -------------------------------------------------------
    🔌 PLUG-IN POINT: Replace with ML-based ranking model
       trained on recommendation acceptance data.
    """
    score_key = SYMPTOM_SCORE_MAP.get(
        request.primary_symptom, "analgesic_score"
    )
    results = []

    for cultivar in request.available_cultivars:
        # Check contraindications first
        if request.contraindications:
            is_safe, reason = _check_contraindications(
                cultivar, request.contraindications
            )
            if not is_safe:
                continue  # Skip unsafe cultivars

        # Build formulation input
        formulation = FormulationInput(
            compound_ids=cultivar.compound_ids or [],
            cannabinoid_profile=cultivar.cannabinoid_profile,
            terpene_profile=cultivar.terpene_profile,
            polysaccharide_profile=cultivar.polysaccharide_profile,
            primary_indication=request.primary_symptom,
            secondary_indications=request.secondary_symptoms,
        )

        # Calculate scores
        scores = _calculate_therapeutic_scores(formulation, request.patient)

        # Apply demographic correction
        scores, adjustment_factor = _apply_demographic_correction(
            scores, request.patient
        )

        # Extract primary efficacy
        primary_score = scores.get(score_key, 0)

        # Generate talking points
        talking_points = _generate_talking_points(
            cultivar, scores, request.primary_symptom
        )

        results.append({
            "cultivar": cultivar,
            "scores": scores,
            "primary_score": primary_score,
            "talking_points": talking_points,
            "adjustment_factor": adjustment_factor,
        })

    # Sort by primary symptom score (descending)
    results.sort(key=lambda x: x["primary_score"], reverse=True)

    return results


# ============================================================================
# ENDPOINTS
# ============================================================================

@router.post("/recommend", response_model=RecommendationResponse)
async def get_recommendations(request: RecommendationRequest):
    """
    Rank available cultivars and return top-N recommendations
    with budtender talking points.
    """
    recommendation_id = f"rec_{uuid.uuid4().hex[:16]}"
    profile_hash = _hash_patient_profile(request.patient)

    # Rank all cultivars
    ranked = _rank_cultivars(request)

    # Build top-N response
    recommendations = []
    for i, result in enumerate(ranked[:request.top_n]):
        cultivar = result["cultivar"]
        scores = result["scores"]

        rec = RankedRecommendation(
            rank=i + 1,
            cultivar_id=cultivar.cultivar_id,
            cultivar_name=cultivar.cultivar_name,
            predicted_efficacy=round(result["primary_score"], 4),
            confidence_score=scores["overall_confidence"],
            evidence_quality=scores["evidence_quality"],
            primary_score=round(result["primary_score"], 4),
            side_effect_summary={
                "psychoactivity_risk": scores["psychoactivity_risk"],
                "dependence_risk": scores["dependence_risk"],
                "sedation_risk": scores["sedation_risk"],
                "anxiety_risk": scores["anxiety_risk"],
            },
            cognitive_summary={
                "memory_impact": scores["memory_impact"],
                "focus_impact": scores["focus_impact"],
                "neuroprotection": scores["neuroprotection_score"],
            },
            talking_points=result["talking_points"],
            prediction_id=f"cp_{uuid.uuid4().hex[:16]}",
            formulation_id=f"form_{uuid.uuid4().hex[:16]}",
            cannabinoid_profile=cultivar.cannabinoid_profile,
            terpene_profile=cultivar.terpene_profile,
        )
        recommendations.append(rec)

    return RecommendationResponse(
        recommendation_id=recommendation_id,
        primary_symptom=request.primary_symptom,
        patient_profile_hash=profile_hash,
        recommendations=recommendations,
        total_cultivars_evaluated=len(request.available_cultivars),
        timestamp=datetime.utcnow().isoformat(),
    )


@router.get("/symptoms", response_model=list[str])
async def list_supported_symptoms():
    """Return list of symptoms the recommendation engine supports."""
    return list(SYMPTOM_SCORE_MAP.keys())
