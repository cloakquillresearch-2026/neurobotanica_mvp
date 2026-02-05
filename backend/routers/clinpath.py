"""
ClinPath API Router — Therapeutic Prediction Endpoints
File: backend/routers/clinpath.py

Integrates with neurobotanica_formulations, clinpath_predictions,


🔌 PLUG-IN POINTS marked where trained ML models replace rule-based MVP logic.
"""

from fastapi import APIRouter, HTTPException
from pydantic import BaseModel, Field
from typing import Optional
from datetime import datetime
import hashlib
import json
import uuid

router = APIRouter(prefix="/api/clinpath", tags=["clinpath"])


# ============================================================================
# REQUEST/RESPONSE MODELS
# ============================================================================

class PatientProfile(BaseModel):
    age_range: Optional[str] = None
    sex: Optional[str] = None
    genetic_ancestry: Optional[str] = "unknown"
    cyp2c9_variant: Optional[str] = None
    faah_variant: Optional[str] = None


class FormulationInput(BaseModel):
    compound_ids: list[str]
    cannabinoid_profile: dict
    terpene_profile: dict
    polysaccharide_profile: Optional[dict] = None
    primary_indication: str
    secondary_indications: Optional[list[str]] = None


class PredictionRequest(BaseModel):
    patient: PatientProfile
    formulation: FormulationInput
    dispensary_id: Optional[str] = None
    customer_id: Optional[str] = None


class TherapeuticPrediction(BaseModel):
    prediction_id: str
    formulation_id: str
    anxiolytic_score: float
    antidepressant_score: float
    sedative_score: float
    analgesic_score: float
    memory_impact: float
    focus_impact: float
    neuroprotection_score: float
    psychoactivity_risk: float
    dependence_risk: float
    sedation_risk: float
    anxiety_risk: float
    overall_confidence: float
    evidence_quality: str
    bias_correction_applied: bool = True
    demographic_adjustment_factor: float = 1.0
    model_version: str = "v1.0-mvp"


# ============================================================================
# PREDICTION ENGINE (MVP — rule-based, placeholder for trained ML models)
# ============================================================================

def _hash_patient_profile(patient: PatientProfile) -> str:
    """SHA-256 hash of demographics for privacy."""
    profile_str = json.dumps({
        "age_range": patient.age_range,
        "sex": patient.sex,
        "ancestry": patient.genetic_ancestry,
    }, sort_keys=True)
    return hashlib.sha256(profile_str.encode()).hexdigest()


def _calculate_therapeutic_scores(
    formulation: FormulationInput,
    patient: PatientProfile,
) -> dict:
    """
    MVP prediction engine using compound-effect mappings.
    -------------------------------------------------------
    🔌 PLUG-IN POINT: Replace this function with trained ML model
       to achieve 88-92% accuracy target from patent spec.
    """
    cp = formulation.cannabinoid_profile
    tp = formulation.terpene_profile

    thc = cp.get("THC", 0)
    cbd = cp.get("CBD", 0)
    cbg = cp.get("CBG", 0)
    cbn = cp.get("CBN", 0)
    myrcene = tp.get("myrcene", 0)
    linalool = tp.get("linalool", 0)
    limonene = tp.get("limonene", 0)
    caryophyllene = tp.get("beta_caryophyllene", 0)
    pinene = tp.get("pinene", 0)

    # --- Therapeutic scores (0-1) ---
    anxiolytic = min(1.0, (cbd * 0.03) + (linalool * 0.15) + (cbg * 0.02))
    antidepressant = min(1.0, (cbd * 0.02) + (limonene * 0.18) + (cbg * 0.03))
    sedative = min(1.0, (thc * 0.015) + (myrcene * 0.2) + (cbn * 0.05) + (linalool * 0.1))
    analgesic = min(1.0, (thc * 0.02) + (cbd * 0.025) + (caryophyllene * 0.15))

    # --- Cognitive effects (-1 to +1) ---
    memory_impact = max(-1.0, min(1.0, -0.02 * thc + 0.01 * cbd + 0.05 * pinene))
    focus_impact = max(-1.0, min(1.0, -0.015 * thc + 0.02 * cbd + 0.04 * pinene + 0.03 * limonene))
    neuroprotection = min(1.0, cbd * 0.03 + cbg * 0.04 + caryophyllene * 0.05)

    # --- Side effects (0-1 probability) ---
    psychoactivity_risk = min(1.0, thc * 0.04)
    dependence_risk = min(1.0, thc * 0.008)
    sedation_risk = min(1.0, (thc * 0.01) + (myrcene * 0.15) + (cbn * 0.04))
    anxiety_risk = min(1.0, max(0, thc * 0.025 - cbd * 0.015))

    # --- Confidence based on evidence coverage ---
    known_compounds = sum(1 for v in [thc, cbd, cbg, cbn] if v > 0)
    known_terpenes = sum(1 for v in [myrcene, linalool, limonene, caryophyllene, pinene] if v > 0)
    confidence = min(1.0, 0.4 + (known_compounds * 0.08) + (known_terpenes * 0.06))

    evidence = "high" if confidence >= 0.75 else ("medium" if confidence >= 0.55 else "low")

    return {
        "anxiolytic_score": round(anxiolytic, 4),
        "antidepressant_score": round(antidepressant, 4),
        "sedative_score": round(sedative, 4),
        "analgesic_score": round(analgesic, 4),
        "memory_impact": round(memory_impact, 4),
        "focus_impact": round(focus_impact, 4),
        "neuroprotection_score": round(neuroprotection, 4),
        "psychoactivity_risk": round(psychoactivity_risk, 4),
        "dependence_risk": round(dependence_risk, 4),
        "sedation_risk": round(sedation_risk, 4),
        "anxiety_risk": round(anxiety_risk, 4),
        "overall_confidence": round(confidence, 4),
        "evidence_quality": evidence,
    }


def _apply_demographic_correction(
    scores: dict,
    patient: PatientProfile,
) -> tuple[dict, float]:
    """
    Apply demographic bias correction per patent Section 3.8.5.
    -------------------------------------------------------
    🔌 PLUG-IN POINT: Replace with BioPath 96% accuracy model.
    """
    factor = 1.0

    # CYP2C9 poor metabolizer — higher cannabinoid exposure
    if patient.cyp2c9_variant in ("*1/*3", "*3/*3"):
        factor *= 1.3
        scores["psychoactivity_risk"] = min(1.0, scores["psychoactivity_risk"] * 1.3)
        scores["sedation_risk"] = min(1.0, scores["sedation_risk"] * 1.2)

    # FAAH A/A variant — elevated endocannabinoid tone
    if patient.faah_variant == "A/A":
        scores["anxiolytic_score"] = min(1.0, scores["anxiolytic_score"] * 1.15)
        factor *= 0.95

    # Age correction (elderly)
    if patient.age_range == "65+":
        scores["sedation_risk"] = min(1.0, scores["sedation_risk"] * 1.4)
        scores["psychoactivity_risk"] = min(1.0, scores["psychoactivity_risk"] * 1.25)
        factor *= 1.2

    return scores, round(factor, 4)


# ============================================================================
# ENDPOINTS
# ============================================================================

@router.post("/predict", response_model=TherapeuticPrediction)
async def predict_therapeutic_profile(request: PredictionRequest):
    """
    Generate ClinPath therapeutic prediction for a formulation + patient.
    This is the primary endpoint the Budtender UI calls.
    """
    prediction_id = f"cp_{uuid.uuid4().hex[:16]}"
    formulation_id = f"form_{uuid.uuid4().hex[:16]}"
    profile_hash = _hash_patient_profile(request.patient)

    # Step 1: Calculate base therapeutic scores
    scores = _calculate_therapeutic_scores(
        request.formulation, request.patient
    )

    # Step 2: Apply demographic bias correction
    scores, adjustment_factor = _apply_demographic_correction(
        scores, request.patient
    )

    # Step 3: Build response
    return TherapeuticPrediction(
        prediction_id=prediction_id,
        formulation_id=formulation_id,
        **scores,
        bias_correction_applied=True,
        demographic_adjustment_factor=adjustment_factor,
    )
            indication=request.indication,
