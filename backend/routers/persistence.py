"""
ClinPath D1 Persistence Layer
File: backend/routers/persistence.py

Writes ClinPath predictions and budtender recommendations
back to D1 tables for audit trail, analytics, and the
budtender_recommendation_details view.

Endpoints:
  POST /api/clinpath/persist/prediction   — Save a prediction to D1
  POST /api/clinpath/persist/recommendation — Save a recommendation to D1
  GET  /api/clinpath/history/{dispensary_id} — Recommendation history
  GET  /api/clinpath/analytics/top-cultivars — Top-performing cultivars
"""

from fastapi import APIRouter, HTTPException, Depends
from pydantic import BaseModel
from typing import Optional
from datetime import datetime
import json
import uuid
import os

router = APIRouter(prefix="/api/clinpath/persist", tags=["persistence"])


# ============================================================================
# D1 DATABASE CONNECTION
# ============================================================================

# When running via Cloudflare Worker proxy, D1 is accessed through
# the worker binding. For local dev, use wrangler d1 local mode.
# This module expects a D1-compatible execute function to be injected.

class D1Connection:
    """
    Abstraction over D1 database access.
    In production: called via Worker binding.
    In dev: called via wrangler d1 execute --local.
    """
    def __init__(self, execute_fn=None):
        self._execute = execute_fn

    async def execute(self, sql: str, params: list = None) -> dict:
        if self._execute:
            return await self._execute(sql, params or [])
        raise HTTPException(
            status_code=503,
            detail="D1 database connection not configured"
        )


# Global connection — injected at startup or by worker proxy
db = D1Connection()


def set_d1_connection(execute_fn):
    """Called by main.py or worker proxy to inject D1 access."""
    global db
    db = D1Connection(execute_fn)


# ============================================================================
# REQUEST MODELS
# ============================================================================

class PersistPredictionRequest(BaseModel):
    prediction_id: str
    formulation_id: str
    patient_profile_hash: str

    # Therapeutic scores
    anxiolytic_score: float
    antidepressant_score: float
    sedative_score: float
    analgesic_score: float

    # Cognitive effects
    memory_impact: Optional[float] = None
    focus_impact: Optional[float] = None
    neuroprotection_score: Optional[float] = None

    # Side effects
    psychoactivity_risk: float
    dependence_risk: float
    sedation_risk: float
    anxiety_risk: float

    # Confidence
    overall_confidence: float
    evidence_quality: str
    model_version: str = "v1.0-mvp"


class PersistFormulationRequest(BaseModel):
    formulation_id: str
    customer_id: Optional[str] = None
    dispensary_id: Optional[str] = None
    compound_ids: list[str]
    cannabinoid_profile: dict
    terpene_profile: dict
    polysaccharide_profile: Optional[dict] = None
    predicted_efficacy: Optional[float] = None
    predicted_safety: Optional[float] = None
    synergy_score: Optional[float] = None
    confidence_score: Optional[float] = None
    primary_indication: Optional[str] = None
    secondary_indications: Optional[list[str]] = None
    tk_flag: bool = False
    consent_id: Optional[str] = None


class PersistRecommendationRequest(BaseModel):
    recommendation_id: str
    patient_age_range: Optional[str] = None
    patient_sex: Optional[str] = None
    primary_symptom: str
    secondary_symptoms: Optional[list[str]] = None
    contraindications: Optional[list[str]] = None
    available_cultivars: list[dict]
    recommended_cultivar_id: str
    recommended_cultivar_name: str
    predicted_efficacy: float
    confidence_score: float
    rank: int = 1
    formulation_id: Optional[str] = None
    clinpath_prediction_id: Optional[str] = None
    budtender_id: Optional[str] = None
    dispensary_id: Optional[str] = None
    accepted: Optional[bool] = None


class PersistDemographicRequest(BaseModel):
    factor_id: str
    prediction_id: str
    cyp2c9_variant: Optional[str] = None
    faah_variant: Optional[str] = None
    cnr1_variant: Optional[str] = None
    fucosidase_variant: Optional[str] = None
    dectin1_variant: Optional[str] = None
    genetic_ancestry: Optional[str] = "unknown"
    bias_correction_applied: bool = True
    demographic_adjustment_factor: float = 1.0
    equalized_odds_score: Optional[float] = None


# ============================================================================
# PERSIST ENDPOINTS
# ============================================================================

@router.post("/formulation")
async def persist_formulation(request: PersistFormulationRequest):
    """
    Save a formulation to the neurobotanica_formulations table.
    Must be called before persisting a prediction that references it.
    """
    sql = """
        INSERT INTO neurobotanica_formulations (
            formulation_id, customer_id, dispensary_id,
            compound_ids, cannabinoid_profile, terpene_profile,
            polysaccharide_profile, predicted_efficacy, predicted_safety,
            synergy_score, confidence_score, primary_indication,
            secondary_indications, tk_flag, consent_id
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
    """
    params = [
        request.formulation_id,
        request.customer_id,
        request.dispensary_id,
        json.dumps(request.compound_ids),
        json.dumps(request.cannabinoid_profile),
        json.dumps(request.terpene_profile),
        json.dumps(request.polysaccharide_profile) if request.polysaccharide_profile else None,
        request.predicted_efficacy,
        request.predicted_safety,
        request.synergy_score,
        request.confidence_score,
        request.primary_indication,
        json.dumps(request.secondary_indications) if request.secondary_indications else None,
        request.tk_flag,
        request.consent_id,
    ]

    try:
        await db.execute(sql, params)
        return {"status": "saved", "formulation_id": request.formulation_id}
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Failed to persist formulation: {str(e)}")


@router.post("/prediction")
async def persist_prediction(request: PersistPredictionRequest):
    """
    Save a ClinPath prediction to the clinpath_predictions table.
    The referenced formulation_id must already exist.
    """
    sql = """
        INSERT INTO clinpath_predictions (
            prediction_id, formulation_id, patient_profile_hash,
            anxiolytic_score, antidepressant_score, sedative_score, analgesic_score,
            memory_impact, focus_impact, neuroprotection_score,
            psychoactivity_risk, dependence_risk, sedation_risk, anxiety_risk,
            overall_confidence, evidence_quality, model_version
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
    """
    params = [
        request.prediction_id,
        request.formulation_id,
        request.patient_profile_hash,
        request.anxiolytic_score,
        request.antidepressant_score,
        request.sedative_score,
        request.analgesic_score,
        request.memory_impact,
        request.focus_impact,
        request.neuroprotection_score,
        request.psychoactivity_risk,
        request.dependence_risk,
        request.sedation_risk,
        request.anxiety_risk,
        request.overall_confidence,
        request.evidence_quality,
        request.model_version,
    ]

    try:
        await db.execute(sql, params)
        return {"status": "saved", "prediction_id": request.prediction_id}
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Failed to persist prediction: {str(e)}")


@router.post("/demographic")
async def persist_demographic_factors(request: PersistDemographicRequest):
    """
    Save demographic factors for bias correction tracking.
    """
    sql = """
        INSERT INTO clinpath_demographic_factors (
            factor_id, prediction_id,
            cyp2c9_variant, faah_variant, cnr1_variant,
            fucosidase_variant, dectin1_variant,
            genetic_ancestry, bias_correction_applied,
            demographic_adjustment_factor, equalized_odds_score
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
    """
    params = [
        request.factor_id,
        request.prediction_id,
        request.cyp2c9_variant,
        request.faah_variant,
        request.cnr1_variant,
        request.fucosidase_variant,
        request.dectin1_variant,
        request.genetic_ancestry,
        request.bias_correction_applied,
        request.demographic_adjustment_factor,
        request.equalized_odds_score,
    ]

    try:
        await db.execute(sql, params)
        return {"status": "saved", "factor_id": request.factor_id}
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Failed to persist demographic factors: {str(e)}")


@router.post("/recommendation")
async def persist_recommendation(request: PersistRecommendationRequest):
    """
    Save a budtender recommendation to the budtender_recommendations table.
    """
    sql = """
        INSERT INTO budtender_recommendations (
            recommendation_id, patient_age_range, patient_sex,
            primary_symptom, secondary_symptoms, contraindications,
            available_cultivars, recommended_cultivar_id, recommended_cultivar_name,
            predicted_efficacy, confidence_score, rank,
            formulation_id, clinpath_prediction_id,
            budtender_id, dispensary_id, accepted
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
    """
    params = [
        request.recommendation_id,
        request.patient_age_range,
        request.patient_sex,
        request.primary_symptom,
        json.dumps(request.secondary_symptoms) if request.secondary_symptoms else None,
        json.dumps(request.contraindications) if request.contraindications else None,
        json.dumps(request.available_cultivars),
        request.recommended_cultivar_id,
        request.recommended_cultivar_name,
        request.predicted_efficacy,
        request.confidence_score,
        request.rank,
        request.formulation_id,
        request.clinpath_prediction_id,
        request.budtender_id,
        request.dispensary_id,
        request.accepted,
    ]

    try:
        await db.execute(sql, params)
        return {"status": "saved", "recommendation_id": request.recommendation_id}
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Failed to persist recommendation: {str(e)}")


@router.put("/recommendation/{recommendation_id}/accept")
async def mark_recommendation_accepted(
    recommendation_id: str,
    accepted: bool = True,
):
    """
    Mark whether the customer accepted (purchased) the recommendation.
    Called by the budtender UI after consultation.
    """
    sql = """
        UPDATE budtender_recommendations
        SET accepted = ?
        WHERE recommendation_id = ?
    """
    try:
        await db.execute(sql, [accepted, recommendation_id])
        return {"status": "updated", "recommendation_id": recommendation_id, "accepted": accepted}
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Failed to update recommendation: {str(e)}")


# ============================================================================
# QUERY ENDPOINTS (Read from D1 views)
# ============================================================================

@router.get("/history/{dispensary_id}")
async def get_recommendation_history(
    dispensary_id: str,
    limit: int = 50,
    offset: int = 0,
):
    """
    Get recommendation history for a dispensary.
    Uses the budtender_recommendation_details view.
    """
    sql = """
        SELECT brd.*
        FROM budtender_recommendation_details brd
        JOIN budtender_recommendations br ON br.recommendation_id = brd.recommendation_id
        WHERE br.dispensary_id = ?
        ORDER BY brd.created_at DESC
        LIMIT ? OFFSET ?
    """
    try:
        result = await db.execute(sql, [dispensary_id, limit, offset])
        return {"dispensary_id": dispensary_id, "recommendations": result}
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Failed to retrieve history: {str(e)}")


@router.get("/analytics/top-cultivars")
async def get_top_cultivars(
    dispensary_id: Optional[str] = None,
    symptom: Optional[str] = None,
    limit: int = 10,
):
    """
    Get top-performing cultivars from the high_confidence_recommendations view.
    Optionally filter by dispensary or symptom.
    """
    conditions = []
    params = []

    base_sql = """
        SELECT
            recommended_cultivar_name,
            primary_symptom,
            COUNT(*) as times_recommended,
            AVG(predicted_efficacy) as avg_efficacy,
            AVG(confidence_score) as avg_confidence
        FROM high_confidence_recommendations hcr
    """

    if dispensary_id or symptom:
        # Join back to full table for filters
        base_sql = """
            SELECT
                hcr.recommended_cultivar_name,
                hcr.primary_symptom,
                COUNT(*) as times_recommended,
                AVG(hcr.predicted_efficacy) as avg_efficacy,
                AVG(hcr.confidence_score) as avg_confidence
            FROM high_confidence_recommendations hcr
            JOIN budtender_recommendations br ON br.recommendation_id = hcr.recommendation_id
        """
        if dispensary_id:
            conditions.append("br.dispensary_id = ?")
            params.append(dispensary_id)
        if symptom:
            conditions.append("hcr.primary_symptom = ?")
            params.append(symptom)

    if conditions:
        base_sql += " WHERE " + " AND ".join(conditions)

    base_sql += """
        GROUP BY recommended_cultivar_name, primary_symptom
        ORDER BY avg_efficacy DESC
        LIMIT ?
    """
    params.append(limit)

    try:
        result = await db.execute(base_sql, params)
        return {"top_cultivars": result}
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Failed to retrieve analytics: {str(e)}")
*** End Patch