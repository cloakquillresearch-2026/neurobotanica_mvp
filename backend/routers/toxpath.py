"""
ToxPath API Router - Toxicity Risk Assessment Endpoints

Exposes ToxPath assessment capabilities while protecting proprietary algorithms.
"""

from fastapi import APIRouter, HTTPException
from pydantic import BaseModel, Field
from typing import Optional, List, Dict, Any
from datetime import datetime
import uuid

from ..services.toxpath import (
    ToxPathAssessor,
    ToxPathRequest as ToxPathRequestInternal,
    ExposureProfile,
    ImpurityEntry,
    ToxPathMemoGenerator,
    MemoFormat,
    MemoConfig
)


router = APIRouter(prefix="/api/v1/toxpath", tags=["ToxPath - Toxicity Assessment"])


# ============================================================================
# Pydantic Request/Response Models
# ============================================================================

class ExposureProfileRequest(BaseModel):
    """Exposure profile for risk assessment."""
    dose_mg: Optional[float] = Field(None, description="Dose in milligrams")
    frequency: Optional[str] = Field(None, description="Dosing frequency: once, daily, multiple_daily")
    duration: Optional[str] = Field(None, description="Duration: acute, subacute, chronic")


class ImpurityEntryRequest(BaseModel):
    """Known impurity entry."""
    name: str = Field(..., description="Impurity name")
    ppm: float = Field(..., description="Concentration in parts per million")
    cas_number: Optional[str] = Field(None, description="CAS registry number")


class ToxPathAssessRequest(BaseModel):
    """Request model for ToxPath assessment."""
    compound_ref: str = Field(..., description="SMILES string or compound identifier")
    compound_name: Optional[str] = Field(None, description="Human-readable compound name")
    route: Optional[str] = Field(None, description="Administration route: oral, inhaled, sublingual, topical, intravenous")
    exposure: Optional[ExposureProfileRequest] = Field(None, description="Expected exposure profile")
    known_impurities: Optional[List[ImpurityEntryRequest]] = Field(None, description="Known impurities")
    chempath_job_id: Optional[str] = Field(None, description="Link to ChemPath analysis job")
    
    class Config:
        json_schema_extra = {
            "example": {
                "compound_ref": "CCO",
                "compound_name": "Ethanol",
                "route": "oral",
                "exposure": {
                    "dose_mg": 100,
                    "frequency": "daily",
                    "duration": "chronic"
                }
            }
        }


class StructuralAlertResponse(BaseModel):
    """Structural alert in response."""
    code: str
    name: str
    severity: str
    description: str
    mechanism: Optional[str]
    affected_organs: List[str]


class TestingStepResponse(BaseModel):
    """Testing step recommendation."""
    order: int
    test_name: str
    test_type: str
    rationale: str
    estimated_cost_range: str
    timeline: str
    required: bool
    gmp_glp_required: bool


class RiskSummaryResponse(BaseModel):
    """Risk summary in response."""
    overall_tier: int
    tier_label: str
    testing_depth: str
    tier_rationale: str
    top_risks: List[str]
    key_unknowns: List[str]
    key_assumptions: List[str]
    route_specific_concerns: List[str]


class ToxPathAssessResponse(BaseModel):
    """Full ToxPath assessment response."""
    toxpath_assessment_id: str
    compound_ref: str
    compound_name: Optional[str]
    route: str
    risk_summary: RiskSummaryResponse
    alerts: List[StructuralAlertResponse]
    testing_plan: List[TestingStepResponse]
    consultation_required: bool
    consultation_flags: List[str]
    properties_snapshot: Dict[str, Any]
    status: str
    error: Optional[str]
    assessed_at: str


class MemoRequest(BaseModel):
    """Request for memo generation."""
    format: str = Field("markdown", description="Output format: markdown, html, json")
    include_testing_costs: bool = Field(True, description="Include cost estimates")
    include_timeline: bool = Field(True, description="Include timeline estimates")
    include_consultation_section: bool = Field(True, description="Include consultation section")
    include_assumptions: bool = Field(True, description="Include assumptions and limitations")
    company_name: Optional[str] = Field(None, description="Organization name for header")
    prepared_by: Optional[str] = Field(None, description="Preparer name for header")


class MemoResponse(BaseModel):
    """Response with generated memo."""
    assessment_id: str
    format: str
    memo_content: str
    generated_at: str


class RiskTierInfo(BaseModel):
    """Risk tier information."""
    tier: int
    label: str
    testing_depth: str


class AlertCodeInfo(BaseModel):
    """Alert code information."""
    code: str
    description: str


# ============================================================================
# In-memory storage (MVP - replace with database in production)
# ============================================================================

_assessments_store: Dict[str, ToxPathAssessResponse] = {}
_assessor = ToxPathAssessor()


# ============================================================================
# Endpoints
# ============================================================================

@router.post("/assess", response_model=ToxPathAssessResponse)
async def assess_compound(request: ToxPathAssessRequest):
    """Run ToxPath toxicity risk assessment on a compound.
    
    Analyzes structural features, detects toxicity alerts, calculates risk tier,
    and generates testing recommendations.
    
    Trade secret: Risk scoring weights, alert patterns, and tier thresholds
    are proprietary and not exposed in the response.
    """
    # Convert request to internal format
    exposure = None
    if request.exposure:
        exposure = ExposureProfile(
            dose_mg=request.exposure.dose_mg,
            frequency=request.exposure.frequency,
            duration=request.exposure.duration
        )
    
    impurities = None
    if request.known_impurities:
        impurities = [
            ImpurityEntry(
                name=imp.name,
                ppm=imp.ppm,
                cas_number=imp.cas_number
            )
            for imp in request.known_impurities
        ]
    
    internal_request = ToxPathRequestInternal(
        compound_ref=request.compound_ref,
        compound_name=request.compound_name,
        route=request.route,
        exposure=exposure,
        known_impurities=impurities,
        chempath_job_id=request.chempath_job_id
    )
    
    # Run assessment
    result = _assessor.assess(internal_request)
    
    # Convert to response format
    response_data = result.to_dict()
    response = ToxPathAssessResponse(**response_data)
    
    # Store for later retrieval
    _assessments_store[response.toxpath_assessment_id] = response
    
    return response


@router.get("/assessments/{assessment_id}", response_model=ToxPathAssessResponse)
async def get_assessment(assessment_id: str):
    """Retrieve a previous ToxPath assessment by ID."""
    if assessment_id not in _assessments_store:
        raise HTTPException(status_code=404, detail="Assessment not found")
    
    return _assessments_store[assessment_id]


@router.post("/assessments/{assessment_id}/memo", response_model=MemoResponse)
async def generate_memo(assessment_id: str, request: MemoRequest):
    """Generate a formatted toxicology memo for an assessment.
    
    Supports markdown, HTML, and JSON output formats.
    """
    if assessment_id not in _assessments_store:
        raise HTTPException(status_code=404, detail="Assessment not found")
    
    # Get stored assessment
    stored = _assessments_store[assessment_id]
    
    # Need to reconstruct internal response object for memo generator
    # This is a simplification - in production, store internal objects
    internal_result = _reconstruct_internal_response(stored)
    
    # Parse format
    try:
        format_enum = MemoFormat(request.format.lower())
    except ValueError:
        raise HTTPException(
            status_code=400, 
            detail=f"Invalid format. Supported: {[f.value for f in MemoFormat]}"
        )
    
    # Configure and generate
    config = MemoConfig(
        include_testing_costs=request.include_testing_costs,
        include_timeline=request.include_timeline,
        include_consultation_section=request.include_consultation_section,
        include_assumptions=request.include_assumptions,
        company_name=request.company_name,
        prepared_by=request.prepared_by
    )
    
    generator = ToxPathMemoGenerator(config)
    memo_content = generator.generate(internal_result, format_enum)
    
    return MemoResponse(
        assessment_id=assessment_id,
        format=format_enum.value,
        memo_content=memo_content,
        generated_at=datetime.utcnow().isoformat()
    )


@router.get("/risk-tiers", response_model=List[RiskTierInfo])
async def get_risk_tiers():
    """Get information about risk tier classifications."""
    return _assessor.get_risk_tiers()


@router.get("/alert-codes", response_model=List[AlertCodeInfo])
async def get_alert_codes():
    """Get available structural alert codes and descriptions."""
    codes = _assessor.get_alert_codes()
    return [AlertCodeInfo(code=k, description=v) for k, v in codes.items()]


@router.get("/routes", response_model=List[str])
async def get_routes():
    """Get supported administration routes."""
    return _assessor.get_routes()


@router.get("/health")
async def health_check():
    """ToxPath service health check."""
    return {
        "service": "toxpath",
        "status": "healthy",
        "rdkit_available": _assessor._rdkit_available,
        "assessments_cached": len(_assessments_store)
    }


# ============================================================================
# Helper Functions
# ============================================================================

def _reconstruct_internal_response(stored: ToxPathAssessResponse):
    """Reconstruct internal ToxPathResponse from stored API response.
    
    Note: In production, store the internal objects directly or use proper serialization.
    """
    from ..services.toxpath import (
        ToxPathResponse as ToxPathResponseInternal,
        RiskSummary,
        RiskTier,
        StructuralAlert,
        AlertSeverity,
        TestingStep
    )
    
    # Reconstruct risk summary
    risk_data = stored.risk_summary
    risk_summary = RiskSummary(
        overall_tier=RiskTier(risk_data.overall_tier),
        tier_rationale=risk_data.tier_rationale,
        top_risks=risk_data.top_risks,
        key_unknowns=risk_data.key_unknowns,
        key_assumptions=risk_data.key_assumptions,
        route_specific_concerns=risk_data.route_specific_concerns
    )
    
    # Reconstruct alerts
    alerts = []
    for a in stored.alerts:
        alerts.append(StructuralAlert(
            code=a.code,
            name=a.name,
            severity=AlertSeverity(a.severity),
            description=a.description,
            mechanism=a.mechanism,
            affected_organs=a.affected_organs
        ))
    
    # Reconstruct testing plan
    testing_plan = []
    for t in stored.testing_plan:
        testing_plan.append(TestingStep(
            order=t.order,
            test_name=t.test_name,
            test_type=t.test_type,
            rationale=t.rationale,
            estimated_cost_range=t.estimated_cost_range,
            timeline=t.timeline,
            required=t.required,
            gmp_glp_required=t.gmp_glp_required
        ))
    
    return ToxPathResponseInternal(
        toxpath_assessment_id=stored.toxpath_assessment_id,
        compound_ref=stored.compound_ref,
        compound_name=stored.compound_name,
        route=stored.route,
        risk_summary=risk_summary,
        alerts=alerts,
        testing_plan=testing_plan,
        consultation_required=stored.consultation_required,
        consultation_flags=stored.consultation_flags,
        properties_snapshot=stored.properties_snapshot,
        status=stored.status,
        error=stored.error,
        assessed_at=stored.assessed_at
    )
