"""
RegPath API Router - Regulatory Pathway Optimization Endpoints

Exposes RegPath strategy generation while protecting proprietary algorithms.
"""

from fastapi import APIRouter, HTTPException
from pydantic import BaseModel, Field
from typing import Optional, List, Dict, Any
from datetime import datetime
import uuid

from ..services.regpath import (
    RegPathStrategist,
    RegPathRequest as RegPathRequestInternal,
    ProductProfile as ProductProfileInternal,
    EvidenceInputs,
    StrategyConstraints,
    RegPathMemoGenerator,
    MemoFormat,
    MemoConfig
)


router = APIRouter(prefix="/api/v1/regpath", tags=["RegPath - Regulatory Strategy"])


# ============================================================================
# Pydantic Request/Response Models
# ============================================================================

class ProductProfileRequest(BaseModel):
    """Product profile for regulatory assessment."""
    product_name: str = Field(..., description="Product name")
    product_type: str = Field(..., description="Product type: isolated_cannabinoid, novel_dimer, formulation, combination, botanical_extract, synthetic")
    intended_use: str = Field(..., description="Intended therapeutic use")
    route: str = Field(..., description="Administration route: oral, inhaled, sublingual, topical")
    target_markets: List[str] = Field(default=["US"], description="Target markets")
    therapeutic_area: Optional[str] = Field(None, description="Therapeutic area")
    active_ingredients: Optional[List[str]] = Field(None, description="Active ingredients")


class EvidenceInputsRequest(BaseModel):
    """Evidence inputs for strategy generation."""
    chempath_job_id: Optional[str] = Field(None, description="ChemPath analysis job ID")
    toxpath_assessment_id: Optional[str] = Field(None, description="ToxPath assessment ID")
    existing_literature_refs: Optional[List[str]] = Field(None, description="Literature references")
    predicate_drug: Optional[str] = Field(None, description="Reference product for 505(b)(2)")


class ConstraintsRequest(BaseModel):
    """Constraints for strategy generation."""
    budget_range: Optional[str] = Field(None, description="Budget range, e.g., '500K-2M'")
    launch_window: Optional[str] = Field(None, description="Target launch window, e.g., '2026 Q4'")
    risk_tolerance: str = Field("moderate", description="Risk tolerance: low, moderate, high")
    preferred_pathway: Optional[str] = Field(None, description="Preferred regulatory pathway")
    excluded_pathways: List[str] = Field(default=[], description="Pathways to exclude")


class RegPathStrategyRequest(BaseModel):
    """Request model for RegPath strategy generation."""
    product_profile: ProductProfileRequest
    evidence_inputs: Optional[EvidenceInputsRequest] = None
    constraints: Optional[ConstraintsRequest] = None
    
    class Config:
        json_schema_extra = {
            "example": {
                "product_profile": {
                    "product_name": "CBD-X Oral Solution",
                    "product_type": "isolated_cannabinoid",
                    "intended_use": "Treatment of epilepsy",
                    "route": "oral",
                    "target_markets": ["US"],
                    "therapeutic_area": "Neurology",
                    "active_ingredients": ["Cannabidiol"]
                },
                "evidence_inputs": {
                    "predicate_drug": "Epidiolex"
                },
                "constraints": {
                    "budget_range": "15M-50M",
                    "risk_tolerance": "moderate"
                }
            }
        }


class PathwayRecommendationResponse(BaseModel):
    """Pathway recommendation in response."""
    pathway: str
    pathway_full_name: str
    confidence: str
    rationale: List[str]
    key_requirements: List[str]
    estimated_cost_range: str
    estimated_timeline_months: Dict[str, int]


class ChecklistItemResponse(BaseModel):
    """Checklist item in response."""
    category: str
    item: str
    description: str
    status: str
    priority: str
    dependencies: List[str]
    estimated_effort: Optional[str]


class TimelineMilestoneResponse(BaseModel):
    """Timeline milestone in response."""
    phase: str
    milestone: str
    description: str
    duration_months: int
    start_month: int
    end_month: int
    dependencies: List[str]
    deliverables: List[str]


class RegPathStrategyResponse(BaseModel):
    """Full RegPath strategy response."""
    regpath_strategy_id: str
    product_profile: Dict[str, Any]
    novelty_assessment: str
    evidence_strength: str
    primary_pathway: PathwayRecommendationResponse
    fallback_pathway: Optional[PathwayRecommendationResponse]
    gating_questions: List[str]
    readiness_checklist: List[ChecklistItemResponse]
    timeline: List[TimelineMilestoneResponse]
    next_actions: List[str]
    key_assumptions: List[str]
    status: str
    error: Optional[str]
    generated_at: str


class MemoRequest(BaseModel):
    """Request for memo generation."""
    format: str = Field("markdown", description="Output format: markdown, html, json")
    include_timeline_table: bool = Field(True, description="Include timeline table")
    include_cost_estimates: bool = Field(True, description="Include cost estimates")
    include_checklist: bool = Field(True, description="Include readiness checklist")
    include_assumptions: bool = Field(True, description="Include assumptions")
    company_name: Optional[str] = Field(None, description="Organization name for header")
    prepared_by: Optional[str] = Field(None, description="Preparer name for header")


class MemoResponse(BaseModel):
    """Response with generated memo."""
    strategy_id: str
    format: str
    memo_content: str
    generated_at: str


class PathwayInfo(BaseModel):
    """Pathway information."""
    pathway: str
    full_name: str
    typical_timeline_min: int
    typical_timeline_max: int


class ProductTypeInfo(BaseModel):
    """Product type information."""
    type: str
    display_name: str


# ============================================================================
# In-memory storage (MVP - replace with database in production)
# ============================================================================

_strategies_store: Dict[str, RegPathStrategyResponse] = {}
_internal_store: Dict[str, Any] = {}  # Store internal objects for memo generation
_strategist = RegPathStrategist()


# ============================================================================
# Endpoints
# ============================================================================

@router.post("/strategy", response_model=RegPathStrategyResponse)
async def generate_strategy(request: RegPathStrategyRequest):
    """Generate regulatory pathway strategy for a product.
    
    Analyzes product profile, evidence base, and constraints to recommend
    primary and fallback regulatory pathways with timelines and checklists.
    
    Trade secret: Pathway decision matrices, scoring weights, and timeline
    templates are proprietary and not exposed in the response.
    """
    # Convert request to internal format
    profile = ProductProfileInternal(
        product_name=request.product_profile.product_name,
        product_type=request.product_profile.product_type,
        intended_use=request.product_profile.intended_use,
        route=request.product_profile.route,
        target_markets=request.product_profile.target_markets,
        therapeutic_area=request.product_profile.therapeutic_area,
        active_ingredients=request.product_profile.active_ingredients
    )
    
    evidence = None
    if request.evidence_inputs:
        evidence = EvidenceInputs(
            chempath_job_id=request.evidence_inputs.chempath_job_id,
            toxpath_assessment_id=request.evidence_inputs.toxpath_assessment_id,
            existing_literature_refs=request.evidence_inputs.existing_literature_refs,
            predicate_drug=request.evidence_inputs.predicate_drug
        )
    
    constraints = None
    if request.constraints:
        constraints = StrategyConstraints(
            budget_range=request.constraints.budget_range,
            launch_window=request.constraints.launch_window,
            risk_tolerance=request.constraints.risk_tolerance,
            preferred_pathway=request.constraints.preferred_pathway,
            excluded_pathways=request.constraints.excluded_pathways
        )
    
    internal_request = RegPathRequestInternal(
        product_profile=profile,
        evidence_inputs=evidence,
        constraints=constraints
    )
    
    # Generate strategy
    result = _strategist.generate_strategy(internal_request)
    
    # Convert to response format
    response_data = result.to_dict()
    response = RegPathStrategyResponse(**response_data)
    
    # Store for later retrieval
    _strategies_store[response.regpath_strategy_id] = response
    _internal_store[response.regpath_strategy_id] = result
    
    return response


@router.get("/strategies/{strategy_id}", response_model=RegPathStrategyResponse)
async def get_strategy(strategy_id: str):
    """Retrieve a previous RegPath strategy by ID."""
    if strategy_id not in _strategies_store:
        raise HTTPException(status_code=404, detail="Strategy not found")
    
    return _strategies_store[strategy_id]


@router.post("/strategies/{strategy_id}/memo", response_model=MemoResponse)
async def generate_memo(strategy_id: str, request: MemoRequest):
    """Generate a formatted regulatory strategy memo.
    
    Supports markdown, HTML, and JSON output formats.
    """
    if strategy_id not in _internal_store:
        raise HTTPException(status_code=404, detail="Strategy not found")
    
    # Get stored internal result
    internal_result = _internal_store[strategy_id]
    
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
        include_timeline_table=request.include_timeline_table,
        include_cost_estimates=request.include_cost_estimates,
        include_checklist=request.include_checklist,
        include_assumptions=request.include_assumptions,
        company_name=request.company_name,
        prepared_by=request.prepared_by
    )
    
    generator = RegPathMemoGenerator(config)
    memo_content = generator.generate(internal_result, format_enum)
    
    return MemoResponse(
        strategy_id=strategy_id,
        format=format_enum.value,
        memo_content=memo_content,
        generated_at=datetime.utcnow().isoformat()
    )


@router.get("/pathways", response_model=List[PathwayInfo])
async def get_pathways():
    """Get available regulatory pathways and their typical timelines."""
    return _strategist.get_pathways()


@router.get("/product-types", response_model=List[ProductTypeInfo])
async def get_product_types():
    """Get supported product types."""
    return _strategist.get_product_types()


@router.get("/checklist-categories", response_model=List[str])
async def get_checklist_categories():
    """Get readiness checklist categories."""
    return _strategist.get_checklist_categories()


@router.get("/health")
async def health_check():
    """RegPath service health check."""
    return {
        "service": "regpath",
        "status": "healthy",
        "strategies_cached": len(_strategies_store)
    }


# ============================================================================
# Integration Endpoints (Pipeline)
# ============================================================================

@router.post("/full-assessment")
async def full_regulatory_assessment(
    product_profile: ProductProfileRequest,
    chempath_job_id: Optional[str] = None,
    toxpath_assessment_id: Optional[str] = None
):
    """Run full regulatory assessment integrating ChemPath and ToxPath results.
    
    This endpoint pulls in prior ChemPath and ToxPath analyses to inform
    the regulatory strategy recommendation.
    """
    # Build evidence from prior analyses
    evidence = EvidenceInputsRequest(
        chempath_job_id=chempath_job_id,
        toxpath_assessment_id=toxpath_assessment_id
    )
    
    # Create full request
    request = RegPathStrategyRequest(
        product_profile=product_profile,
        evidence_inputs=evidence
    )
    
    # Generate strategy
    return await generate_strategy(request)
