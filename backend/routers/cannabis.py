"""
Cannabis API Router

REST API endpoints for cannabis-specific analysis:
- Cannabinoid/terpene entourage effect analysis
- Safety assessment and abuse liability
- Schedule III regulatory strategy

Trade Secret: Exposes capabilities while protecting proprietary algorithms.
"""

from fastapi import APIRouter, HTTPException
from pydantic import BaseModel, Field
from typing import Dict, List, Optional
from datetime import date

# Import cannabis modules
from backend.services.chempath import (
    CannabinoidAnalyzer,
    EntourageProfile,
    FormulationRecommendation,
)
from backend.services.toxpath import (
    CannabisSafetyAssessor,
    SafetyProfile,
    AbuseLiabilityProfile,
    DrugInteraction,
    SafetyRecommendation,
)
from backend.services.regpath import (
    Schedule3Strategist,
    DEARegistrationPlan,
    BotanicalDrugPathway,
    DEASchedule,
    DEALicenseType,
)


# =============================================================================
# Request/Response Models
# =============================================================================

class EntourageAnalysisRequest(BaseModel):
    """Request for entourage effect analysis."""
    cannabinoids: Dict[str, float] = Field(
        ...,
        description="Cannabinoid concentrations (%). Keys: THC, CBD, CBG, CBN, etc.",
        example={"THC": 18.5, "CBD": 0.5, "CBG": 0.3, "CBN": 0.1}
    )
    terpenes: Dict[str, float] = Field(
        ...,
        description="Terpene concentrations (%). Keys: myrcene, limonene, etc.",
        example={"myrcene": 0.8, "limonene": 0.4, "beta-caryophyllene": 0.3}
    )
    formulation_name: Optional[str] = Field(
        None,
        description="Optional name for the formulation"
    )


class EntourageAnalysisResponse(BaseModel):
    """Response with entourage effect analysis."""
    profile_id: str
    formulation_name: str
    overall_entourage_score: float = Field(
        ..., description="0-100 scale, higher = stronger entourage effect"
    )
    synergy_count: int
    primary_indication: Optional[str]
    secondary_indications: List[str]
    receptor_activity: Dict[str, float]
    therapeutic_predictions: Dict[str, float]
    tk_attributed: bool
    tk_communities: List[str]
    confidence_score: float


class FormulationOptimizationRequest(BaseModel):
    """Request for formulation optimization."""
    profile_id: str = Field(..., description="Profile ID from previous analysis")
    target_indication: str = Field(
        ...,
        description="Target therapeutic indication",
        example="analgesic"
    )


class SafetyAssessmentRequest(BaseModel):
    """Request for cannabis safety assessment."""
    cannabinoids: Dict[str, float] = Field(
        ...,
        description="Cannabinoid concentrations (%)",
        example={"THC": 15.0, "CBD": 5.0}
    )
    formulation_name: str = Field(..., description="Name of formulation")
    formulation_type: str = Field(
        "oral",
        description="Delivery method: oral, inhalation, sublingual, topical"
    )
    concomitant_medications: Optional[List[str]] = Field(
        None,
        description="List of current medications to check for interactions"
    )
    populations: Optional[List[str]] = Field(
        None,
        description="Populations to assess: pediatric, geriatric, pregnant, etc."
    )


class SafetyAssessmentResponse(BaseModel):
    """Response with safety assessment."""
    profile_id: str
    formulation_name: str
    
    # Abuse liability
    proposed_schedule: str
    overall_abuse_liability: float
    schedule_justification: str
    
    # Drug interactions
    drug_interactions: List[Dict]
    high_risk_interactions: int
    
    # Adverse events
    adverse_event_count: int
    severe_events: List[str]
    
    # Population risks
    population_risks: Dict[str, str]
    
    # Overall safety
    overall_safety_score: float
    monitoring_required: List[str]
    
    # Contraindications
    absolute_contraindications: List[str]
    relative_contraindications: List[str]


class AbuseLiabilityRequest(BaseModel):
    """Request for abuse liability assessment only."""
    cannabinoids: Dict[str, float]
    formulation_type: str = "oral"


class DEARegistrationRequest(BaseModel):
    """Request for DEA registration plan."""
    proposed_schedule: str = Field(
        "Schedule III",
        description="Target schedule: Schedule II, Schedule III, etc."
    )
    license_types: Optional[List[str]] = Field(
        None,
        description="License types: Manufacturer, Distributor, Researcher"
    )
    estimated_annual_production_kg: float = Field(
        100.0,
        description="Projected annual production in kilograms"
    )


class DEARegistrationResponse(BaseModel):
    """Response with DEA registration plan."""
    plan_id: str
    proposed_schedule: str
    license_types: List[str]
    requirements_count: int
    total_timeline_months: int
    estimated_cost_min: int
    estimated_cost_max: int
    annual_compliance_cost_min: int
    annual_compliance_cost_max: int
    requires_quota: bool
    quota_deadline: Optional[str]
    critical_path_items: List[str]
    security_requirements: List[str]


class BotanicalPathwayRequest(BaseModel):
    """Request for botanical drug pathway design."""
    product_name: str = Field(..., description="Product name")
    botanical_source: str = Field(
        "Cannabis sativa L.",
        description="Scientific name of botanical source"
    )
    active_constituents: Optional[List[str]] = Field(
        None,
        description="Key active compounds"
    )
    target_indication: str = Field(
        "chronic pain",
        description="Primary therapeutic indication"
    )
    start_date: Optional[str] = Field(
        None,
        description="Development start date (YYYY-MM-DD)"
    )


class BotanicalPathwayResponse(BaseModel):
    """Response with botanical drug pathway."""
    pathway_id: str
    product_description: str
    botanical_source: str
    total_development_years: float
    estimated_cost_min_millions: float
    estimated_cost_max_millions: float
    probability_of_success: float
    milestones_count: int
    clinical_trials_count: int
    risk_factors: List[str]
    mitigation_strategies: List[str]


class TimeToMarketRequest(BaseModel):
    """Request for time to market estimate."""
    current_phase: str = Field(
        "Pre-clinical",
        description="Current development phase"
    )


# =============================================================================
# Router
# =============================================================================

router = APIRouter(
    prefix="/api/v1/cannabis",
    tags=["Cannabis Analysis"],
    responses={404: {"description": "Not found"}},
)

# Initialize analyzers
cannabinoid_analyzer = CannabinoidAnalyzer()
safety_assessor = CannabisSafetyAssessor()
schedule3_strategist = Schedule3Strategist()


# =============================================================================
# Entourage Effect Endpoints
# =============================================================================

@router.post(
    "/entourage/analyze",
    response_model=EntourageAnalysisResponse,
    summary="Analyze entourage effect",
    description="""
    Predict synergistic effects based on cannabinoid-terpene interactions.
    
    Trade Secret: Uses proprietary synergy algorithms derived from indigenous
    cannabis knowledge and receptor binding predictions.
    """
)
async def analyze_entourage_effect(request: EntourageAnalysisRequest):
    """Analyze cannabinoid/terpene entourage effect."""
    try:
        profile = cannabinoid_analyzer.analyze_entourage_effect(
            cannabinoids=request.cannabinoids,
            terpenes=request.terpenes,
            formulation_name=request.formulation_name,
        )
        
        return EntourageAnalysisResponse(
            profile_id=profile.profile_id,
            formulation_name=profile.formulation_name,
            overall_entourage_score=profile.overall_entourage_score,
            synergy_count=len(profile.synergy_interactions),
            primary_indication=profile.primary_indication,
            secondary_indications=profile.secondary_indications,
            receptor_activity=profile.receptor_activity,
            therapeutic_predictions=profile.therapeutic_predictions,
            tk_attributed=profile.tk_attributed,
            tk_communities=profile.tk_communities,
            confidence_score=profile.confidence_score,
        )
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router.post(
    "/entourage/optimize",
    response_model=Dict,
    summary="Optimize formulation",
    description="Generate recommendations to optimize formulation for target indication."
)
async def optimize_formulation(request: FormulationOptimizationRequest):
    """Get formulation optimization recommendations."""
    profile = cannabinoid_analyzer.get_analysis(request.profile_id)
    if not profile:
        raise HTTPException(status_code=404, detail="Profile not found")
    
    try:
        recommendation = cannabinoid_analyzer.get_formulation_recommendations(
            current_profile=profile,
            target_indication=request.target_indication,
        )
        
        return {
            "recommendation_id": recommendation.recommendation_id,
            "target_indication": recommendation.target_indication,
            "current_entourage_score": recommendation.current_entourage_score,
            "projected_entourage_score": recommendation.projected_entourage_score,
            "improvement_percentage": recommendation.improvement_percentage,
            "cannabinoid_adjustments": recommendation.cannabinoid_adjustments,
            "terpene_adjustments": recommendation.terpene_adjustments,
            "rationale": recommendation.rationale,
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router.get(
    "/entourage/{profile_id}",
    response_model=Dict,
    summary="Get entourage analysis",
    description="Retrieve a previously generated entourage analysis."
)
async def get_entourage_analysis(profile_id: str):
    """Get cached entourage analysis."""
    profile = cannabinoid_analyzer.get_analysis(profile_id)
    if not profile:
        raise HTTPException(status_code=404, detail="Profile not found")
    
    return {
        "profile_id": profile.profile_id,
        "formulation_name": profile.formulation_name,
        "overall_entourage_score": profile.overall_entourage_score,
        "synergy_interactions": [
            {
                "compound_a": i.compound_a,
                "compound_b": i.compound_b,
                "interaction_type": i.interaction_type,
                "synergy_coefficient": i.synergy_coefficient,
            }
            for i in profile.synergy_interactions
        ],
        "therapeutic_predictions": profile.therapeutic_predictions,
        "primary_indication": profile.primary_indication,
        "tk_attributed": profile.tk_attributed,
        "tk_communities": profile.tk_communities,
    }


# =============================================================================
# Safety Assessment Endpoints
# =============================================================================

@router.post(
    "/safety/assess",
    response_model=SafetyAssessmentResponse,
    summary="Comprehensive safety assessment",
    description="""
    Generate complete safety profile including abuse liability,
    drug interactions, and population-specific risks.
    
    Trade Secret: Proprietary safety scoring algorithms.
    """
)
async def assess_safety(request: SafetyAssessmentRequest):
    """Generate comprehensive safety assessment."""
    try:
        profile = safety_assessor.generate_safety_profile(
            cannabinoids=request.cannabinoids,
            formulation_name=request.formulation_name,
            formulation_type=request.formulation_type,
            concomitant_medications=request.concomitant_medications,
            populations=request.populations,
        )
        
        # Identify severe events
        severe_events = [
            e.name for e in profile.adverse_events
            if e.severity.value in ["high", "severe"]
        ]
        
        return SafetyAssessmentResponse(
            profile_id=profile.profile_id,
            formulation_name=profile.formulation_name,
            proposed_schedule=profile.abuse_liability.proposed_schedule.value,
            overall_abuse_liability=profile.abuse_liability.overall_abuse_liability,
            schedule_justification=profile.abuse_liability.schedule_justification,
            drug_interactions=[
                {
                    "drug": i.drug_name,
                    "class": i.drug_class,
                    "risk_level": i.risk_level.value,
                    "recommendation": i.recommendation,
                }
                for i in profile.drug_interactions
            ],
            high_risk_interactions=profile.high_risk_interactions,
            adverse_event_count=len(profile.adverse_events),
            severe_events=severe_events,
            population_risks={k: v.value for k, v in profile.population_risks.items()},
            overall_safety_score=profile.overall_safety_score,
            monitoring_required=profile.monitoring_parameters,
            absolute_contraindications=profile.absolute_contraindications,
            relative_contraindications=profile.relative_contraindications,
        )
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router.post(
    "/safety/abuse-liability",
    response_model=Dict,
    summary="Abuse liability assessment",
    description="Assess abuse liability for DEA scheduling determination."
)
async def assess_abuse_liability(request: AbuseLiabilityRequest):
    """Assess abuse liability only."""
    try:
        profile = safety_assessor.assess_abuse_liability(
            cannabinoids=request.cannabinoids,
            formulation_type=request.formulation_type,
        )
        
        return {
            "profile_id": profile.profile_id,
            "proposed_schedule": profile.proposed_schedule.value,
            "schedule_justification": profile.schedule_justification,
            "overall_abuse_liability": profile.overall_abuse_liability,
            "components": {
                "dependence_potential": profile.dependence_potential,
                "reinforcement_potential": profile.reinforcement_potential,
                "euphoria_potential": profile.euphoria_potential,
                "withdrawal_severity": profile.withdrawal_severity,
                "tolerance_development": profile.tolerance_development,
            },
            "risk_factors": profile.risk_factors,
            "protective_factors": profile.protective_factors,
            "dea_requirements": profile.controlled_substance_act_requirements,
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router.post(
    "/safety/interactions",
    response_model=Dict,
    summary="Drug interaction check",
    description="Check for drug-drug interactions with cannabis formulation."
)
async def check_drug_interactions(
    cannabinoids: Dict[str, float],
    medications: List[str],
):
    """Check drug interactions."""
    try:
        interactions = safety_assessor.assess_drug_interactions(
            cannabinoids=cannabinoids,
            concomitant_medications=medications,
        )
        
        return {
            "total_interactions": len(interactions),
            "high_risk_count": sum(
                1 for i in interactions 
                if i.risk_level.value in ["high", "severe"]
            ),
            "interactions": [
                {
                    "drug": i.drug_name,
                    "class": i.drug_class,
                    "type": i.interaction_type.value,
                    "mechanism": i.mechanism.value,
                    "risk_level": i.risk_level.value,
                    "significance": i.clinical_significance,
                    "recommendation": i.recommendation,
                    "monitoring": i.monitoring_parameters,
                }
                for i in interactions
            ],
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router.get(
    "/safety/{profile_id}/recommendations",
    response_model=Dict,
    summary="Safety recommendations",
    description="Get safety recommendations for a assessed formulation."
)
async def get_safety_recommendations(profile_id: str):
    """Get safety recommendations."""
    profile = safety_assessor.get_assessment(profile_id)
    if not profile:
        raise HTTPException(status_code=404, detail="Profile not found")
    
    recommendations = safety_assessor.get_safety_recommendations(profile)
    
    return {
        "profile_id": profile_id,
        "formulation_name": profile.formulation_name,
        "overall_safety_score": profile.overall_safety_score,
        "recommendations": [
            {
                "category": r.category,
                "priority": r.priority,
                "recommendation": r.recommendation,
                "rationale": r.rationale,
            }
            for r in recommendations
        ],
    }


# =============================================================================
# Schedule III Regulatory Endpoints
# =============================================================================

@router.post(
    "/regulatory/dea-registration",
    response_model=DEARegistrationResponse,
    summary="DEA registration plan",
    description="""
    Generate comprehensive DEA registration plan for Schedule III cannabis.
    
    Trade Secret: Regulatory timeline optimization, cost estimation algorithms.
    """
)
async def generate_dea_registration_plan(request: DEARegistrationRequest):
    """Generate DEA registration plan."""
    try:
        # Convert schedule string to enum
        schedule_map = {
            "Schedule II": DEASchedule.SCHEDULE_II,
            "Schedule III": DEASchedule.SCHEDULE_III,
            "Schedule IV": DEASchedule.SCHEDULE_IV,
            "Schedule V": DEASchedule.SCHEDULE_V,
        }
        schedule = schedule_map.get(request.proposed_schedule, DEASchedule.SCHEDULE_III)
        
        # Convert license types
        license_types = None
        if request.license_types:
            license_map = {
                "Manufacturer": DEALicenseType.MANUFACTURER,
                "Distributor": DEALicenseType.DISTRIBUTOR,
                "Researcher": DEALicenseType.RESEARCHER,
                "Pharmacy": DEALicenseType.PHARMACY,
            }
            license_types = [license_map[lt] for lt in request.license_types if lt in license_map]
        
        plan = schedule3_strategist.generate_dea_registration_plan(
            proposed_schedule=schedule,
            license_types=license_types,
            estimated_annual_production_kg=request.estimated_annual_production_kg,
        )
        
        return DEARegistrationResponse(
            plan_id=plan.plan_id,
            proposed_schedule=plan.proposed_schedule.value,
            license_types=[lt.value for lt in plan.license_types_needed],
            requirements_count=len(plan.requirements),
            total_timeline_months=plan.total_timeline_months,
            estimated_cost_min=plan.estimated_total_cost[0],
            estimated_cost_max=plan.estimated_total_cost[1],
            annual_compliance_cost_min=plan.annual_compliance_cost[0],
            annual_compliance_cost_max=plan.annual_compliance_cost[1],
            requires_quota=plan.requires_quota,
            quota_deadline=plan.quota_application_deadline.isoformat() if plan.quota_application_deadline else None,
            critical_path_items=plan.critical_path_items,
            security_requirements=plan.security_requirements,
        )
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router.post(
    "/regulatory/botanical-pathway",
    response_model=BotanicalPathwayResponse,
    summary="FDA botanical drug pathway",
    description="""
    Design complete FDA botanical drug development pathway.
    
    Trade Secret: Pathway optimization, timeline prediction, cost estimation.
    """
)
async def design_botanical_pathway(request: BotanicalPathwayRequest):
    """Design FDA botanical drug pathway."""
    try:
        start_date = None
        if request.start_date:
            start_date = date.fromisoformat(request.start_date)
        
        pathway = schedule3_strategist.design_botanical_drug_pathway(
            product_name=request.product_name,
            botanical_source=request.botanical_source,
            active_constituents=request.active_constituents,
            target_indication=request.target_indication,
            start_date=start_date,
        )
        
        return BotanicalPathwayResponse(
            pathway_id=pathway.pathway_id,
            product_description=pathway.product_description,
            botanical_source=pathway.botanical_source,
            total_development_years=pathway.total_development_years,
            estimated_cost_min_millions=pathway.estimated_total_cost_millions[0],
            estimated_cost_max_millions=pathway.estimated_total_cost_millions[1],
            probability_of_success=pathway.probability_of_success,
            milestones_count=len(pathway.milestones),
            clinical_trials_count=len(pathway.clinical_trials),
            risk_factors=pathway.risk_factors,
            mitigation_strategies=pathway.mitigation_strategies,
        )
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router.get(
    "/regulatory/dea/{plan_id}",
    response_model=Dict,
    summary="Get DEA registration plan",
    description="Retrieve a previously generated DEA registration plan."
)
async def get_dea_plan(plan_id: str):
    """Get cached DEA registration plan."""
    plan = schedule3_strategist.get_registration_plan(plan_id)
    if not plan:
        raise HTTPException(status_code=404, detail="Plan not found")
    
    return {
        "plan_id": plan.plan_id,
        "proposed_schedule": plan.proposed_schedule.value,
        "license_types": [lt.value for lt in plan.license_types_needed],
        "requirements": [
            {
                "category": r.category,
                "requirement": r.requirement,
                "cfr_reference": r.cfr_reference,
                "timeline": r.timeline,
                "complexity": r.complexity,
                "estimated_cost": r.estimated_cost,
                "documentation_needed": r.documentation_needed,
            }
            for r in plan.requirements
        ],
        "total_timeline_months": plan.total_timeline_months,
        "estimated_total_cost": plan.estimated_total_cost,
        "requires_quota": plan.requires_quota,
        "security_requirements": plan.security_requirements,
    }


@router.get(
    "/regulatory/botanical/{pathway_id}",
    response_model=Dict,
    summary="Get botanical pathway",
    description="Retrieve a previously designed botanical drug pathway."
)
async def get_botanical_pathway(pathway_id: str):
    """Get cached botanical drug pathway."""
    pathway = schedule3_strategist.get_botanical_pathway(pathway_id)
    if not pathway:
        raise HTTPException(status_code=404, detail="Pathway not found")
    
    return {
        "pathway_id": pathway.pathway_id,
        "product_description": pathway.product_description,
        "botanical_source": pathway.botanical_source,
        "milestones": [
            {
                "name": m.name,
                "submission_type": m.submission_type.value,
                "target_date": m.target_date.isoformat() if m.target_date else None,
                "estimated_duration_days": m.estimated_duration_days,
                "required_documents": m.required_documents,
                "estimated_cost": m.estimated_cost,
            }
            for m in pathway.milestones
        ],
        "clinical_trials": [
            {
                "phase": t.phase,
                "design_type": t.design_type,
                "primary_endpoint": t.primary_endpoint,
                "sample_size": t.sample_size,
                "duration_weeks": t.duration_weeks,
                "estimated_cost": t.estimated_cost,
            }
            for t in pathway.clinical_trials
        ],
        "total_development_years": pathway.total_development_years,
        "probability_of_success": pathway.probability_of_success,
        "risk_factors": pathway.risk_factors,
    }


@router.post(
    "/regulatory/time-to-market",
    response_model=Dict,
    summary="Time to market estimate",
    description="Estimate time to market from current development phase."
)
async def estimate_time_to_market(request: TimeToMarketRequest):
    """Estimate time to market."""
    try:
        estimate = schedule3_strategist.estimate_time_to_market(
            current_phase=request.current_phase
        )
        
        return {
            "current_phase": estimate["current_phase"],
            "estimated_months_to_approval": estimate["estimated_months_to_approval"],
            "estimated_years_to_approval": round(estimate["estimated_years_to_approval"], 1),
            "remaining_phases": estimate["remaining_phases"],
            "estimated_cost_range_millions": estimate["estimated_cost_range_millions"],
            "success_probability": round(estimate["success_probability"], 4),
            "target_approval_date": estimate["target_approval_date"].isoformat(),
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


# =============================================================================
# Health Check
# =============================================================================

@router.get("/health", summary="Health check")
async def health_check():
    """Check cannabis module health."""
    return {
        "status": "healthy",
        "modules": {
            "cannabinoid_analyzer": "operational",
            "safety_assessor": "operational",
            "schedule3_strategist": "operational",
        },
        "endpoints": {
            "entourage_analysis": 3,
            "safety_assessment": 4,
            "regulatory_planning": 5,
        },
    }
