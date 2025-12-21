"""
PatentPath Lite API Router - NeuroBotanica Week 7-8
API endpoints for patent analysis features.

Endpoints:
- POST /prior-art/search - Search prior art
- POST /novelty/assess - Assess novelty
- POST /fto/check - FTO analysis
- POST /claims/generate - Generate patent claims
- POST /tk/check - TK attribution check (Week 8)
- POST /costs/estimate - Cost estimation (Week 8)
- POST /analyze-compound - Full compound analysis
"""
from typing import Dict, List, Optional
from fastapi import APIRouter, HTTPException, Query
from pydantic import BaseModel, Field
import logging

from ..services.patentpath import (
    PriorArtSearcher,
    NoveltyScorer,
    FTOChecker,
    ClaimGenerator,
    ClaimType,
    FTO_DISCLAIMER,
    CLAIM_DISCLAIMER,
    TKAttributionChecker,
    CostEstimator,
    ProtectionStrategy
)

logger = logging.getLogger(__name__)
router = APIRouter(prefix="/patentpath", tags=["PatentPath Lite"])


# ============================================================================
# Request/Response Models
# ============================================================================

class PriorArtSearchRequest(BaseModel):
    """Request model for prior art search."""
    keywords: List[str] = Field(..., min_length=1, description="Search keywords")
    cpc_codes: Optional[List[str]] = Field(None, description="CPC classification codes")
    date_from: Optional[str] = Field(None, description="Start date (YYYY-MM-DD)")
    date_to: Optional[str] = Field(None, description="End date (YYYY-MM-DD)")
    template: Optional[str] = Field(None, description="Search template name")
    max_results: int = Field(30, ge=1, le=100)


class NoveltyAssessRequest(BaseModel):
    """Request model for novelty assessment."""
    title: str = Field(..., min_length=3, description="Innovation title")
    description: str = Field(..., min_length=10, description="Innovation description")
    keywords: Optional[List[str]] = Field(None, description="Key terms")
    claims: Optional[List[str]] = Field(None, description="Draft claims")
    smiles: Optional[str] = Field(None, description="SMILES structure")
    cpc_codes: Optional[List[str]] = Field(None, description="Relevant CPC codes")


class FTOCheckRequest(BaseModel):
    """Request model for FTO analysis."""
    compound_description: str = Field(..., min_length=10)
    keywords: List[str] = Field(..., min_length=1)
    target_market: str = Field("US", description="Target market (US only for MVP)")
    intended_use: Optional[str] = Field(None, description="Therapeutic indication")
    smiles: Optional[str] = Field(None, description="SMILES structure")


class ClaimGenerationRequest(BaseModel):
    """Request model for claim generation."""
    compound_name: str = Field(..., min_length=2)
    compound_smiles: Optional[str] = Field(None)
    therapeutic_use: Optional[str] = Field(None)
    key_features: Optional[List[str]] = Field(None)
    claim_types: Optional[List[str]] = Field(None)
    custom_params: Optional[Dict[str, str]] = Field(None)


class DimerClaimRequest(BaseModel):
    """Request model for dimer-specific claims."""
    dimer_name: str = Field(..., min_length=2)
    parent_1: str = Field(..., description="First parent cannabinoid")
    parent_2: str = Field(..., description="Second parent cannabinoid")
    linker_type: str = Field(..., description="Linker type (ester, ether, amide, alkyl, peg)")
    therapeutic_use: Optional[str] = Field(None)
    smiles: Optional[str] = Field(None)


class TKCheckRequest(BaseModel):
    """Request model for TK attribution check."""
    compound_description: str = Field(..., min_length=10, description="Description of compound")
    derivation_source: Optional[str] = Field(None, description="Source type (e.g., ethnobotanical_literature)")
    source_description: Optional[str] = Field(None, description="Detailed source description")
    source_region: Optional[str] = Field(None, description="Geographic region of source")
    community_name: Optional[str] = Field(None, description="Name of community if known")


class CostEstimateRequest(BaseModel):
    """Request model for cost estimation."""
    num_compounds: int = Field(..., ge=1, le=100, description="Number of compounds")
    jurisdictions: List[str] = Field(default=["US"], description="Filing jurisdictions")
    claim_count_avg: int = Field(15, ge=1, le=100, description="Average claims per patent")
    independent_claims_avg: int = Field(3, ge=1, le=20, description="Average independent claims")
    complexity: str = Field("moderate", description="simple, moderate, or complex")
    strategy: str = Field("patent", description="patent, trade_secret, or hybrid")
    entity_size: str = Field("large", description="large, small, or micro")
    expected_office_actions: int = Field(2, ge=0, le=10, description="Expected office actions")


class ComprehensiveAnalysisRequest(BaseModel):
    """Request model for full compound IP analysis."""
    compound_name: str = Field(..., min_length=2)
    compound_description: str = Field(..., min_length=10)
    keywords: List[str] = Field(..., min_length=1)
    therapeutic_use: Optional[str] = Field(None)
    smiles: Optional[str] = Field(None)
    derivation_source: Optional[str] = Field(None)
    source_description: Optional[str] = Field(None)
    source_region: Optional[str] = Field(None)
    community_name: Optional[str] = Field(None)
    generate_claims: bool = Field(True)
    include_cost_estimate: bool = Field(True)


# ============================================================================
# Prior Art Endpoints
# ============================================================================

@router.post("/prior-art/search")
async def search_prior_art(request: PriorArtSearchRequest):
    """Search for prior art patents.
    
    Searches USPTO PatentsView API for relevant prior art.
    Supports keyword search, CPC code filtering, and date ranges.
    """
    try:
        searcher = PriorArtSearcher()
        
        if request.template:
            results = await searcher.search(
                template=request.template,
                max_results=request.max_results
            )
        else:
            results = await searcher.search(
                keywords=request.keywords,
                cpc_codes=request.cpc_codes,
                date_from=request.date_from,
                date_to=request.date_to,
                max_results=request.max_results
            )
        
        return {
            "status": "success",
            "count": len(results),
            "results": [r.to_dict() for r in results],
            "search_params": {
                "keywords": request.keywords,
                "template": request.template,
                "max_results": request.max_results
            }
        }
        
    except Exception as e:
        logger.error(f"Prior art search error: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/prior-art/templates")
async def get_search_templates():
    """Get available prior art search templates."""
    searcher = PriorArtSearcher()
    templates = searcher.get_template_names()
    
    template_details = []
    for name in templates:
        template = searcher.get_template(name)
        if template:
            template_details.append({
                "name": name,
                "keywords": template.keywords,
                "cpc_codes": template.cpc_codes,
                "max_results": template.max_results
            })
    
    return {
        "templates": template_details,
        "count": len(template_details)
    }


@router.post("/prior-art/compound/{compound_name}")
async def search_compound_prior_art(
    compound_name: str,
    smiles: Optional[str] = None,
    include_related: bool = True
):
    """Search prior art for a specific compound."""
    try:
        searcher = PriorArtSearcher()
        results = await searcher.search_for_compound(
            compound_name=compound_name,
            smiles=smiles,
            include_related=include_related
        )
        
        return {
            "status": "success",
            "compound": compound_name,
            "count": len(results),
            "results": [r.to_dict() for r in results]
        }
        
    except Exception as e:
        logger.error(f"Compound prior art search error: {e}")
        raise HTTPException(status_code=500, detail=str(e))


# ============================================================================
# Novelty Assessment Endpoints
# ============================================================================

@router.post("/novelty/assess")
async def assess_novelty(request: NoveltyAssessRequest):
    """Assess novelty of an innovation against prior art.
    
    Returns novelty scores, risk assessment, and recommendations.
    """
    try:
        scorer = NoveltyScorer()
        report = await scorer.assess_novelty(
            title=request.title,
            description=request.description,
            keywords=request.keywords,
            claims=request.claims,
            smiles=request.smiles,
            cpc_codes=request.cpc_codes
        )
        
        return {
            "status": "success",
            "report": report.to_dict()
        }
        
    except Exception as e:
        logger.error(f"Novelty assessment error: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/novelty/summary/{title}")
async def quick_novelty_summary(
    title: str,
    description: str = Query(..., min_length=10)
):
    """Quick novelty check for an innovation."""
    try:
        from ..services.patentpath.novelty import quick_novelty_check
        
        result = await quick_novelty_check(title, description)
        return {
            "status": "success",
            **result
        }
        
    except Exception as e:
        logger.error(f"Quick novelty check error: {e}")
        raise HTTPException(status_code=500, detail=str(e))


# ============================================================================
# FTO (Freedom-to-Operate) Endpoints
# ============================================================================

@router.post("/fto/check")
async def check_fto(request: FTOCheckRequest):
    """Perform Freedom-to-Operate analysis.
    
    Analyzes blocking patents, estimates licensing costs,
    and provides risk assessment with mitigation strategies.
    
    **DISCLAIMER**: This is decision support only. Not legal advice.
    Consult qualified patent counsel before filing or commercial decisions.
    """
    try:
        checker = FTOChecker()
        report = await checker.check_fto(
            compound_description=request.compound_description,
            keywords=request.keywords,
            target_market=request.target_market,
            intended_use=request.intended_use,
            smiles=request.smiles
        )
        
        return {
            "status": "success",
            "report": report.to_dict()
        }
        
    except Exception as e:
        logger.error(f"FTO check error: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/fto/quick-check")
async def quick_fto_check(
    keywords: List[str] = Query(..., min_length=1),
    intended_use: Optional[str] = None
):
    """Quick FTO assessment for a compound."""
    try:
        from ..services.patentpath.fto_checker import quick_fto_check
        
        result = await quick_fto_check(keywords, intended_use)
        return {
            "status": "success",
            **result
        }
        
    except Exception as e:
        logger.error(f"Quick FTO check error: {e}")
        raise HTTPException(status_code=500, detail=str(e))


# ============================================================================
# Claim Generation Endpoints
# ============================================================================

@router.post("/claims/generate")
async def generate_claims(request: ClaimGenerationRequest):
    """Generate patent claims from templates.
    
    Creates USPTO-formatted claims for various claim types:
    - Composition of matter
    - Method of use
    - Pharmaceutical composition
    - Method of synthesis
    
    **DISCLAIMER**: These claims are template-generated for planning purposes only.
    Consult a patent attorney for filing-ready claims.
    """
    try:
        generator = ClaimGenerator()
        report = generator.generate_claims(
            compound_name=request.compound_name,
            compound_smiles=request.compound_smiles,
            therapeutic_use=request.therapeutic_use,
            key_features=request.key_features,
            claim_types=request.claim_types,
            custom_params=request.custom_params
        )
        
        return {
            "status": "success",
            "report": report.to_dict()
        }
        
    except Exception as e:
        logger.error(f"Claim generation error: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/claims/generate-dimer")
async def generate_dimer_claims(request: DimerClaimRequest):
    """Generate patent claims specifically for dimer compounds.
    
    Creates specialized claims for dimeric cannabinoid compounds
    with linker-specific claim language.
    """
    try:
        generator = ClaimGenerator()
        report = generator.generate_dimer_claims(
            dimer_name=request.dimer_name,
            parent_1=request.parent_1,
            parent_2=request.parent_2,
            linker_type=request.linker_type,
            therapeutic_use=request.therapeutic_use,
            smiles=request.smiles
        )
        
        return {
            "status": "success",
            "report": report.to_dict()
        }
        
    except Exception as e:
        logger.error(f"Dimer claim generation error: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/claims/types")
async def get_claim_types():
    """Get available claim types for generation."""
    generator = ClaimGenerator()
    claim_types = generator.get_available_claim_types()
    
    type_descriptions = {
        "composition_of_matter": "Claims covering the compound itself",
        "method_of_use": "Claims covering therapeutic methods",
        "pharmaceutical_composition": "Claims covering formulations",
        "method_of_synthesis": "Claims covering synthesis processes",
        "dimer_composition": "Specialized claims for dimeric compounds"
    }
    
    return {
        "claim_types": [
            {
                "type": ct,
                "description": type_descriptions.get(ct, "Patent claim type")
            }
            for ct in claim_types
        ],
        "disclaimer": CLAIM_DISCLAIMER
    }


# ============================================================================
# TK Attribution Endpoints (Week 8 - Feature 5)
# ============================================================================

@router.post("/tk/check")
async def check_tk_attribution(request: TKCheckRequest):
    """Check if compound requires traditional knowledge attribution.
    
    Analyzes compound description and sources for traditional knowledge
    indicators that may require:
    - Nagoya Protocol compliance
    - Benefit-sharing agreements
    - Community consent (Prior Informed Consent)
    - Attribution in patent specification
    
    Also detects sacred knowledge (absolute bar to patenting).
    """
    try:
        checker = TKAttributionChecker()
        result = checker.check_tk_attribution(
            compound_description=request.compound_description,
            derivation_source=request.derivation_source,
            source_description=request.source_description,
            source_region=request.source_region,
            community_name=request.community_name
        )
        
        return {
            "status": "success",
            "result": result.to_dict()
        }
        
    except Exception as e:
        logger.error(f"TK check error: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/tk/source-types")
async def get_tk_source_types():
    """Get list of recognized TK source types."""
    checker = TKAttributionChecker()
    return {
        "source_types": checker.get_tk_source_types(),
        "description": "Source types that indicate traditional knowledge derivation"
    }


@router.get("/tk/high-abs-regions")
async def get_high_abs_regions():
    """Get list of regions with strong ABS (Access and Benefit Sharing) requirements."""
    checker = TKAttributionChecker()
    return {
        "regions": checker.get_high_abs_regions(),
        "description": "Regions requiring national ABS permits for TK-derived inventions"
    }


@router.post("/tk/validate-completeness")
async def validate_tk_completeness(
    community_name: str,
    consent_obtained: bool = False,
    benefit_sharing_configured: bool = False,
    attribution_drafted: bool = False,
    abs_permit_obtained: bool = False,
    source_region: Optional[str] = None
):
    """Validate that all TK attribution requirements are met.
    
    Checks whether all necessary steps have been completed before
    filing a patent application for TK-derived compounds.
    """
    checker = TKAttributionChecker()
    result = checker.validate_attribution_completeness(
        community_name=community_name,
        consent_obtained=consent_obtained,
        benefit_sharing_configured=benefit_sharing_configured,
        attribution_drafted=attribution_drafted,
        abs_permit_obtained=abs_permit_obtained,
        source_region=source_region
    )
    
    return {
        "status": "success",
        "validation": result
    }


# ============================================================================
# Cost Estimation Endpoints (Week 8 - Feature 6)
# ============================================================================

@router.post("/costs/estimate")
async def estimate_costs(request: CostEstimateRequest):
    """Estimate patent filing and maintenance costs.
    
    Provides detailed cost projections for:
    - USPTO filing fees (with entity size discounts)
    - Attorney costs (drafting and prosecution)
    - Maintenance fees over 20-year patent life
    - Trade secret protection alternatives
    
    Supports patent, trade secret, and hybrid strategies.
    """
    try:
        estimator = CostEstimator()
        estimate = estimator.estimate_costs(
            num_compounds=request.num_compounds,
            jurisdictions=request.jurisdictions,
            claim_count_avg=request.claim_count_avg,
            independent_claims_avg=request.independent_claims_avg,
            complexity=request.complexity,
            strategy=request.strategy,
            entity_size=request.entity_size,
            expected_office_actions=request.expected_office_actions
        )
        
        return {
            "status": "success",
            "estimate": estimate.to_dict()
        }
        
    except Exception as e:
        logger.error(f"Cost estimation error: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/costs/compare-strategies")
async def compare_protection_strategies(
    num_compounds: int,
    complexity: str = "moderate",
    entity_size: str = "large"
):
    """Compare patent vs. trade secret protection costs.
    
    Returns side-by-side comparison of:
    - Patent protection costs
    - Trade secret costs
    - Hybrid approach costs
    - Recommendations based on portfolio size
    """
    try:
        estimator = CostEstimator()
        comparison = estimator.compare_strategies(
            num_compounds=num_compounds,
            complexity=complexity,
            entity_size=entity_size
        )
        
        return {
            "status": "success",
            "comparison": comparison
        }
        
    except Exception as e:
        logger.error(f"Strategy comparison error: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/costs/entity-discounts")
async def get_entity_discounts():
    """Get information about USPTO entity size discounts.
    
    Small and micro entities receive significant fee reductions:
    - Small entity: 50% discount
    - Micro entity: 75% discount
    """
    estimator = CostEstimator()
    return {
        "entity_discounts": estimator.get_entity_discount_info(),
        "note": "Verify eligibility with patent counsel before claiming small/micro entity status"
    }


@router.get("/costs/strategies")
async def get_protection_strategies():
    """Get available IP protection strategies."""
    return {
        "strategies": [
            {
                "strategy": "patent",
                "description": "Full patent protection with 20-year exclusivity",
                "best_for": "Novel compound structures, therapeutic uses"
            },
            {
                "strategy": "trade_secret",
                "description": "Indefinite protection without disclosure",
                "best_for": "Manufacturing processes, formulations"
            },
            {
                "strategy": "hybrid",
                "description": "Patent core innovations, trade secret for processes",
                "best_for": "Large portfolios seeking cost optimization"
            }
        ]
    }


# ============================================================================
# Combined Analysis Endpoints
# ============================================================================

@router.post("/analyze/comprehensive")
async def comprehensive_patent_analysis(request: ComprehensiveAnalysisRequest):
    """Comprehensive patent analysis combining all 6 PatentPath features.
    
    Performs complete IP analysis:
    1. Prior art search
    2. Novelty assessment
    3. FTO analysis
    4. Claim generation
    5. TK attribution check
    6. Cost estimation
    
    Returns complete patent landscape analysis with recommendations.
    """
    try:
        results = {
            "compound_name": request.compound_name,
            "analysis_type": "comprehensive_6_feature"
        }
        
        # 1. Prior Art Search
        searcher = PriorArtSearcher()
        prior_art = await searcher.search(keywords=request.keywords, max_results=30)
        results["prior_art"] = {
            "count": len(prior_art),
            "top_results": [r.to_dict() for r in prior_art[:5]]
        }
        
        # 2. Novelty Assessment
        scorer = NoveltyScorer(searcher=searcher)
        novelty_report = await scorer.assess_novelty(
            title=request.compound_name,
            description=request.compound_description,
            keywords=request.keywords,
            smiles=request.smiles
        )
        results["novelty"] = novelty_report.to_dict()
        
        # 3. FTO Analysis
        fto_checker = FTOChecker(searcher=searcher)
        fto_report = await fto_checker.check_fto(
            compound_description=request.compound_description,
            keywords=request.keywords,
            intended_use=request.therapeutic_use,
            smiles=request.smiles
        )
        results["fto"] = fto_report.to_dict()
        
        # 4. Claim Generation (optional)
        if request.generate_claims:
            generator = ClaimGenerator()
            claims_report = generator.generate_claims(
                compound_name=request.compound_name,
                compound_smiles=request.smiles,
                therapeutic_use=request.therapeutic_use
            )
            results["claims"] = claims_report.to_dict()
        
        # 5. TK Attribution Check
        tk_checker = TKAttributionChecker()
        tk_result = tk_checker.check_tk_attribution(
            compound_description=request.compound_description,
            derivation_source=request.derivation_source,
            source_description=request.source_description,
            source_region=request.source_region,
            community_name=request.community_name
        )
        results["tk_attribution"] = tk_result.to_dict()
        
        # 6. Cost Estimation (optional)
        if request.include_cost_estimate:
            estimator = CostEstimator()
            claim_count = len(results.get("claims", {}).get("generated_claims", {})) * 3 or 15
            cost_estimate = estimator.estimate_costs(
                num_compounds=1,
                claim_count_avg=claim_count
            )
            results["cost_estimate"] = cost_estimate.to_dict()
        
        # Generate comprehensive summary
        results["summary"] = {
            "novelty_score": novelty_report.overall_novelty_score,
            "novelty_level": novelty_report.novelty_level.value,
            "fto_risk_level": fto_report.fto_risk_level.value,
            "blocking_patents": fto_report.num_blocking_patents,
            "tk_derived": tk_result.tk_derived,
            "tk_attribution_required": tk_result.attribution_required,
            "sacred_knowledge_detected": tk_result.sacred_knowledge_detected,
            "can_proceed": tk_result.can_proceed,
            "recommendation": _generate_full_recommendation(
                novelty_report, fto_report, tk_result
            )
        }
        
        if request.include_cost_estimate:
            results["summary"]["estimated_total_cost"] = cost_estimate.total_cost_estimate.total_cost
        
        results["disclaimer"] = FTO_DISCLAIMER
        
        return {
            "status": "success",
            "analysis": results
        }
        
    except Exception as e:
        logger.error(f"Comprehensive analysis error: {e}")
        raise HTTPException(status_code=500, detail=str(e))


def _generate_full_recommendation(novelty_report, fto_report, tk_result) -> str:
    """Generate recommendation based on all 6 features."""
    # Check TK blocking first
    if not tk_result.can_proceed:
        return "DO NOT PROCEED: Sacred traditional knowledge detected - cannot be patented."
    
    # Check TK attribution requirement
    tk_warning = ""
    if tk_result.attribution_required:
        tk_warning = " TK ATTRIBUTION REQUIRED: Establish consent and benefit-sharing before filing."
    
    novelty_level = novelty_report.novelty_level.value
    fto_risk = fto_report.fto_risk_level.value
    
    if novelty_level == "high" and fto_risk == "low":
        return f"STRONG POSITION: High novelty with low FTO risk.{tk_warning}"
    elif novelty_level in ["high", "moderate"] and fto_risk in ["low", "moderate"]:
        return f"FAVORABLE POSITION: Consider patent filing with attorney review.{tk_warning}"
    elif fto_risk in ["high", "critical"]:
        return f"CAUTION: Significant FTO risks. Address blocking patents first.{tk_warning}"
    elif novelty_level in ["low", "blocked"]:
        return f"LIMITED POTENTIAL: Low novelty may impact patentability.{tk_warning}"
    else:
        return f"MIXED SIGNALS: Recommend detailed attorney consultation.{tk_warning}"


def _generate_comprehensive_recommendation(novelty_report, fto_report) -> str:
    """Generate recommendation based on combined analysis (legacy 4-feature)."""
    novelty_level = novelty_report.novelty_level.value
    fto_risk = fto_report.fto_risk_level.value
    
    if novelty_level == "high" and fto_risk == "low":
        return "STRONG POSITION: High novelty with low FTO risk. Recommend proceeding with patent filing."
    elif novelty_level in ["high", "moderate"] and fto_risk in ["low", "moderate"]:
        return "FAVORABLE POSITION: Consider patent filing with attorney review of FTO concerns."
    elif fto_risk in ["high", "critical"]:
        return "CAUTION: Significant FTO risks identified. Address blocking patents before proceeding."
    elif novelty_level in ["low", "blocked"]:
        return "LIMITED POTENTIAL: Low novelty may impact patentability. Consider design modifications."
    else:
        return "MIXED SIGNALS: Recommend detailed attorney consultation before proceeding."

