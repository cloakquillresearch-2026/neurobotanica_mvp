"""
Terpene Analysis API Router - NeuroBotanica Week 7
API endpoints for terpene-cannabinoid synergy analysis.

Endpoints:
- POST /synergy/analyze - Analyze terpene-cannabinoid synergy
- GET /terpenes - List all terpenes in database
- GET /terpenes/{name} - Get terpene details
- POST /profiles/compare - Compare two terpene profiles
- POST /profiles/optimize - Get optimal profile for therapeutic goal
"""
from typing import Dict, List, Optional
from fastapi import APIRouter, HTTPException, Query
from pydantic import BaseModel, Field
import logging

from ..services.terpene_analyzer import (
    TerpeneAnalyzer,
    EvidenceTier
)

logger = logging.getLogger(__name__)
router = APIRouter(prefix="/terpenes", tags=["Terpene Analysis"])


# ============================================================================
# Request/Response Models
# ============================================================================

class SynergyAnalysisRequest(BaseModel):
    """Request model for synergy analysis."""
    cannabinoid: str = Field(..., description="Cannabinoid name (CBD, THC, CBG, etc.)")
    terpene_profile: Dict[str, float] = Field(
        ..., 
        description="Terpene concentrations as percentage (e.g., {'myrcene': 1.0, 'limonene': 0.5})"
    )
    min_evidence_tier: int = Field(
        3, 
        ge=1, 
        le=5, 
        description="Minimum evidence tier (1=very low, 5=high)"
    )
    include_all_synergies: bool = Field(
        False, 
        description="Include synergies below evidence threshold"
    )


class ProfileComparisonRequest(BaseModel):
    """Request model for profile comparison."""
    cannabinoid: str = Field(..., description="Cannabinoid to compare for")
    profile_a: Dict[str, float] = Field(..., description="First terpene profile")
    profile_b: Dict[str, float] = Field(..., description="Second terpene profile")
    name_a: str = Field("Profile A", description="Name for first profile")
    name_b: str = Field("Profile B", description="Name for second profile")
    min_evidence_tier: int = Field(3, ge=1, le=5)


class OptimalProfileRequest(BaseModel):
    """Request model for optimal profile recommendation."""
    cannabinoid: str = Field(..., description="Cannabinoid base")
    therapeutic_goal: str = Field(..., description="Therapeutic goal (pain, anxiety, etc.)")
    target_synergy: str = Field(
        ..., 
        description="Target synergy type (enhanced_analgesic, anxiolytic, etc.)"
    )
    max_terpenes: int = Field(5, ge=1, le=10, description="Maximum terpenes to include")
    min_evidence_tier: int = Field(3, ge=1, le=5)


class TerpeneAddRequest(BaseModel):
    """Request model for adding a terpene."""
    name: str = Field(..., min_length=2)
    cas_number: Optional[str] = Field(None)
    molecular_weight: Optional[float] = Field(None)
    effects: Optional[List[str]] = Field(None)
    cannabinoid_synergies: Optional[Dict[str, Dict]] = Field(None)
    sources: Optional[List[str]] = Field(None)


# ============================================================================
# Terpene Database Endpoints
# ============================================================================

@router.get("/")
async def list_terpenes(
    include_synergies: bool = Query(False, description="Include synergy data")
):
    """List all terpenes in the database."""
    analyzer = TerpeneAnalyzer()
    terpene_names = analyzer.terpene_database.list_all()
    
    result = []
    for name in terpene_names:
        terpene = analyzer.terpene_database.get(name)
        if terpene:
            if include_synergies:
                result.append(terpene.to_dict())
            else:
                result.append({
                    "name": terpene.name,
                    "mechanism": terpene.mechanism,
                    "therapeutic_effects": terpene.therapeutic_effects,
                    "boiling_point_c": terpene.boiling_point_c,
                    "aroma_profile": terpene.aroma_profile
                })
    
    return {
        "count": len(result),
        "terpenes": result
    }


@router.get("/therapeutic-targets")
async def get_therapeutic_targets():
    """Get list of supported therapeutic targets for optimization."""
    return {
        "therapeutic_targets": [
            {
                "target": "pain",
                "related_synergies": ["analgesic", "anti-inflammatory"],
                "key_terpenes": ["myrcene", "beta_caryophyllene", "linalool"]
            },
            {
                "target": "anxiety",
                "related_synergies": ["anxiolytic", "sedative"],
                "key_terpenes": ["linalool", "limonene", "myrcene"]
            },
            {
                "target": "inflammation",
                "related_synergies": ["anti-inflammatory", "immunomodulatory"],
                "key_terpenes": ["beta_caryophyllene", "pinene", "limonene"]
            },
            {
                "target": "neuroprotection",
                "related_synergies": ["neuroprotective", "antioxidant"],
                "key_terpenes": ["pinene", "limonene", "beta_caryophyllene"]
            },
            {
                "target": "sleep",
                "related_synergies": ["sedative", "anxiolytic"],
                "key_terpenes": ["myrcene", "linalool", "terpinolene"]
            },
            {
                "target": "focus",
                "related_synergies": ["cognitive_enhancement", "bronchodilation"],
                "key_terpenes": ["pinene", "limonene", "terpinolene"]
            }
        ]
    }


@router.get("/{terpene_name}")
async def get_terpene(terpene_name: str):
    """Get detailed information about a specific terpene."""
    analyzer = TerpeneAnalyzer()
    terpene = analyzer.terpene_database.get(terpene_name)
    
    if not terpene:
        raise HTTPException(
            status_code=404, 
            detail=f"Terpene '{terpene_name}' not found"
        )
    
    return {
        "terpene": terpene.to_dict()
    }


@router.get("/{terpene_name}/synergies/{cannabinoid}")
async def get_terpene_synergies(
    terpene_name: str,
    cannabinoid: str,
    min_evidence_tier: int = Query(1, ge=1, le=5)
):
    """Get synergy data for a terpene-cannabinoid pair."""
    analyzer = TerpeneAnalyzer()
    terpene = analyzer.terpene_database.get(terpene_name)
    
    if not terpene:
        raise HTTPException(
            status_code=404, 
            detail=f"Terpene '{terpene_name}' not found"
        )
    
    cannabinoid_upper = cannabinoid.upper()
    synergies = terpene.cannabinoid_synergy.get(cannabinoid_upper, {})
    
    if not synergies:
        return {
            "terpene": terpene_name,
            "cannabinoid": cannabinoid_upper,
            "synergies": {},
            "message": f"No synergy data for {terpene_name} + {cannabinoid_upper}"
        }
    
    # Filter by evidence tier
    min_tier = EvidenceTier(min_evidence_tier)
    filtered_synergies = {}
    
    # The synergy structure is {enhancement_factor, evidence_tier}
    if synergies.get("evidence_tier", 1) >= min_tier:
        filtered_synergies = synergies
    
    return {
        "terpene": terpene_name,
        "cannabinoid": cannabinoid_upper,
        "synergies": filtered_synergies,
        "evidence_filter": {
            "min_tier": min_evidence_tier,
            "tier_name": min_tier.name
        }
    }


# ============================================================================
# Synergy Analysis Endpoints
# ============================================================================

@router.post("/synergy/analyze")
async def analyze_synergy(request: SynergyAnalysisRequest):
    """Analyze terpene-cannabinoid synergies with evidence gating.
    
    Analyzes the synergistic effects between terpenes and cannabinoids,
    filtering results by evidence tier quality.
    
    Evidence Tiers:
    - 1: Very Low (theoretical only)
    - 2: Low (in vitro studies)
    - 3: Moderate (animal studies)
    - 4: Good (pilot human studies)
    - 5: High (RCTs, clinical trials)
    """
    try:
        analyzer = TerpeneAnalyzer()
        result = analyzer.analyze_synergy(
            cannabinoid=request.cannabinoid,
            terpene_profile=request.terpene_profile,
            min_evidence_tier=request.min_evidence_tier
        )
        
        response = result
        
        # Optionally include filtered-out synergies
        if request.include_all_synergies:
            all_result = analyzer.analyze_synergy(
                cannabinoid=request.cannabinoid,
                terpene_profile=request.terpene_profile,
                min_evidence_tier=1
            )
            response["all_synergies"] = all_result["synergies"]
            response["filtered_count"] = all_result["num_synergies_included"] - result["num_synergies_included"]
        
        return {
            "status": "success",
            "analysis": response
        }
        
    except Exception as e:
        logger.error(f"Synergy analysis error: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/profiles/compare")
async def compare_profiles(request: ProfileComparisonRequest):
    """Compare two terpene profiles for synergistic potential.
    
    Analyzes both profiles and provides a comparison showing
    which is likely to produce stronger synergistic effects.
    """
    try:
        analyzer = TerpeneAnalyzer()
        comparison = analyzer.compare_profiles(
            cannabinoid=request.cannabinoid,
            profile_a=request.profile_a,
            profile_b=request.profile_b,
            min_evidence_tier=request.min_evidence_tier
        )
        
        # Add names to the response
        comparison["profile_a"]["name"] = request.name_a
        comparison["profile_b"]["name"] = request.name_b
        
        return {
            "status": "success",
            "comparison": comparison
        }
        
    except Exception as e:
        logger.error(f"Profile comparison error: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/profiles/optimize")
async def get_optimal_profile(request: OptimalProfileRequest):
    """Get optimal terpene profile for a therapeutic goal.
    
    Recommends terpene combinations that maximize synergistic
    effects for a specific therapeutic target.
    """
    try:
        analyzer = TerpeneAnalyzer()
        optimal = analyzer.get_optimal_profile(
            cannabinoid=request.cannabinoid,
            therapeutic_goal=request.therapeutic_goal,
            min_evidence_tier=request.min_evidence_tier
        )
        
        # Limit to max_terpenes
        if request.max_terpenes and len(optimal.get("recommended_terpenes", [])) > request.max_terpenes:
            optimal["recommended_terpenes"] = optimal["recommended_terpenes"][:request.max_terpenes]
        
        # Add target synergy to response
        optimal["target_synergy"] = request.target_synergy
        
        return {
            "status": "success",
            "recommendation": optimal
        }
        
    except Exception as e:
        logger.error(f"Optimal profile error: {e}")
        raise HTTPException(status_code=500, detail=str(e))


# ============================================================================
# Evidence Information Endpoints
# ============================================================================

@router.get("/evidence/tiers")
async def get_evidence_tiers():
    """Get evidence tier definitions and descriptions."""
    tier_descriptions = {
        1: {
            "name": "Very Low",
            "label": EvidenceTier.TIER_1.name,
            "description": "Theoretical or computational predictions only",
            "example": "Molecular docking simulations"
        },
        2: {
            "name": "Low",
            "label": EvidenceTier.TIER_2.name,
            "description": "In vitro (cell culture) studies",
            "example": "Cell line experiments showing receptor binding"
        },
        3: {
            "name": "Moderate",
            "label": EvidenceTier.TIER_3.name,
            "description": "Animal model studies",
            "example": "Mouse/rat studies showing behavioral effects"
        },
        4: {
            "name": "Good",
            "label": EvidenceTier.TIER_4.name,
            "description": "Pilot human studies or observational data",
            "example": "Small open-label human trials"
        },
        5: {
            "name": "High",
            "label": EvidenceTier.TIER_5.name,
            "description": "Randomized controlled trials (RCTs)",
            "example": "Double-blind placebo-controlled clinical trials"
        }
    }
    
    return {
        "evidence_tiers": tier_descriptions,
        "recommendation": "Use tier 3+ for regulatory-relevant claims",
        "default_filter": 3
    }


@router.get("/evidence/summary")
async def get_evidence_summary():
    """Get summary of evidence coverage in the database."""
    analyzer = TerpeneAnalyzer()
    terpene_names = analyzer.terpene_database.list_all()
    
    tier_counts = {i: 0 for i in range(1, 6)}
    cannabinoid_coverage = {}
    total_synergies = 0
    
    for name in terpene_names:
        terpene = analyzer.terpene_database.get(name)
        if terpene:
            for cannabinoid, synergy_data in terpene.cannabinoid_synergy.items():
                if cannabinoid not in cannabinoid_coverage:
                    cannabinoid_coverage[cannabinoid] = 0
                
                tier = synergy_data.get("evidence_tier", 1)
                if hasattr(tier, 'value'):
                    tier = tier.value
                tier_counts[tier] += 1
                cannabinoid_coverage[cannabinoid] += 1
                total_synergies += 1
    
    return {
        "total_terpenes": len(terpene_names),
        "total_synergies": total_synergies,
        "evidence_tier_distribution": tier_counts,
        "cannabinoid_coverage": cannabinoid_coverage,
        "high_evidence_synergies": tier_counts[4] + tier_counts[5],
        "percentage_moderate_or_better": round(
            (tier_counts[3] + tier_counts[4] + tier_counts[5]) / max(total_synergies, 1) * 100, 1
        )
    }


# ============================================================================
# Strain/Product Analysis Endpoints
# ============================================================================

@router.post("/strain/analyze")
async def analyze_strain(
    strain_name: str,
    cannabinoid_profile: Dict[str, float],
    terpene_profile: Dict[str, float],
    min_evidence_tier: int = Query(3, ge=1, le=5)
):
    """Analyze a cannabis strain's potential synergies.
    
    Provides comprehensive synergy analysis for a complete strain profile.
    """
    try:
        analyzer = TerpeneAnalyzer()
        
        # Analyze synergies for each cannabinoid in the profile
        all_synergies = {}
        total_enhancement = 0.0
        
        for cannabinoid, concentration in cannabinoid_profile.items():
            if concentration > 0:
                result = analyzer.analyze_synergy(
                    cannabinoid=cannabinoid,
                    terpene_profile=terpene_profile,
                    min_evidence_tier=min_evidence_tier
                )
                all_synergies[cannabinoid] = {
                    "concentration": concentration,
                    "synergies": result["synergies"],
                    "total_enhancement": result["total_enhancement_factor"]
                }
                total_enhancement += result["total_enhancement_factor"] * (concentration / 100)
        
        # Determine dominant effects
        effect_scores = {}
        for cannabinoid, data in all_synergies.items():
            for synergy in data["synergies"]:
                effect = synergy.get("mechanism", "unknown")[:20]
                score = synergy.get("enhancement_factor", 1.0)
                if effect not in effect_scores:
                    effect_scores[effect] = 0
                effect_scores[effect] += score
        
        dominant_effects = sorted(
            effect_scores.items(), 
            key=lambda x: x[1], 
            reverse=True
        )[:5]
        
        return {
            "status": "success",
            "strain_name": strain_name,
            "analysis": {
                "cannabinoid_synergies": all_synergies,
                "weighted_enhancement": round(total_enhancement, 3),
                "dominant_effects": [
                    {"effect": e[0], "score": round(e[1], 2)}
                    for e in dominant_effects
                ],
                "evidence_filter": {
                    "min_tier": min_evidence_tier,
                    "tier_name": EvidenceTier(min_evidence_tier).name
                }
            }
        }
        
    except Exception as e:
        logger.error(f"Strain analysis error: {e}")
        raise HTTPException(status_code=500, detail=str(e))


