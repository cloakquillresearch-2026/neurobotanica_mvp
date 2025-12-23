"""
Dimer API Endpoints
FastAPI routes for dimeric cannabinoid prediction and analysis

Supports NeuroBotanica patent claims:
- Dimeric cannabinoid prediction
- Conformer generation for dimers
- Triangulation scoring
- Synergy prediction

Reference: NeuroBotanica MVP Development Plan - Week 4 Task 4.4
"""
from typing import Optional, List, Dict, Any
from datetime import datetime
import logging

from fastapi import APIRouter, HTTPException, Depends, Query
from pydantic import BaseModel, Field
from sqlalchemy.ext.asyncio import AsyncSession

from backend.models.database import get_db
from backend.services.dimer_predictor import (
    DimericPredictor,
    DimerPrediction,
    LinkageType,
    DimerType
)
from backend.services.dimer_conformer_generator import (
    DimericConformerGenerator,
    DimericConformerResult
)
from backend.services.triangulation_scorer import (
    TriangulationScorer,
    TriangulationResult,
    ExperimentalStatus
)

logger = logging.getLogger(__name__)
router = APIRouter(prefix="/api/dimers", tags=["dimers"])


# -----------------------------------------------------------------------------
# Request/Response Models
# -----------------------------------------------------------------------------

class MonomerInput(BaseModel):
    """Input model for a single monomer."""
    name: str
    smiles: str
    
    class Config:
        json_schema_extra = {
            "example": {
                "name": "THC",
                "smiles": "CCCCCc1cc(O)c2C3C=C(C)CCC3C(C)(C)Oc2c1"
            }
        }


class HomodimerRequest(BaseModel):
    """Request for homodimer prediction."""
    monomer: MonomerInput
    linkage_type: Optional[str] = None  # METHYLENE, ETHER, etc.
    
    class Config:
        json_schema_extra = {
            "example": {
                "monomer": {
                    "name": "CBD",
                    "smiles": "CCCCCc1cc(O)c(C2C=C(C)CCC2C(C)=C)c(O)c1"
                },
                "linkage_type": "METHYLENE"
            }
        }


class HeterodimerRequest(BaseModel):
    """Request for heterodimer prediction."""
    monomer_a: MonomerInput
    monomer_b: MonomerInput
    linkage_type: Optional[str] = None
    
    class Config:
        json_schema_extra = {
            "example": {
                "monomer_a": {
                    "name": "THC",
                    "smiles": "CCCCCc1cc(O)c2C3C=C(C)CCC3C(C)(C)Oc2c1"
                },
                "monomer_b": {
                    "name": "CBD",
                    "smiles": "CCCCCc1cc(O)c(C2C=C(C)CCC2C(C)=C)c(O)c1"
                },
                "linkage_type": "METHYLENE"
            }
        }


class BulkPredictionRequest(BaseModel):
    """Request for bulk dimer predictions."""
    monomers: List[MonomerInput]
    include_homodimers: bool = True
    include_heterodimers: bool = True
    top_n: int = Field(default=50, ge=1, le=500)


class ConformerRequest(BaseModel):
    """Request for conformer generation."""
    dimer_smiles: str
    monomer_a_smiles: Optional[str] = None
    monomer_b_smiles: Optional[str] = None
    num_conformers: int = Field(default=20, ge=1, le=100)
    optimize: bool = True


class TriangulationRequest(BaseModel):
    """Request for triangulation scoring."""
    formation_probability: Optional[float] = Field(None, ge=0, le=1)
    structural_similarity: Optional[float] = Field(None, ge=0, le=1)
    ml_therapeutic_score: Optional[float] = Field(None, ge=0, le=1)
    experimental_validation: Optional[str] = None  # confirmed/tentative/predicted/contradicted
    literature_score: Optional[float] = Field(None, ge=0, le=1)
    database_score: Optional[float] = Field(None, ge=0, le=1)


class DimerPredictionResponse(BaseModel):
    """Response model for dimer prediction."""
    dimer_name: str
    dimer_smiles: str
    dimer_type: str
    linkage_type: str
    monomer_a_name: str
    monomer_b_name: str
    formation_probability: float
    synergy_score: float
    novelty_score: float
    structural_validity: float
    predicted_properties: Dict[str, Any]
    therapeutic_potential: Dict[str, float]
    
    class Config:
        json_schema_extra = {
            "example": {
                "dimer_name": "THC-CBD dimer",
                "dimer_smiles": "CCCCCc1cc(O)c2...",
                "dimer_type": "heterodimer",
                "linkage_type": "methylene",
                "monomer_a_name": "THC",
                "monomer_b_name": "CBD",
                "formation_probability": 0.72,
                "synergy_score": 0.85,
                "novelty_score": 0.45,
                "structural_validity": 0.90,
                "predicted_properties": {
                    "molecular_weight": 625.8,
                    "logp": 7.2,
                    "cb1_affinity": 0.75,
                    "cb2_affinity": 0.68
                },
                "therapeutic_potential": {
                    "pain": 0.85,
                    "inflammation": 0.78,
                    "anxiety": 0.70
                }
            }
        }


class ConformerResponse(BaseModel):
    """Response model for conformer generation."""
    smiles: str
    num_conformers: int
    successful_conformers: int
    lowest_energy: float
    conformer_energies: List[float]
    dimer_analysis: Optional[Dict[str, Any]] = None
    best_for_binding: Optional[int] = None
    best_for_stability: Optional[int] = None


class TriangulationResponse(BaseModel):
    """Response model for triangulation scoring."""
    triangulation_score: float
    uncertainty: float
    confidence_interval: Dict[str, float]
    validation_grade: str
    num_sources: int
    sources_used: List[str]
    source_contributions: Dict[str, float]
    experimental_status: str
    source_agreement: float
    consensus_strength: str
    confidence_level_label: str
    recommended_actions: List[str]


# -----------------------------------------------------------------------------
# API Endpoints
# -----------------------------------------------------------------------------

@router.post("/predict/homodimer", response_model=DimerPredictionResponse)
async def predict_homodimer(
    request: HomodimerRequest,
    db: AsyncSession = Depends(get_db)
) -> DimerPredictionResponse:
    """Predict homodimer (A-A dimer) from a single monomer.
    
    Args:
        request: Monomer information and optional linkage type
        
    Returns:
        Dimer prediction with properties and therapeutic potential
    """
    try:
        predictor = DimericPredictor()
        linkage = LinkageType.METHYLENE  # Default
        if request.linkage_type:
            try:
                linkage = LinkageType(request.linkage_type.lower())
            except ValueError:
                pass  # Use default
        # Prepare dimer entry for next-gen model
        dimer_entry = {
            'compound_1': request.monomer.name,
            'compound_2': request.monomer.name,
            'dimer_type': 'homodimer',
            'effect_size': None,
            'confidence': None,
            'regulatory_status': None,
            'benefit_risk_ratio': 1.0,
            'population_group': None,
            'omics_signature': None,
            'pathway': None,
            'patient_stratification': None
        }
        # Call next-gen model
        try:
            nextgen_result = predictor.predict_nextgen_dimer(dimer_entry)
        except Exception as e:
            logger.warning(f"NextGen model prediction failed: {e}")
            nextgen_result = {}
        # Legacy prediction for structure, etc.
        prediction = predictor.predict_homodimer(
            monomer_smiles=request.monomer.smiles,
            monomer_name=request.monomer.name,
            linkage_type=linkage
        )
        if prediction is None:
            raise HTTPException(
                status_code=400,
                detail=f"Could not generate homodimer for {request.monomer.name}"
            )
        # Merge next-gen results into predicted_properties
        predicted_properties = prediction.to_dict().get('predicted_properties', {})
        predicted_properties.update(nextgen_result)
        return DimerPredictionResponse(
            dimer_name=prediction.dimer_name,
            dimer_smiles=prediction.predicted_smiles,
            dimer_type=prediction.dimer_type.value,
            linkage_type=prediction.linkage_type.value,
            monomer_a_name=prediction.parent_1_name,
            monomer_b_name=prediction.parent_2_name,
            formation_probability=prediction.formation_probability,
            synergy_score=nextgen_result.get('synergy_score', prediction.synergy_prediction),
            novelty_score=prediction.novelty_score,
            structural_validity=prediction.structural_validity,
            predicted_properties=predicted_properties,
            therapeutic_potential=prediction.therapeutic_potential
        )
    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Error predicting homodimer: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/predict/heterodimer", response_model=DimerPredictionResponse)
async def predict_heterodimer(
    request: HeterodimerRequest,
    db: AsyncSession = Depends(get_db)
) -> DimerPredictionResponse:
    """Predict heterodimer (A-B dimer) from two monomers.
    
    Args:
        request: Two monomers and optional linkage type
        
    Returns:
        Dimer prediction with properties and therapeutic potential
    """
    try:
        predictor = DimericPredictor()
        linkage = None
        if request.linkage_type:
            try:
                linkage = LinkageType(request.linkage_type.upper())
            except ValueError:
                pass  # Use default
        dimer_entry = {
            'compound_1': request.monomer_a.name,
            'compound_2': request.monomer_b.name,
            'dimer_type': 'heterodimer',
            'effect_size': None,
            'confidence': None,
            'regulatory_status': None,
            'benefit_risk_ratio': 1.0,
            'population_group': None,
            'omics_signature': None,
            'pathway': None,
            'patient_stratification': None
        }
        try:
            nextgen_result = predictor.predict_nextgen_dimer(dimer_entry)
        except Exception as e:
            logger.warning(f"NextGen model prediction failed: {e}")
            nextgen_result = {}
        prediction = predictor.predict_heterodimer(
            monomer1_name=request.monomer_a.name,
            monomer1_smiles=request.monomer_a.smiles,
            monomer2_name=request.monomer_b.name,
            monomer2_smiles=request.monomer_b.smiles,
            linkage_type=linkage
        )
        if prediction is None:
            raise HTTPException(
                status_code=400,
                detail=f"Could not generate heterodimer for {request.monomer_a.name}-{request.monomer_b.name}"
            )
        predicted_properties = prediction.to_dict().get('predicted_properties', {})
        predicted_properties.update(nextgen_result)
        return DimerPredictionResponse(
            dimer_name=prediction.dimer_name,
            dimer_smiles=prediction.predicted_smiles,
            dimer_type=prediction.dimer_type.value,
            linkage_type=prediction.linkage_type.value,
            monomer_a_name=prediction.parent_1_name,
            monomer_b_name=prediction.parent_2_name,
            formation_probability=prediction.formation_probability,
            synergy_score=nextgen_result.get('synergy_score', prediction.synergy_prediction),
            novelty_score=prediction.novelty_score,
            structural_validity=prediction.structural_validity,
            predicted_properties=predicted_properties,
            therapeutic_potential=prediction.therapeutic_potential
        )
    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Error predicting heterodimer: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/predict/bulk")
async def predict_bulk_dimers(
    request: BulkPredictionRequest,
    db: AsyncSession = Depends(get_db)
) -> Dict[str, Any]:
    """Generate bulk dimer predictions from a list of monomers.
    
    Args:
        request: List of monomers and prediction options
        
    Returns:
        Ranked list of dimer predictions
    """
    try:
        predictor = DimericPredictor()
        
        # Build monomer dictionary
        monomers = {
            m.name: m.smiles for m in request.monomers
        }
        
        # Generate all combinations
        all_predictions = predictor.generate_all_combinations(
            monomers,
            include_homodimers=request.include_homodimers,
            include_heterodimers=request.include_heterodimers
        )
        
        # Rank predictions
        ranked = predictor.rank_predictions(all_predictions)[:request.top_n]
        
        # Format response
        results = []
        for pred in ranked:
            results.append({
                "dimer_name": pred.dimer_name,
                "dimer_smiles": pred.dimer_smiles,
                "dimer_type": pred.dimer_type.value,
                "linkage_type": pred.linkage_type.value,
                "monomer_a": pred.monomer_a_name,
                "monomer_b": pred.monomer_b_name,
                "formation_probability": pred.formation_probability,
                "synergy_score": pred.synergy_score,
                "novelty_score": pred.novelty_score,
                "therapeutic_potential": pred.therapeutic_potential
            })
        
        return {
            "total_predictions": len(all_predictions),
            "returned_predictions": len(results),
            "predictions": results
        }
        
    except Exception as e:
        logger.error(f"Error in bulk prediction: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/conformers", response_model=ConformerResponse)
async def generate_dimer_conformers(
    request: ConformerRequest,
    db: AsyncSession = Depends(get_db)
) -> ConformerResponse:
    """Generate 3D conformers for a dimer structure.
    
    Args:
        request: Dimer SMILES and generation parameters
        
    Returns:
        Conformer analysis with dimer-specific geometry
    """
    try:
        generator = DimericConformerGenerator()
        
        result = generator.generate_dimer_conformers(
            dimer_smiles=request.dimer_smiles,
            monomer_a_smiles=request.monomer_a_smiles,
            monomer_b_smiles=request.monomer_b_smiles,
            num_conformers=request.num_conformers,
            optimize=request.optimize
        )
        
        if result is None:
            raise HTTPException(
                status_code=400,
                detail="Could not generate conformers for the provided structure"
            )
        
        # Format dimer analysis
        dimer_analysis = None
        if result.dimer_analysis:
            agg = result.dimer_analysis.get("aggregate", {})
            dimer_analysis = {
                "avg_intermonomer_distance": agg.get("avg_intermonomer_distance"),
                "distance_range": agg.get("distance_range"),
                "avg_planarity": agg.get("avg_planarity"),
                "conformational_flexibility": agg.get("conformational_flexibility"),
                "best_for_binding": result.best_conformer_binding,
                "best_for_stability": result.best_conformer_stability
            }
        
        return ConformerResponse(
            smiles=result.smiles,
            num_conformers=result.num_conformers,
            successful_conformers=result.successful_conformers,
            lowest_energy=result.lowest_energy,
            conformer_energies=result.conformer_energies,
            dimer_analysis=dimer_analysis,
            best_for_binding=result.best_conformer_binding,
            best_for_stability=result.best_conformer_stability
        )
        
    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Error generating conformers: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/triangulate", response_model=TriangulationResponse)
async def triangulate_prediction(
    request: TriangulationRequest,
    db: AsyncSession = Depends(get_db)
) -> TriangulationResponse:
    """Calculate triangulation score for a dimer prediction.
    
    Combines multiple validation sources to calculate confidence:
    - Computational formation probability
    - Structural similarity to known dimers
    - ML-predicted therapeutic potential
    - Experimental validation (if available)
    - Literature evidence (if available)
    
    Args:
        request: Validation scores from different sources
        
    Returns:
        Triangulation result with confidence and recommendations
    """
    try:
        scorer = TriangulationScorer()
        
        result = scorer.calculate_triangulation_score(
            formation_probability=request.formation_probability,
            structural_similarity=request.structural_similarity,
            ml_therapeutic_score=request.ml_therapeutic_score,
            experimental_validation=request.experimental_validation,
            literature_score=request.literature_score,
            database_score=request.database_score
        )
        
        return TriangulationResponse(
            triangulation_score=result.triangulation_score,
            uncertainty=result.uncertainty,
            confidence_interval={
                "lower": result.confidence_interval_lower,
                "upper": result.confidence_interval_upper,
                "level": result.confidence_level
            },
            validation_grade=result.validation_grade,
            num_sources=result.num_sources,
            sources_used=result.sources_used,
            source_contributions=result.source_contributions,
            experimental_status=result.experimental_status.value,
            source_agreement=result.source_agreement,
            consensus_strength=result.consensus_strength,
            confidence_level_label=result.confidence_level_label,
            recommended_actions=result.recommended_actions
        )
        
    except Exception as e:
        logger.error(f"Error calculating triangulation: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/{dimer_id}/triangulation", response_model=TriangulationResponse)
async def get_dimer_triangulation(
    dimer_id: int,
    db: AsyncSession = Depends(get_db)
) -> TriangulationResponse:
    """Get triangulation score for a stored dimer prediction.
    
    Args:
        dimer_id: ID of the dimer prediction in database
        
    Returns:
        Triangulation result with confidence and recommendations
    """
    # TODO: Implement database lookup for stored dimer predictions
    # For now, return a placeholder response
    raise HTTPException(
        status_code=501,
        detail="Database storage for dimer predictions not yet implemented"
    )


@router.get("/linkage-types")
async def list_linkage_types() -> Dict[str, Any]:
    """List available linkage types for dimer prediction."""
    return {
        "linkage_types": [
            {
                "name": "METHYLENE",
                "description": "CH2 bridge between monomers",
                "example": "R-CH2-R'"
            },
            {
                "name": "ETHER",
                "description": "Oxygen bridge between monomers",
                "example": "R-O-R'"
            },
            {
                "name": "ESTER",
                "description": "Ester linkage between monomers",
                "example": "R-COO-R'"
            },
            {
                "name": "DIRECT",
                "description": "Direct carbon-carbon bond",
                "example": "R-R'"
            },
            {
                "name": "CARBON_CHAIN",
                "description": "Alkyl chain bridge",
                "example": "R-(CH2)n-R'"
            }
        ]
    }


@router.get("/statistics")
async def get_dimer_statistics(
    db: AsyncSession = Depends(get_db)
) -> Dict[str, Any]:
    """Get statistics about dimer predictions in the database."""
    # TODO: Implement database statistics
    return {
        "message": "Database statistics not yet implemented",
        "placeholder_stats": {
            "total_predictions": 0,
            "homodimers": 0,
            "heterodimers": 0,
            "avg_formation_probability": 0.0,
            "avg_synergy_score": 0.0,
            "grade_distribution": {
                "A": 0,
                "B": 0,
                "C": 0,
                "D": 0
            }
        }
    }
