"""
Cannabinoid Compounds API Routes
Access to 63-compound dataset with molecular descriptors
"""
from fastapi import APIRouter, Depends, HTTPException, Query
from sqlalchemy.orm import Session
from sqlalchemy import func
from typing import List, Optional
from pydantic import BaseModel

from backend.models.database import get_db
from backend.models.compound import Cannabinoid, DimericPrediction


router = APIRouter()


class CompoundResponse(BaseModel):
    """Compound response model."""
    name: str
    abbreviation: Optional[str]
    compound_class: Optional[str]
    molecular_weight: Optional[float]
    is_natural: bool
    fda_approved: bool
    
    class Config:
        from_attributes = True


class CompoundDetail(BaseModel):
    """Detailed compound response."""
    name: str
    abbreviation: Optional[str]
    smiles: str
    compound_class: Optional[str]
    molecular_weight: Optional[float]
    logp: Optional[float]
    tpsa: Optional[float]
    is_natural: bool
    is_dimeric: bool
    fda_approved: bool
    fda_drug_name: Optional[str]
    schedule_classification: Optional[str]
    mechanism_of_action: Optional[str]
    therapeutic_categories: Optional[list]
    cb1_affinity_ki: Optional[float]
    cb2_affinity_ki: Optional[float]
    
    class Config:
        from_attributes = True


@router.get("/", response_model=List[CompoundResponse])
async def list_compounds(
    compound_class: Optional[str] = Query(None, description="Filter by class (Phytocannabinoid, Synthetic, etc.)"),
    fda_approved: Optional[bool] = Query(None, description="Filter by FDA approval status"),
    is_dimeric: Optional[bool] = Query(None, description="Filter dimeric compounds"),
    has_conformers: Optional[bool] = Query(None, description="Filter by 3D conformer availability"),
    limit: int = Query(50, le=100),
    offset: int = Query(0),
    db: Session = Depends(get_db)
):
    """List cannabinoid compounds with optional filters.
    
    Returns 63 compounds with molecular descriptors.
    """
    query = db.query(Cannabinoid)
    
    if compound_class:
        query = query.filter(Cannabinoid.compound_class == compound_class)
    if fda_approved is not None:
        query = query.filter(Cannabinoid.fda_approved == fda_approved)
    if is_dimeric is not None:
        query = query.filter(Cannabinoid.is_dimeric == is_dimeric)
    if has_conformers is not None:
        query = query.filter(Cannabinoid.has_conformers == has_conformers)
    
    compounds = query.offset(offset).limit(limit).all()
    return compounds


@router.get("/fda-approved")
async def list_fda_approved(db: Session = Depends(get_db)):
    """Get FDA-approved cannabinoid drugs.
    
    Supports Schedule III documentation with approved drug precedents.
    """
    approved = db.query(Cannabinoid).filter(
        Cannabinoid.fda_approved == True
    ).all()
    
    return {
        "fda_approved_compounds": [
            {
                "name": c.name,
                "drug_name": c.fda_drug_name,
                "indications": c.fda_indications,
                "schedule": c.schedule_classification,
                "mechanism": c.mechanism_of_action
            }
            for c in approved
        ],
        "regulatory_precedent": "These approvals establish pathway for cannabis Schedule III transition"
    }


@router.get("/receptor-affinities")
async def list_receptor_affinities(
    receptor: str = Query("CB1", description="Receptor type (CB1, CB2, TRPV1, GPR55)"),
    min_affinity: Optional[float] = Query(None, description="Minimum affinity (Ki in nM)"),
    max_affinity: Optional[float] = Query(None, description="Maximum affinity (Ki in nM)"),
    db: Session = Depends(get_db)
):
    """Get compounds by receptor affinity.
    
    Supports pharmacology data packages for FDA submissions.
    """
    query = db.query(Cannabinoid)
    
    if receptor.upper() == "CB1":
        query = query.filter(Cannabinoid.cb1_affinity_ki.isnot(None))
        if min_affinity:
            query = query.filter(Cannabinoid.cb1_affinity_ki >= min_affinity)
        if max_affinity:
            query = query.filter(Cannabinoid.cb1_affinity_ki <= max_affinity)
        compounds = query.order_by(Cannabinoid.cb1_affinity_ki).all()
        
    elif receptor.upper() == "CB2":
        query = query.filter(Cannabinoid.cb2_affinity_ki.isnot(None))
        if min_affinity:
            query = query.filter(Cannabinoid.cb2_affinity_ki >= min_affinity)
        if max_affinity:
            query = query.filter(Cannabinoid.cb2_affinity_ki <= max_affinity)
        compounds = query.order_by(Cannabinoid.cb2_affinity_ki).all()
    else:
        compounds = query.all()
    
    return {
        "receptor": receptor,
        "compounds": [
            {
                "name": c.name,
                "cb1_ki": c.cb1_affinity_ki,
                "cb2_ki": c.cb2_affinity_ki,
                "selectivity_ratio": c.cb2_affinity_ki / c.cb1_affinity_ki if c.cb1_affinity_ki and c.cb2_affinity_ki else None
            }
            for c in compounds
        ]
    }


@router.get("/dimeric-predictions")
async def list_dimeric_predictions(
    min_score: Optional[float] = Query(None, description="Minimum triangulation score"),
    limit: int = Query(20, le=50),
    db: Session = Depends(get_db)
):
    """Get dimeric cannabinoid predictions.
    
    Novel compound predictions from triangulation scoring.
    """
    query = db.query(DimericPrediction)
    
    if min_score:
        query = query.filter(DimericPrediction.triangulation_score >= min_score)
    
    predictions = query.order_by(
        DimericPrediction.triangulation_score.desc()
    ).limit(limit).all()
    
    return {
        "total_predictions": len(predictions),
        "predictions": [
            {
                "dimer_name": p.dimer_name,
                "parents": [p.parent_1_name, p.parent_2_name],
                "triangulation_score": p.triangulation_score,
                "synergy_prediction": p.synergy_prediction,
                "novelty_score": p.novelty_score,
                "therapeutic_potential": p.therapeutic_potential,
                "fto_status": p.fto_status
            }
            for p in predictions
        ]
    }


@router.get("/{compound_name}", response_model=CompoundDetail)
async def get_compound(compound_name: str, db: Session = Depends(get_db)):
    """Get detailed compound information by name."""
    compound = db.query(Cannabinoid).filter(
        Cannabinoid.name.ilike(f"%{compound_name}%")
    ).first()
    
    if not compound:
        raise HTTPException(status_code=404, detail=f"Compound {compound_name} not found")
    
    return compound


@router.get("/{compound_name}/cmc")
async def get_compound_cmc(compound_name: str, db: Session = Depends(get_db)):
    """Get CMC (Chemistry Manufacturing Controls) profile.
    
    Patent Claim 1(j): CMC templates for FDA pharmaceutical approval.
    """
    compound = db.query(Cannabinoid).filter(
        Cannabinoid.name.ilike(f"%{compound_name}%")
    ).first()
    
    if not compound:
        raise HTTPException(status_code=404, detail=f"Compound {compound_name} not found")
    
    return compound.to_cmc_profile()


@router.get("/{compound_name}/pharmacology")
async def get_compound_pharmacology(compound_name: str, db: Session = Depends(get_db)):
    """Get pharmacology profile for FDA submission.
    
    Patent Claim 1(j): pharmacology data packages demonstrating mechanism of action.
    """
    compound = db.query(Cannabinoid).filter(
        Cannabinoid.name.ilike(f"%{compound_name}%")
    ).first()
    
    if not compound:
        raise HTTPException(status_code=404, detail=f"Compound {compound_name} not found")
    
    return compound.to_pharmacology_profile()
