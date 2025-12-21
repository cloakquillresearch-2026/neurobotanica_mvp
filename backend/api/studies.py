"""
Clinical Studies API Routes
Access to 320-study evidence database
"""
from fastapi import APIRouter, Depends, HTTPException, Query
from sqlalchemy.orm import Session
from sqlalchemy import func, or_
from typing import List, Optional
from pydantic import BaseModel

from backend.models.database import get_db
from backend.models.study import ClinicalStudy


router = APIRouter()


class StudyResponse(BaseModel):
    """Study response model."""
    study_id: str
    condition: str
    study_type: str
    study_title: str
    cannabinoid: Optional[str]
    year: Optional[int]
    effect_size: Optional[str]
    regulatory_relevance: Optional[str]
    
    class Config:
        from_attributes = True


class StudyDetail(BaseModel):
    """Detailed study response."""
    study_id: str
    condition: str
    study_type: str
    study_title: str
    citation: Optional[str]
    authors: Optional[str]
    year: Optional[int]
    journal: Optional[str]
    sample_size: Optional[str]
    cannabinoid: Optional[str]
    dosage: Optional[str]
    duration: Optional[str]
    primary_measure: Optional[str]
    results_summary: Optional[str]
    effect_size: Optional[str]
    adverse_events: Optional[str]
    regulatory_relevance: Optional[str]
    evidence_grade: Optional[str]
    mechanism_of_action: Optional[str]
    
    class Config:
        from_attributes = True


@router.get("/", response_model=List[StudyResponse])
async def list_studies(
    condition: Optional[str] = Query(None, description="Filter by condition"),
    cannabinoid: Optional[str] = Query(None, description="Filter by cannabinoid"),
    study_type: Optional[str] = Query(None, description="Filter by study type (RCT, Observational, etc.)"),
    year_min: Optional[int] = Query(None, description="Minimum publication year"),
    year_max: Optional[int] = Query(None, description="Maximum publication year"),
    fda_relevant: Optional[bool] = Query(None, description="Filter FDA-relevant studies only"),
    limit: int = Query(50, le=100),
    offset: int = Query(0),
    db: Session = Depends(get_db)
):
    """List clinical studies with optional filters.
    
    Returns paginated list of 320 studies across 16 conditions.
    """
    query = db.query(ClinicalStudy)
    
    if condition:
        query = query.filter(ClinicalStudy.condition.ilike(f"%{condition}%"))
    if cannabinoid:
        query = query.filter(ClinicalStudy.cannabinoid.ilike(f"%{cannabinoid}%"))
    if study_type:
        query = query.filter(ClinicalStudy.study_type == study_type)
    if year_min:
        query = query.filter(ClinicalStudy.year >= year_min)
    if year_max:
        query = query.filter(ClinicalStudy.year <= year_max)
    if fda_relevant:
        query = query.filter(
            or_(
                ClinicalStudy.is_pivotal_trial == True,
                ClinicalStudy.fda_drug_reference.isnot(None)
            )
        )
    
    studies = query.offset(offset).limit(limit).all()
    return studies


@router.get("/conditions")
async def list_conditions(db: Session = Depends(get_db)):
    """Get list of all conditions with study counts."""
    conditions = db.query(
        ClinicalStudy.condition,
        func.count(ClinicalStudy.id).label("study_count")
    ).group_by(ClinicalStudy.condition).all()
    
    return {
        "total_conditions": len(conditions),
        "conditions": [
            {"condition": c[0], "study_count": c[1]}
            for c in conditions
        ]
    }


@router.get("/cannabinoids")
async def list_study_cannabinoids(db: Session = Depends(get_db)):
    """Get list of cannabinoids studied with counts."""
    cannabinoids = db.query(
        ClinicalStudy.cannabinoid,
        func.count(ClinicalStudy.id).label("study_count")
    ).filter(
        ClinicalStudy.cannabinoid.isnot(None)
    ).group_by(ClinicalStudy.cannabinoid).all()
    
    return {
        "total_cannabinoids": len(cannabinoids),
        "cannabinoids": [
            {"cannabinoid": c[0], "study_count": c[1]}
            for c in cannabinoids
        ]
    }


@router.get("/fda-drugs")
async def list_fda_drug_studies(db: Session = Depends(get_db)):
    """Get studies involving FDA-approved cannabinoid drugs.
    
    Supports Patent Claim 1(j): FDA Schedule III compliance documentation.
    """
    fda_drugs = ["Epidiolex", "Marinol", "Cesamet", "Sativex", "Dronabinol", "Nabilone", "CBD", "Nabiximols"]
    
    results = {}
    for drug in fda_drugs:
        studies = db.query(ClinicalStudy).filter(
            or_(
                ClinicalStudy.fda_drug_reference.ilike(f"%{drug}%"),
                ClinicalStudy.cannabinoid.ilike(f"%{drug}%")
            )
        ).all()
        
        if studies:
            results[drug] = {
                "study_count": len(studies),
                "conditions": list(set(s.condition for s in studies)),
                "pivotal_trials": sum(1 for s in studies if s.is_pivotal_trial)
            }
    
    return {
        "fda_approved_drugs": list(results.keys()),
        "drug_evidence": results,
        "regulatory_note": "Evidence supports FDA Schedule III reclassification"
    }


@router.get("/{study_id}", response_model=StudyDetail)
async def get_study(study_id: str, db: Session = Depends(get_db)):
    """Get detailed study information by study ID."""
    study = db.query(ClinicalStudy).filter(
        ClinicalStudy.study_id == study_id
    ).first()
    
    if not study:
        raise HTTPException(status_code=404, detail=f"Study {study_id} not found")
    
    return study


@router.get("/{study_id}/pharmacology")
async def get_study_pharmacology(study_id: str, db: Session = Depends(get_db)):
    """Get pharmacology data package for a study.
    
    Patent Claim 1(j): pharmacology data packages demonstrating mechanism of action.
    """
    study = db.query(ClinicalStudy).filter(
        ClinicalStudy.study_id == study_id
    ).first()
    
    if not study:
        raise HTTPException(status_code=404, detail=f"Study {study_id} not found")
    
    return study.to_pharmacology_package()


@router.get("/{study_id}/efficacy")
async def get_study_efficacy(study_id: str, db: Session = Depends(get_db)):
    """Get comparative efficacy data for a study.
    
    Patent Claim 1(j): comparative efficacy analysis.
    """
    study = db.query(ClinicalStudy).filter(
        ClinicalStudy.study_id == study_id
    ).first()
    
    if not study:
        raise HTTPException(status_code=404, detail=f"Study {study_id} not found")
    
    return study.to_efficacy_comparison()
