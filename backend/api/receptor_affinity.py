"""
Receptor Affinity API Endpoints
Manage receptor binding data with provenance and heterogeneity analysis

Supports NeuroBotanica patent claims:
- Provenance-rich receptor affinity data
- Assay heterogeneity detection
- Quality-weighted aggregation

Reference: NeuroBotanica MVP Development Plan - Week 3 Tasks 3.2-3.3
"""
from typing import Optional, List
from datetime import datetime
from pydantic import BaseModel, Field
from fastapi import APIRouter, HTTPException, status, Depends, Query
from sqlalchemy.orm import Session
from sqlalchemy import func, and_, or_

from ..models.database import get_db
from ..models.receptor_affinity import (
    ReceptorAffinity,
    ReceptorAffinityAggregate,
    ReceptorAffinityCreate,
    ReceptorAffinityResponse,
    AggregatedAffinityResponse,
    AffinityMeasurement,
    AssayDetails,
    SourceAttribution,
    QualityMetrics,
    calculate_confidence_score,
    normalize_affinity_unit,
    detect_assay_heterogeneity
)
from ..models.compound import Cannabinoid
from ..services.assay_analyzer import AssayAnalyzer, HeterogeneitReport

router = APIRouter(prefix="/receptor-affinity", tags=["Receptor Affinity"])

# Initialize analyzer
analyzer = AssayAnalyzer(strict_mode=False)


# ==================== Response Models ====================

class HeterogeneityResponse(BaseModel):
    """Response for heterogeneity analysis."""
    compound_name: str
    receptor: str
    total_measurements: int
    heterogeneity_score: float
    can_aggregate: bool
    overall_severity: str
    issues: List[dict]
    warnings: List[str]
    recommendations: List[str]


class AggregationResponse(BaseModel):
    """Response for weighted aggregation."""
    compound_name: str
    receptor: str
    weighted_mean: float
    geometric_mean: Optional[float]
    median: float
    min_value: float
    max_value: float
    n_measurements: int
    unit: str
    weighting_method: str
    heterogeneity_warning: Optional[str]


class ReceptorProfile(BaseModel):
    """Complete receptor binding profile for a compound."""
    compound_name: str
    receptors: dict
    total_measurements: int
    primary_receptor: Optional[str]
    selectivity_ratio: Optional[float]
    data_quality_summary: dict


# ==================== CRUD Endpoints ====================

@router.post("/", response_model=ReceptorAffinityResponse, status_code=status.HTTP_201_CREATED)
async def create_receptor_affinity(
    data: ReceptorAffinityCreate,
    db: Session = Depends(get_db)
):
    """Create new receptor affinity record with provenance."""
    # Calculate confidence if not provided
    quality = data.quality or QualityMetrics(
        confidence_score=0.5,
        data_quality="unknown"
    )
    
    if quality.confidence_score == 0.5:
        quality.confidence_score = calculate_confidence_score(
            source_type=data.source.type,
            has_pubmed=data.source.pubmed_id is not None,
            has_chembl=data.source.chembl_assay_id is not None,
            data_quality=quality.data_quality
        )
    
    # Create record
    affinity = ReceptorAffinity(
        cannabinoid_id=data.cannabinoid_id,
        compound_name=data.compound_name,
        receptor=data.receptor,
        affinity_value=data.affinity.value,
        affinity_unit=data.affinity.unit,
        affinity_type=data.affinity.type,
        affinity_modifier=data.affinity.modifier,
        affinity_error=data.affinity.error,
        error_type=data.affinity.error_type,
        assay_type=data.assay.type,
        assay_description=data.assay.description,
        target_organism=data.assay.organism,
        target_cell_type=data.assay.cell_type,
        target_tissue=data.assay.tissue,
        radioligand=data.assay.radioligand,
        temperature=data.assay.temperature,
        source_type=data.source.type,
        pubmed_id=data.source.pubmed_id,
        doi=data.source.doi,
        source_database=data.source.database,
        chembl_assay_id=data.source.chembl_assay_id,
        source_url=data.source.url,
        confidence_score=quality.confidence_score,
        data_quality=quality.data_quality,
        curation_level=quality.curation_level,
        assay_heterogeneity_flag=quality.heterogeneity_flag,
        heterogeneity_notes=quality.heterogeneity_notes,
        extraction_timestamp=datetime.utcnow()
    )
    
    db.add(affinity)
    db.commit()
    db.refresh(affinity)
    
    return _format_affinity_response(affinity)


@router.get("/compound/{compound_name}", response_model=List[ReceptorAffinityResponse])
async def get_affinities_by_compound(
    compound_name: str,
    receptor: Optional[str] = None,
    min_confidence: float = Query(0.0, ge=0.0, le=1.0),
    limit: int = Query(100, le=500),
    db: Session = Depends(get_db)
):
    """Get all receptor affinities for a compound."""
    query = db.query(ReceptorAffinity).filter(
        ReceptorAffinity.compound_name.ilike(f"%{compound_name}%")
    )
    
    if receptor:
        query = query.filter(ReceptorAffinity.receptor.ilike(f"%{receptor}%"))
    
    if min_confidence > 0:
        query = query.filter(ReceptorAffinity.confidence_score >= min_confidence)
    
    affinities = query.order_by(
        ReceptorAffinity.confidence_score.desc()
    ).limit(limit).all()
    
    return [_format_affinity_response(a) for a in affinities]


@router.get("/receptor/{receptor}", response_model=List[ReceptorAffinityResponse])
async def get_affinities_by_receptor(
    receptor: str,
    min_confidence: float = Query(0.0, ge=0.0, le=1.0),
    limit: int = Query(100, le=500),
    db: Session = Depends(get_db)
):
    """Get all affinities for a specific receptor."""
    query = db.query(ReceptorAffinity).filter(
        ReceptorAffinity.receptor.ilike(f"%{receptor}%")
    )
    
    if min_confidence > 0:
        query = query.filter(ReceptorAffinity.confidence_score >= min_confidence)
    
    affinities = query.order_by(
        ReceptorAffinity.affinity_value.asc()  # Lowest Ki = highest affinity
    ).limit(limit).all()
    
    return [_format_affinity_response(a) for a in affinities]


# ==================== Heterogeneity Analysis Endpoints ====================

@router.get("/analyze/heterogeneity/{compound_name}/{receptor}", response_model=HeterogeneityResponse)
async def analyze_heterogeneity(
    compound_name: str,
    receptor: str,
    strict_mode: bool = Query(False, description="Use stricter heterogeneity thresholds"),
    db: Session = Depends(get_db)
):
    """Analyze heterogeneity in receptor affinity data for a compound/receptor pair."""
    # Get all measurements
    affinities = db.query(ReceptorAffinity).filter(
        and_(
            ReceptorAffinity.compound_name.ilike(f"%{compound_name}%"),
            ReceptorAffinity.receptor.ilike(f"%{receptor}%")
        )
    ).all()
    
    if not affinities:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"No affinity data found for {compound_name} at {receptor}"
        )
    
    # Convert to analyzer format
    measurements = [
        {
            "id": a.id,
            "value": a.affinity_value,
            "unit": a.affinity_unit,
            "assay_type": a.assay_type,
            "source_type": a.source_type or "unknown",
            "confidence": a.confidence_score or 0.5,
            "organism": a.target_organism,
            "target_organism": a.target_organism
        }
        for a in affinities
    ]
    
    # Run analysis
    local_analyzer = AssayAnalyzer(strict_mode=strict_mode)
    report = local_analyzer.analyze(measurements, compound_name, receptor)
    
    return HeterogeneityResponse(
        compound_name=report.compound_name,
        receptor=report.receptor,
        total_measurements=report.total_measurements,
        heterogeneity_score=report.heterogeneity_score,
        can_aggregate=report.can_aggregate,
        overall_severity=report.overall_severity.value,
        issues=[
            {
                "type": i.issue_type,
                "severity": i.severity.value,
                "description": i.description,
                "affected_measurements": i.affected_measurements,
                "remediation": i.remediation
            }
            for i in report.issues
        ],
        warnings=report.aggregation_warnings,
        recommendations=report.recommended_actions
    )


@router.get("/aggregate/{compound_name}/{receptor}", response_model=AggregationResponse)
async def aggregate_affinity(
    compound_name: str,
    receptor: str,
    weighting: str = Query("confidence", description="Weighting method: confidence, source, equal"),
    check_heterogeneity: bool = Query(True, description="Check heterogeneity before aggregation"),
    db: Session = Depends(get_db)
):
    """Aggregate receptor affinity with quality weighting."""
    # Get measurements
    affinities = db.query(ReceptorAffinity).filter(
        and_(
            ReceptorAffinity.compound_name.ilike(f"%{compound_name}%"),
            ReceptorAffinity.receptor.ilike(f"%{receptor}%")
        )
    ).all()
    
    if not affinities:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"No affinity data found for {compound_name} at {receptor}"
        )
    
    measurements = [
        {
            "id": a.id,
            "value": a.affinity_value,
            "unit": a.affinity_unit,
            "assay_type": a.assay_type,
            "source_type": a.source_type or "unknown",
            "confidence": a.confidence_score or 0.5,
            "organism": a.target_organism
        }
        for a in affinities
    ]
    
    # Check heterogeneity if requested
    warning = None
    if check_heterogeneity:
        report = analyzer.analyze(measurements, compound_name, receptor)
        if not report.can_aggregate:
            warning = f"High heterogeneity detected (score: {report.heterogeneity_score:.2f}). Results may be unreliable."
    
    # Perform aggregation
    result = analyzer.aggregate_with_weights(measurements, weighting)
    
    if "error" in result:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=result["error"]
        )
    
    return AggregationResponse(
        compound_name=compound_name,
        receptor=receptor,
        weighted_mean=result["weighted_mean"],
        geometric_mean=result.get("geometric_mean"),
        median=result["median"],
        min_value=result["min"],
        max_value=result["max"],
        n_measurements=result["n_measurements"],
        unit=result["unit"],
        weighting_method=result["weighting_method"],
        heterogeneity_warning=warning
    )


# ==================== Profile Endpoints ====================

@router.get("/profile/{compound_name}", response_model=ReceptorProfile)
async def get_receptor_profile(
    compound_name: str,
    db: Session = Depends(get_db)
):
    """Get complete receptor binding profile for a compound."""
    # Get all affinities
    affinities = db.query(ReceptorAffinity).filter(
        ReceptorAffinity.compound_name.ilike(f"%{compound_name}%")
    ).all()
    
    if not affinities:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"No receptor affinity data found for {compound_name}"
        )
    
    # Group by receptor
    receptors = {}
    for a in affinities:
        receptor = a.receptor
        if receptor not in receptors:
            receptors[receptor] = {
                "measurements": [],
                "best_ki": None,
                "sources": set()
            }
        
        receptors[receptor]["measurements"].append(a.affinity_value)
        receptors[receptor]["sources"].add(a.source_type or "unknown")
        
        # Track best (lowest) Ki
        if receptors[receptor]["best_ki"] is None or a.affinity_value < receptors[receptor]["best_ki"]:
            receptors[receptor]["best_ki"] = a.affinity_value
    
    # Calculate stats per receptor
    for receptor, data in receptors.items():
        measurements = data["measurements"]
        data["n_measurements"] = len(measurements)
        data["mean_ki"] = sum(measurements) / len(measurements) if measurements else None
        data["sources"] = list(data["sources"])
        del data["measurements"]  # Don't return raw measurements
    
    # Determine primary receptor and selectivity
    primary_receptor = None
    selectivity_ratio = None
    
    if "CB1" in receptors and "CB2" in receptors:
        cb1_ki = receptors["CB1"]["best_ki"]
        cb2_ki = receptors["CB2"]["best_ki"]
        
        if cb1_ki and cb2_ki and cb1_ki > 0 and cb2_ki > 0:
            selectivity_ratio = cb2_ki / cb1_ki
            primary_receptor = "CB1" if selectivity_ratio > 1 else "CB2"
    
    # Data quality summary
    total = len(affinities)
    high_confidence = sum(1 for a in affinities if (a.confidence_score or 0) >= 0.7)
    
    return ReceptorProfile(
        compound_name=compound_name,
        receptors=receptors,
        total_measurements=total,
        primary_receptor=primary_receptor,
        selectivity_ratio=round(selectivity_ratio, 2) if selectivity_ratio else None,
        data_quality_summary={
            "total_measurements": total,
            "high_confidence": high_confidence,
            "high_confidence_percentage": round(high_confidence / total * 100, 1) if total > 0 else 0,
            "unique_receptors": len(receptors)
        }
    )


# ==================== Statistics Endpoints ====================

@router.get("/stats/overview")
async def get_affinity_stats(
    db: Session = Depends(get_db)
):
    """Get overview statistics for receptor affinity database."""
    total = db.query(ReceptorAffinity).count()
    
    # By receptor
    receptor_counts = db.query(
        ReceptorAffinity.receptor,
        func.count(ReceptorAffinity.id)
    ).group_by(ReceptorAffinity.receptor).all()
    
    # By source type
    source_counts = db.query(
        ReceptorAffinity.source_type,
        func.count(ReceptorAffinity.id)
    ).group_by(ReceptorAffinity.source_type).all()
    
    # By assay type
    assay_counts = db.query(
        ReceptorAffinity.assay_type,
        func.count(ReceptorAffinity.id)
    ).group_by(ReceptorAffinity.assay_type).all()
    
    # Confidence distribution
    avg_confidence = db.query(
        func.avg(ReceptorAffinity.confidence_score)
    ).scalar() or 0.5
    
    high_conf = db.query(ReceptorAffinity).filter(
        ReceptorAffinity.confidence_score >= 0.7
    ).count()
    
    # Unique compounds
    unique_compounds = db.query(
        func.count(func.distinct(ReceptorAffinity.compound_name))
    ).scalar() or 0
    
    return {
        "total_measurements": total,
        "unique_compounds": unique_compounds,
        "receptors": {r[0]: r[1] for r in receptor_counts if r[0]},
        "sources": {s[0] or "unknown": s[1] for s in source_counts},
        "assay_types": {a[0] or "unknown": a[1] for a in assay_counts},
        "average_confidence": round(avg_confidence, 3),
        "high_confidence_count": high_conf,
        "high_confidence_percentage": round(high_conf / total * 100, 1) if total > 0 else 0
    }


# ==================== Helper Functions ====================

def _format_affinity_response(affinity: ReceptorAffinity) -> ReceptorAffinityResponse:
    """Format database record to response model."""
    return ReceptorAffinityResponse(
        id=affinity.id,
        compound_name=affinity.compound_name,
        receptor=affinity.receptor,
        affinity=AffinityMeasurement(
            value=affinity.affinity_value,
            unit=affinity.affinity_unit,
            type=affinity.affinity_type or "Ki",
            modifier=affinity.affinity_modifier,
            error=affinity.affinity_error,
            error_type=affinity.error_type
        ),
        assay=AssayDetails(
            type=affinity.assay_type,
            description=affinity.assay_description,
            organism=affinity.target_organism or "Homo sapiens",
            cell_type=affinity.target_cell_type,
            tissue=affinity.target_tissue,
            radioligand=affinity.radioligand,
            temperature=affinity.temperature
        ),
        source=SourceAttribution(
            type=affinity.source_type or "unknown",
            pubmed_id=affinity.pubmed_id,
            doi=affinity.doi,
            database=affinity.source_database,
            chembl_assay_id=affinity.chembl_assay_id,
            url=affinity.source_url
        ),
        quality=QualityMetrics(
            confidence_score=affinity.confidence_score or 0.5,
            data_quality=affinity.data_quality or "unknown",
            curation_level=affinity.curation_level,
            heterogeneity_flag=affinity.assay_heterogeneity_flag or False,
            heterogeneity_notes=affinity.heterogeneity_notes
        ),
        created_at=affinity.created_at or datetime.utcnow()
    )
