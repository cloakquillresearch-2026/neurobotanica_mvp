"""
Clinical Evidence API Endpoints
Query and aggregate clinical evidence from study database

Supports NeuroBotanica patent claims:
- Evidence-based therapeutic recommendations
- Confidence-weighted evidence aggregation
- Citation and provenance tracking

Reference: NeuroBotanica MVP Development Plan - Week 3 Task 3.1
"""
from typing import Optional, List, Dict, Any
from datetime import datetime
from pydantic import BaseModel, Field
from fastapi import APIRouter, HTTPException, status, Depends, Query
from sqlalchemy.orm import Session
from sqlalchemy import func, and_, or_

from ..models.database import get_db
from ..models.study import ClinicalStudy
from ..models.compound import Cannabinoid

router = APIRouter(prefix="/evidence", tags=["Clinical Evidence"])


# ==================== Request/Response Models ====================

class EvidenceSummary(BaseModel):
    """Summary of evidence for a condition or compound."""
    total_studies: int
    favorable_outcomes: int
    neutral_outcomes: int
    negative_outcomes: int
    average_confidence: float
    evidence_strength: str  # "strong", "moderate", "limited", "insufficient"
    conditions: List[str]
    cannabinoids: List[str]


class StudyEvidence(BaseModel):
    """Individual study evidence record."""
    study_id: int
    title: str
    doi: Optional[str]
    pubmed_id: Optional[str]
    condition: str
    study_design: str
    sample_size: Optional[int]
    cannabinoids_studied: List[str]
    outcome_direction: str  # "favorable", "neutral", "negative"
    outcome_summary: Optional[str]
    confidence_weight: float
    year: Optional[int]


class AggregatedEvidence(BaseModel):
    """Aggregated evidence with confidence weighting."""
    condition: str
    total_studies: int
    weighted_efficacy_score: float
    confidence_interval: Dict[str, float]
    evidence_grade: str  # A, B, C, D based on strength
    studies: List[StudyEvidence]
    cannabinoid_breakdown: Dict[str, Dict[str, Any]]
    recommendation: str


class CannabinoidEvidence(BaseModel):
    """Evidence summary for a specific cannabinoid."""
    cannabinoid_name: str
    cannabinoid_id: int
    total_studies: int
    conditions_studied: List[str]
    efficacy_by_condition: Dict[str, Dict[str, Any]]
    overall_evidence_grade: str


# ==================== Helper Functions ====================

def calculate_evidence_strength(
    total_studies: int,
    favorable_ratio: float,
    average_confidence: float
) -> str:
    """Calculate evidence strength classification.
    
    Based on:
    - Number of studies
    - Favorable outcome ratio
    - Average confidence weighting
    
    Returns: "strong", "moderate", "limited", "insufficient"
    """
    if total_studies == 0:
        return "insufficient"
    
    # Weight by study count and quality
    score = (
        (min(total_studies, 20) / 20) * 0.3 +  # Max 20 studies for full points
        favorable_ratio * 0.4 +
        average_confidence * 0.3
    )
    
    if score >= 0.7 and total_studies >= 5:
        return "strong"
    elif score >= 0.5 and total_studies >= 3:
        return "moderate"
    elif score >= 0.3 or total_studies >= 1:
        return "limited"
    else:
        return "insufficient"


def calculate_evidence_grade(
    weighted_score: float,
    num_studies: int,
    has_rct: bool = False
) -> str:
    """Calculate evidence grade (A-D).
    
    A: Strong evidence from multiple high-quality studies
    B: Moderate evidence from good-quality studies
    C: Limited evidence from lower-quality studies
    D: Insufficient or conflicting evidence
    """
    if weighted_score >= 0.75 and num_studies >= 5 and has_rct:
        return "A"
    elif weighted_score >= 0.6 and num_studies >= 3:
        return "B"
    elif weighted_score >= 0.4 or num_studies >= 2:
        return "C"
    else:
        return "D"


def parse_outcome_direction(outcome_text: str) -> str:
    """Parse outcome text to determine direction.
    
    Returns: "favorable", "neutral", "negative"
    """
    if not outcome_text:
        return "neutral"
    
    outcome_lower = outcome_text.lower()
    
    # Negative indicators (check first to avoid false positives like "ineffective")
    negative_terms = [
        "no effect", "no significant", "ineffective", "worsened",
        "adverse", "failed", "negative", "no benefit", "no improvement",
        "not effective", "did not improve", "no difference"
    ]
    
    for term in negative_terms:
        if term in outcome_lower:
            return "negative"
    
    # Favorable indicators
    favorable_terms = [
        "effective", "significant improvement", "reduced", "decreased",
        "beneficial", "positive", "improved", "ameliorated", "efficacy",
        "successful", "promising", "well-tolerated"
    ]
    
    for term in favorable_terms:
        if term in outcome_lower:
            return "favorable"
    
    return "neutral"


def generate_recommendation(
    evidence_grade: str,
    condition: str,
    top_cannabinoid: Optional[str],
    weighted_score: float
) -> str:
    """Generate clinical recommendation based on evidence."""
    if evidence_grade == "A":
        return f"Strong evidence supports cannabinoid therapy for {condition}. {top_cannabinoid or 'Cannabis-based treatments'} shows consistent efficacy across multiple high-quality studies."
    elif evidence_grade == "B":
        return f"Moderate evidence suggests potential benefit of cannabinoid therapy for {condition}. Consider {top_cannabinoid or 'cannabinoid treatment'} with appropriate monitoring."
    elif evidence_grade == "C":
        return f"Limited evidence available for cannabinoid therapy in {condition}. Further research needed. May consider as adjunct therapy with careful evaluation."
    else:
        return f"Insufficient evidence to recommend cannabinoid therapy for {condition}. Clinical trials may be warranted."


# ==================== Evidence Query Endpoints ====================

@router.get("/cannabinoid/{compound_name}", response_model=CannabinoidEvidence)
async def get_evidence_by_cannabinoid(
    compound_name: str,
    db: Session = Depends(get_db)
):
    """Get all clinical evidence for a specific cannabinoid.
    
    Returns studies, conditions studied, and efficacy breakdown.
    """
    # Find cannabinoid
    cannabinoid = db.query(Cannabinoid).filter(
        or_(
            Cannabinoid.name.ilike(f"%{compound_name}%"),
            Cannabinoid.abbreviation.ilike(f"%{compound_name}%")
        )
    ).first()
    
    if not cannabinoid:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"Cannabinoid '{compound_name}' not found"
        )
    
    # Find studies mentioning this cannabinoid
    studies = db.query(ClinicalStudy).filter(
        ClinicalStudy.cannabinoids_studied.contains([cannabinoid.name])
    ).all()
    
    if not studies:
        # Try searching in outcome text
        studies = db.query(ClinicalStudy).filter(
            or_(
                ClinicalStudy.cannabinoids_studied.cast(str).ilike(f"%{compound_name}%"),
                ClinicalStudy.outcome_summary.ilike(f"%{compound_name}%")
            )
        ).all()
    
    # Aggregate by condition
    conditions = set()
    efficacy_by_condition = {}
    
    for study in studies:
        condition = study.condition
        conditions.add(condition)
        
        if condition not in efficacy_by_condition:
            efficacy_by_condition[condition] = {
                "study_count": 0,
                "favorable": 0,
                "neutral": 0,
                "negative": 0,
                "total_confidence": 0.0
            }
        
        efficacy_by_condition[condition]["study_count"] += 1
        efficacy_by_condition[condition]["total_confidence"] += study.confidence_weight or 0.5
        
        direction = parse_outcome_direction(study.outcome_summary)
        efficacy_by_condition[condition][direction] += 1
    
    # Calculate scores per condition
    for condition, data in efficacy_by_condition.items():
        if data["study_count"] > 0:
            data["favorable_ratio"] = data["favorable"] / data["study_count"]
            data["avg_confidence"] = data["total_confidence"] / data["study_count"]
            data["evidence_grade"] = calculate_evidence_grade(
                data["favorable_ratio"] * data["avg_confidence"],
                data["study_count"]
            )
    
    # Calculate overall grade
    if len(studies) >= 5:
        overall_favorable = sum(d["favorable"] for d in efficacy_by_condition.values())
        overall_ratio = overall_favorable / len(studies) if studies else 0
        overall_grade = calculate_evidence_grade(overall_ratio, len(studies))
    else:
        overall_grade = "D"
    
    return CannabinoidEvidence(
        cannabinoid_name=cannabinoid.name,
        cannabinoid_id=cannabinoid.id,
        total_studies=len(studies),
        conditions_studied=list(conditions),
        efficacy_by_condition=efficacy_by_condition,
        overall_evidence_grade=overall_grade
    )


@router.get("/condition/{condition}", response_model=EvidenceSummary)
async def get_evidence_by_condition(
    condition: str,
    min_confidence: float = Query(0.0, ge=0.0, le=1.0),
    db: Session = Depends(get_db)
):
    """Get clinical evidence summary for a specific condition.
    
    Returns study counts, outcome distribution, and evidence strength.
    """
    # Search for studies by condition
    studies = db.query(ClinicalStudy).filter(
        ClinicalStudy.condition.ilike(f"%{condition}%")
    ).all()
    
    if not studies:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"No studies found for condition '{condition}'"
        )
    
    # Filter by confidence if specified
    if min_confidence > 0:
        studies = [s for s in studies if (s.confidence_weight or 0.5) >= min_confidence]
    
    # Count outcomes
    favorable = 0
    neutral = 0
    negative = 0
    total_confidence = 0.0
    cannabinoids = set()
    
    for study in studies:
        direction = parse_outcome_direction(study.results_summary)
        if direction == "favorable":
            favorable += 1
        elif direction == "negative":
            negative += 1
        else:
            neutral += 1
        
        total_confidence += study.confidence_weight or 0.5
        
        if study.cannabinoid:
            cannabinoids.add(study.cannabinoid)
    
    avg_confidence = total_confidence / len(studies) if studies else 0
    favorable_ratio = favorable / len(studies) if studies else 0
    
    evidence_strength = calculate_evidence_strength(
        len(studies),
        favorable_ratio,
        avg_confidence
    )
    
    return EvidenceSummary(
        total_studies=len(studies),
        favorable_outcomes=favorable,
        neutral_outcomes=neutral,
        negative_outcomes=negative,
        average_confidence=round(avg_confidence, 3),
        evidence_strength=evidence_strength,
        conditions=[condition],
        cannabinoids=list(cannabinoids)
    )


@router.get("/aggregate", response_model=AggregatedEvidence)
async def get_aggregated_evidence(
    condition: str = Query(..., description="Condition to aggregate evidence for"),
    include_studies: bool = Query(True, description="Include individual study details"),
    limit: int = Query(50, le=200, description="Maximum studies to include"),
    db: Session = Depends(get_db)
):
    """Get aggregated evidence with confidence weighting for a condition.
    
    Implements evidence aggregation with:
    - Confidence-weighted efficacy scoring
    - Cannabinoid breakdown
    - Evidence grade calculation
    - Clinical recommendation
    """
    # Query studies
    studies = db.query(ClinicalStudy).filter(
        ClinicalStudy.condition.ilike(f"%{condition}%")
    ).order_by(ClinicalStudy.confidence_weight.desc()).limit(limit).all()
    
    if not studies:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"No studies found for condition '{condition}'"
        )
    
    # Process studies
    study_list = []
    cannabinoid_data = {}
    total_weighted_score = 0.0
    total_weight = 0.0
    has_rct = False
    
    for study in studies:
        confidence = study.confidence_weight or 0.5
        direction = parse_outcome_direction(study.outcome_summary)
        
        # Score: favorable=1.0, neutral=0.5, negative=0.0
        outcome_score = {"favorable": 1.0, "neutral": 0.5, "negative": 0.0}[direction]
        
        total_weighted_score += outcome_score * confidence
        total_weight += confidence
        
        # Check for RCT
        if study.study_design and "RCT" in study.study_design.upper():
            has_rct = True
        
        # Track cannabinoid data
        if study.cannabinoids_studied:
            for cb in study.cannabinoids_studied:
                if cb not in cannabinoid_data:
                    cannabinoid_data[cb] = {
                        "study_count": 0,
                        "favorable": 0,
                        "total_confidence": 0.0
                    }
                cannabinoid_data[cb]["study_count"] += 1
                cannabinoid_data[cb]["total_confidence"] += confidence
                if direction == "favorable":
                    cannabinoid_data[cb]["favorable"] += 1
        
        if include_studies:
            study_list.append(StudyEvidence(
                study_id=study.id,
                title=study.title,
                doi=study.doi,
                pubmed_id=study.pubmed_id,
                condition=study.condition,
                study_design=study.study_design or "Unknown",
                sample_size=study.sample_size,
                cannabinoids_studied=study.cannabinoids_studied or [],
                outcome_direction=direction,
                outcome_summary=study.outcome_summary,
                confidence_weight=confidence,
                year=study.year
            ))
    
    # Calculate weighted efficacy score
    weighted_efficacy = total_weighted_score / total_weight if total_weight > 0 else 0
    
    # Calculate confidence interval (simplified)
    # Using Wilson score interval approximation
    n = len(studies)
    p = weighted_efficacy
    z = 1.96  # 95% confidence
    
    denominator = 1 + z**2 / n
    center = (p + z**2 / (2 * n)) / denominator
    spread = z * ((p * (1 - p) / n + z**2 / (4 * n**2)) ** 0.5) / denominator
    
    confidence_interval = {
        "lower": max(0, round(center - spread, 3)),
        "upper": min(1, round(center + spread, 3))
    }
    
    # Evidence grade
    evidence_grade = calculate_evidence_grade(weighted_efficacy, n, has_rct)
    
    # Find top cannabinoid
    top_cannabinoid = None
    if cannabinoid_data:
        top_cannabinoid = max(
            cannabinoid_data.keys(),
            key=lambda x: cannabinoid_data[x]["favorable"] / cannabinoid_data[x]["study_count"]
            if cannabinoid_data[x]["study_count"] > 0 else 0
        )
    
    # Generate recommendation
    recommendation = generate_recommendation(
        evidence_grade,
        condition,
        top_cannabinoid,
        weighted_efficacy
    )
    
    # Calculate cannabinoid efficacy ratios
    for cb, data in cannabinoid_data.items():
        if data["study_count"] > 0:
            data["favorable_ratio"] = round(data["favorable"] / data["study_count"], 3)
            data["avg_confidence"] = round(data["total_confidence"] / data["study_count"], 3)
    
    return AggregatedEvidence(
        condition=condition,
        total_studies=len(studies),
        weighted_efficacy_score=round(weighted_efficacy, 3),
        confidence_interval=confidence_interval,
        evidence_grade=evidence_grade,
        studies=study_list if include_studies else [],
        cannabinoid_breakdown=cannabinoid_data,
        recommendation=recommendation
    )


@router.get("/search")
async def search_evidence(
    query: str = Query(..., min_length=2, description="Search query"),
    conditions: Optional[List[str]] = Query(None, description="Filter by conditions"),
    cannabinoids: Optional[List[str]] = Query(None, description="Filter by cannabinoids"),
    min_year: Optional[int] = Query(None, description="Minimum publication year"),
    min_confidence: float = Query(0.0, ge=0.0, le=1.0),
    limit: int = Query(50, le=200),
    db: Session = Depends(get_db)
):
    """Search clinical evidence across multiple fields.
    
    Searches in: title, condition, outcome summary, cannabinoids studied.
    """
    # Build query
    base_query = db.query(ClinicalStudy)
    
    # Text search
    search_filter = or_(
        ClinicalStudy.title.ilike(f"%{query}%"),
        ClinicalStudy.condition.ilike(f"%{query}%"),
        ClinicalStudy.outcome_summary.ilike(f"%{query}%")
    )
    base_query = base_query.filter(search_filter)
    
    # Apply filters
    if conditions:
        condition_filters = [ClinicalStudy.condition.ilike(f"%{c}%") for c in conditions]
        base_query = base_query.filter(or_(*condition_filters))
    
    if min_year:
        base_query = base_query.filter(ClinicalStudy.year >= min_year)
    
    if min_confidence > 0:
        base_query = base_query.filter(ClinicalStudy.confidence_weight >= min_confidence)
    
    # Execute
    studies = base_query.order_by(
        ClinicalStudy.confidence_weight.desc()
    ).limit(limit).all()
    
    # Format results
    results = []
    for study in studies:
        results.append({
            "id": study.id,
            "title": study.title,
            "condition": study.condition,
            "year": study.year,
            "cannabinoids": study.cannabinoids_studied or [],
            "outcome_direction": parse_outcome_direction(study.outcome_summary),
            "confidence": study.confidence_weight or 0.5,
            "doi": study.doi
        })
    
    return {
        "query": query,
        "total_results": len(results),
        "results": results
    }


@router.get("/stats/summary")
async def get_evidence_statistics(
    db: Session = Depends(get_db)
):
    """Get overall evidence database statistics.
    
    Returns counts by condition, study design, and evidence quality.
    """
    total_studies = db.query(ClinicalStudy).count()
    
    # Count by condition
    condition_counts = db.query(
        ClinicalStudy.condition,
        func.count(ClinicalStudy.id)
    ).group_by(ClinicalStudy.condition).all()
    
    # Count by study design
    design_counts = db.query(
        ClinicalStudy.study_design,
        func.count(ClinicalStudy.id)
    ).group_by(ClinicalStudy.study_design).all()
    
    # Average confidence
    avg_confidence = db.query(
        func.avg(ClinicalStudy.confidence_weight)
    ).scalar() or 0.5
    
    # Studies with high confidence
    high_confidence = db.query(ClinicalStudy).filter(
        ClinicalStudy.confidence_weight >= 0.7
    ).count()
    
    return {
        "total_studies": total_studies,
        "conditions": {c[0]: c[1] for c in condition_counts if c[0]},
        "study_designs": {d[0] or "Unknown": d[1] for d in design_counts},
        "average_confidence": round(avg_confidence, 3),
        "high_confidence_studies": high_confidence,
        "high_confidence_percentage": round(high_confidence / total_studies * 100, 1) if total_studies > 0 else 0
    }


@router.get("/citation/{study_id}")
async def get_study_citation(
    study_id: int,
    format: str = Query("apa", description="Citation format: apa, mla, bibtex"),
    db: Session = Depends(get_db)
):
    """Generate formatted citation for a study.
    
    Supports APA, MLA, and BibTeX formats.
    """
    study = db.query(ClinicalStudy).filter(ClinicalStudy.id == study_id).first()
    
    if not study:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"Study {study_id} not found"
        )
    
    authors = study.authors or "Unknown Authors"
    title = study.title
    year = study.year or "n.d."
    journal = study.journal or "Unknown Journal"
    doi = study.doi
    
    if format.lower() == "apa":
        citation = f"{authors} ({year}). {title}. {journal}."
        if doi:
            citation += f" https://doi.org/{doi}"
    elif format.lower() == "mla":
        citation = f'{authors}. "{title}." {journal}, {year}.'
        if doi:
            citation += f" doi:{doi}"
    elif format.lower() == "bibtex":
        key = f"study{study_id}"
        citation = f"""@article{{{key},
  author = {{{authors}}},
  title = {{{title}}},
  journal = {{{journal}}},
  year = {{{year}}},
  doi = {{{doi or ""}}}
}}"""
    else:
        citation = f"{authors} ({year}). {title}. {journal}."
    
    return {
        "study_id": study_id,
        "format": format,
        "citation": citation
    }
