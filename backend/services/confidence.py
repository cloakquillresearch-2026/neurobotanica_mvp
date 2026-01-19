"""Confidence scoring utilities for clinical studies.

Provides:
- compute_confidence_for_study(study) -> float (0.0-1.0)
- compute_overall_confidence(condition: Optional[str]) -> float (0.0-1.0)

These are lightweight heuristics used by the recommendation engine and tests.
"""
from typing import Optional
from math import sqrt
from datetime import datetime

from backend.models.study import ClinicalStudy
from backend.models.database import SessionLocal


def compute_confidence_for_study(study: ClinicalStudy) -> float:
    """Compute a heuristic confidence (0.0-1.0) for a single ClinicalStudy.

    Heuristic basis:
    - Study design weighting (RCTs highest, preclinical lowest)
    - Sample size scaling (sqrt to reduce extreme influence)
    - Recency penalty (older studies slightly down-weighted)
    - Pivotal trials get a small boost
    """
    if study is None:
        return 0.0

    # Base score by study_type / evidence
    stype = (study.study_type or "").upper()
    base_map = {
        "RCT": 0.95,
        "SYSTEMATIC REVIEW": 0.9,
        "META-ANALYSIS": 0.9,
        "OBSERVATIONAL": 0.65,
        "CASE SERIES": 0.4,
        "PRECLINICAL": 0.2,
    }
    base = base_map.get(stype, 0.5)

    # Sample size factor (try to parse ints, tolerate strings)
    sample = 30
    try:
        if study.sample_size:
            # Some records store ranges or text like 'n=30' or '30 (mean)'
            s = str(study.sample_size)
            digits = "".join([c for c in s if c.isdigit()])
            if digits:
                sample = max(1, int(digits))
    except Exception:
        sample = 30

    size_factor = min(1.0, sqrt(sample) / sqrt(100))  # cap at 1.0

    # Recency factor (newer = slightly higher)
    year = study.year or 2000
    current = datetime.now().year
    age = max(0, current - int(year))
    recency = max(0.5, 1.0 - (age * 0.01))  # lose 1% per year, floor at 0.5

    # Pivotal trial boost
    pivotal = 1.05 if getattr(study, "is_pivotal_trial", False) else 1.0

    raw = base * size_factor * recency * pivotal
    # Bound and return
    return max(0.0, min(1.0, float(raw)))


def compute_overall_confidence(condition: Optional[str] = None) -> float:
    """Compute an overall confidence for a condition (average of studies).

    If `condition` is None, compute across all studies. Returns 0.0-1.0.
    """
    db = SessionLocal()
    try:
        q = db.query(ClinicalStudy)
        if condition:
            q = q.filter(ClinicalStudy.condition == condition.upper())
        studies = q.limit(500).all()
        if not studies:
            # default conservative baseline
            return 0.45

        weights = []
        for s in studies:
            try:
                w = compute_confidence_for_study(s)
                # allow per-record confidence_weight override
                if getattr(s, "confidence_weight", None) is not None:
                    w = (w + float(s.confidence_weight)) / 2.0
                weights.append(max(0.0, min(1.0, w)))
            except Exception:
                continue

        if not weights:
            return 0.45

        return sum(weights) / len(weights)
    finally:
        db.close()
