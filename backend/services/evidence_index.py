"""Evidence index utilities.

Builds condition -> evidence lookups and provides summary APIs used
by the recommendation engine and admin tooling.
"""
from typing import Dict, Any, List, Optional
from collections import defaultdict

from backend.models.study import ClinicalStudy
from backend.models.database import SessionLocal
from backend.services.confidence import compute_confidence_for_study


def build_condition_index(limit_per_condition: int = 500) -> Dict[str, Dict[str, Any]]:
    """Build an in-memory index mapping CONDITION -> summary and top studies.

    Returns a dict keyed by uppercase condition name with values:
      - total_studies
      - by_type: {study_type: count}
      - average_confidence
      - top_studies: list of study dicts (limited)
    """
    db = SessionLocal()
    try:
        conditions = defaultdict(list)
        for s in db.query(ClinicalStudy).all():
            cond = (s.condition or "").upper()
            conditions[cond].append(s)

        result = {}
        for cond, studies in conditions.items():
            by_type = defaultdict(int)
            confidences = []
            for st in studies:
                stype = (st.study_type or "UNKNOWN").upper()
                by_type[stype] += 1
                confidences.append(compute_confidence_for_study(st))

            avg_conf = sum(confidences) / len(confidences) if confidences else 0.0

            # pick top studies by confidence
            top = sorted(studies, key=lambda x: compute_confidence_for_study(x), reverse=True)[:limit_per_condition]
            top_serialized = [
                {
                    "study_id": t.study_id,
                    "title": t.study_title,
                    "year": t.year,
                    "study_type": t.study_type,
                    "confidence": compute_confidence_for_study(t),
                }
                for t in top
            ]

            result[cond] = {
                "total_studies": len(studies),
                "by_type": dict(by_type),
                "average_confidence": round(avg_conf, 4),
                "top_studies": top_serialized,
            }

        return result
    finally:
        db.close()


def get_condition_summary(condition: str, limit: int = 20) -> Optional[Dict[str, Any]]:
    """Return summary for a single condition (case-insensitive).

    If no studies found, returns None.
    """
    idx = build_condition_index(limit_per_condition=limit)
    return idx.get(condition.upper())
