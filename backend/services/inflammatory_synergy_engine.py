from __future__ import annotations

from typing import Any, Dict, List, Optional
import logging

from backend.routers import persistence as d1_persistence

logger = logging.getLogger(__name__)


class InflammatorySynergyEngine:
    """TS-PS-001 engine backed by D1 with heuristic fallback."""

    DEFAULT_KINGDOMS = ["cannabis", "fungal", "marine", "plant"]

    async def predict_inflammatory_synergy(self, *args, **kwargs) -> Dict[str, Any]:
        if args:
            biomarkers = args[0] if len(args) > 0 else kwargs.get("biomarkers", {})
            condition_profile = args[1] if len(args) > 1 else kwargs.get("condition_profile", {})
            kingdoms = args[2] if len(args) > 2 else kwargs.get("available_kingdoms")
        else:
            biomarkers = kwargs.get("biomarkers", {})
            condition_profile = kwargs.get("condition_profile", {})
            kingdoms = kwargs.get("available_kingdoms")

        available_kingdoms = self._normalize_kingdoms(kingdoms)
        primary_condition = self._primary_condition(condition_profile)

        try:
            condition_evidence = await self._fetch_condition_evidence(primary_condition)
            recommended_tokens = self._tokenize_recommendations(
                condition_evidence.get("recommended_cannabinoids") if condition_evidence else None
            )
            compound_rows = await self._resolve_compounds(recommended_tokens or ["CBD"])
            primary_kingdom = self._select_primary_kingdom(compound_rows, available_kingdoms)
            recommended_compounds = self._pick_recommended_compounds(compound_rows, primary_kingdom)

            synergy = await self._predict_synergy(compound_rows)
            confidence_level = self._combine_confidence(
                synergy.get("synergy_score", 0.5),
                condition_evidence.get("avg_confidence") if condition_evidence else None,
            )
            expected_reduction = self._expected_reduction(biomarkers, synergy.get("synergy_score", 0.5))

            return {
                "primary_kingdom": primary_kingdom,
                "secondary_kingdoms": [k for k in available_kingdoms if k != primary_kingdom][:2],
                "synergy_score": round(synergy.get("synergy_score", 0.5), 4),
                "confidence_level": round(confidence_level, 4),
                "recommended_compounds": recommended_compounds or ["CBD"],
                "dosing_guidance": self._dosing_guidance(primary_kingdom),
                "expected_reduction": expected_reduction,
                "warning": None,
            }
        except Exception as exc:
            logger.warning("D1-backed TS-PS-001 failed, using heuristic fallback: %s", exc)
            return self._fallback_prediction(biomarkers, condition_profile, available_kingdoms)

    async def _fetch_condition_evidence(self, condition: str) -> Dict[str, Any]:
        if not condition:
            return {}

        normalized = condition.upper()
        condition_row = await self._execute(
            "SELECT condition_name, category, recommended_cannabinoids, evidence_count FROM conditions WHERE UPPER(condition_name) = ? OR condition_name LIKE ?",
            [normalized, f"%{condition}%"],
        )
        condition_data = (condition_row.get("results") or [None])[0]

        studies = await self._execute(
            "SELECT study_id, study_type, citation, confidence_score FROM clinical_studies WHERE UPPER(condition) = ? OR condition LIKE ? ORDER BY confidence_score DESC LIMIT 10",
            [normalized, f"%{condition}%"],
        )
        study_rows = studies.get("results") or []
        avg_confidence = sum(row.get("confidence_score", 0) for row in study_rows) / max(1, len(study_rows))

        return {
            "condition": (condition_data or {}).get("condition_name", condition),
            "category": (condition_data or {}).get("category"),
            "recommended_cannabinoids": (condition_data or {}).get("recommended_cannabinoids"),
            "evidence_count": (condition_data or {}).get("evidence_count", len(study_rows)),
            "avg_confidence": round(avg_confidence, 2),
        }

    async def _resolve_compounds(self, tokens: List[str]) -> List[Dict[str, Any]]:
        if not tokens:
            return []
        placeholders = ",".join(["?"] * len(tokens))
        result = await self._execute(
            f"SELECT compound_id, compound_name, kingdom FROM neurobotanica_compounds WHERE compound_id IN ({placeholders}) OR compound_name IN ({placeholders}) LIMIT 20",
            tokens + tokens,
        )
        return result.get("results") or []

    async def _predict_synergy(self, compounds: List[Dict[str, Any]]) -> Dict[str, Any]:
        compound_ids = [c.get("compound_id") for c in compounds if c.get("compound_id")]
        if not compound_ids:
            return {"synergy_score": 0.5, "evidence": "Default prediction"}

        a = compound_ids[0]
        b = compound_ids[1] if len(compound_ids) > 1 else compound_ids[0]
        result = await self._execute(
            "SELECT synergy_score, clinical_evidence FROM neurobotanica_synergy_predictions WHERE (compound_a_id = ? AND compound_b_id = ?) OR (compound_a_id = ? AND compound_b_id = ?) ORDER BY confidence_score DESC LIMIT 1",
            [a, b, b, a],
        )
        row = (result.get("results") or [None])[0]
        if row:
            return {
                "synergy_score": row.get("synergy_score", 0.5),
                "evidence": row.get("clinical_evidence") or "Database prediction",
            }

        return {"synergy_score": 0.5, "evidence": "No database match"}

    def _combine_confidence(self, synergy_score: float, evidence_confidence: Optional[float]) -> float:
        base = evidence_confidence or 0.5
        return max(0.45, min(0.95, base + (synergy_score * 0.2)))

    def _expected_reduction(self, biomarkers: Dict[str, Any], synergy_score: float) -> Dict[str, float]:
        multiplier = 0.25 + (synergy_score * 0.5)
        reduction: Dict[str, float] = {}
        for marker, value in biomarkers.items():
            if not value:
                continue
            reduction[marker] = round(float(value) * multiplier, 2)
        if not reduction:
            reduction = {
                "tnf_alpha": round(10 * synergy_score, 2),
                "crp": round(5 * synergy_score, 2),
            }
        return reduction

    def _dosing_guidance(self, kingdom: str) -> Dict[str, str]:
        return {
            "cannabis": "25-50mg CBD daily",
            "fungal": "500mg extract daily",
            "marine": "200-400mg daily",
            "plant": "500-1000mg Curcumin daily",
        }.get(kingdom, "25-50mg CBD daily")

    def _tokenize_recommendations(self, raw: Optional[str]) -> List[str]:
        if not raw:
            return []
        return [token.strip() for token in raw.split(",") if token.strip()]

    def _select_primary_kingdom(self, compounds: List[Dict[str, Any]], available: List[str]) -> str:
        counts: Dict[str, int] = {}
        for compound in compounds:
            kingdom = (compound.get("kingdom") or "").lower()
            if kingdom:
                counts[kingdom] = counts.get(kingdom, 0) + 1
        if not counts:
            return available[0] if available else "cannabis"
        ordered = sorted(counts.items(), key=lambda item: item[1], reverse=True)
        for kingdom, _count in ordered:
            if kingdom in available:
                return kingdom
        return ordered[0][0]

    def _pick_recommended_compounds(self, compounds: List[Dict[str, Any]], kingdom: str) -> List[str]:
        return [
            c.get("compound_name") or c.get("compound_id")
            for c in compounds
            if (c.get("kingdom") or "").lower() == kingdom
        ][:3]

    def _normalize_kingdoms(self, kingdoms: Optional[List[str]]) -> List[str]:
        if not kingdoms:
            return self.DEFAULT_KINGDOMS
        seen: List[str] = []
        for k in kingdoms:
            if not isinstance(k, str):
                continue
            normalized = k.lower().strip()
            if normalized and normalized not in seen:
                seen.append(normalized)
        return seen or self.DEFAULT_KINGDOMS

    def _primary_condition(self, profile: Dict[str, Any]) -> str:
        conditions = profile.get("conditions") if isinstance(profile, dict) else None
        if not isinstance(conditions, list) or not conditions:
            return ""
        first = conditions[0]
        if isinstance(first, dict):
            return str(first.get("name") or "")
        if isinstance(first, str):
            return first
        return ""

    async def _execute(self, sql: str, params: List[Any]) -> Dict[str, Any]:
        try:
            return await d1_persistence.db.execute(sql, params)
        except Exception as exc:
            logger.warning("D1 execute failed: %s", exc)
            return {"results": []}

    def _fallback_prediction(
        self,
        biomarkers: Dict[str, Any],
        condition_profile: Dict[str, Any],
        available_kingdoms: List[str],
    ) -> Dict[str, Any]:
        primary_kingdom = available_kingdoms[0] if available_kingdoms else "cannabis"
        synergy_score = 0.55
        confidence_level = 0.6
        expected_reduction = self._expected_reduction(biomarkers, synergy_score)
        return {
            "primary_kingdom": primary_kingdom,
            "secondary_kingdoms": [k for k in available_kingdoms if k != primary_kingdom][:2],
            "synergy_score": round(synergy_score, 4),
            "confidence_level": round(confidence_level, 4),
            "recommended_compounds": ["CBD"],
            "dosing_guidance": self._dosing_guidance(primary_kingdom),
            "expected_reduction": expected_reduction,
            "warning": "D1 unavailable; using heuristic fallback.",
        }


def get_inflammatory_synergy_engine() -> InflammatorySynergyEngine:
    return InflammatorySynergyEngine()
