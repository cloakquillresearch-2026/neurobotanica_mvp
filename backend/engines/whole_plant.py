"""Whole-plant analysis engine for the NeuroBotanica platform."""
from __future__ import annotations

import json
import sqlite3
from typing import Any, Dict, List, Optional, Sequence, Tuple

from fastapi import APIRouter, Depends

from backend.models.compound_profile import CompoundProfile, WholePlantAnalysisResult
from backend.database import get_db_connection
from backend.services.terpene_synergy import TerpeneSynergyScorer


router = APIRouter(prefix="/api/v1/analysis", tags=["analysis"])


@router.post("/whole-plant", response_model=WholePlantAnalysisResult)
async def analyze_whole_plant(
    profile: CompoundProfile,
    db: sqlite3.Connection = Depends(get_db_connection),
) -> WholePlantAnalysisResult:
    """API endpoint wrapper for the whole-plant analyzer."""
    analyzer = WholePlantAnalyzer(db)
    return analyzer.analyze(profile)


class WholePlantAnalyzer:
    """Analyzes cannabinoid and terpene profiles to generate recommendations."""

    def __init__(self, db_connection: Optional[sqlite3.Connection]):
        self.db = db_connection
        self.terpene_scorer = TerpeneSynergyScorer(db_connection)
        self.compound_metadata: Dict[str, Dict[str, Any]] = {}
        self.condition_targets = self._build_condition_targets()
        self._load_compound_data()

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------
    def analyze(self, profile: CompoundProfile) -> WholePlantAnalysisResult:
        cannabinoid_analysis = self._analyze_cannabinoids(
            profile.cannabinoid_profile, profile.target_condition
        )
        terpene_synergy, terpene_confidence = self.terpene_scorer.calculate_synergy_score(
            profile.terpene_profile, profile.target_condition
        )
        entourage_compounds = self._get_entourage_effects(
            profile.cannabinoid_profile, profile.terpene_profile, profile.target_condition
        )
        clinical_evidence = self._get_clinical_evidence(
            profile.target_condition, profile.cannabinoid_profile
        )
        overall_synergy = self._calculate_overall_synergy(
            cannabinoid_analysis["score"], terpene_synergy, entourage_compounds
        )
        confidence = self._calculate_confidence(
            terpene_confidence,
            len(clinical_evidence),
            entourage_compounds,
            cannabinoid_analysis.get("profile_completeness", 0.5),
        )
        recommended_conditions = self._generate_recommendations(
            profile, cannabinoid_analysis["normalized_profile"], terpene_synergy
        )

        return WholePlantAnalysisResult(
            synergy_score=overall_synergy,
            confidence=confidence,
            recommended_conditions=recommended_conditions,
            entourage_compounds=entourage_compounds,
            clinical_evidence=clinical_evidence,
        )

    # ------------------------------------------------------------------
    # Data loading
    # ------------------------------------------------------------------
    def _load_compound_data(self) -> None:
        if not self.db:
            return
        try:
            cursor = self.db.cursor()
            cursor.execute(
                "SELECT name FROM sqlite_master WHERE type='table' AND name='neurobotanica_compounds'"
            )
            if not cursor.fetchone():
                raise RuntimeError("Missing required table: neurobotanica_compounds")
            cursor.execute(
                "SELECT compound_id, compound_name, therapeutic_targets, primary_mechanisms FROM neurobotanica_compounds"
            )
            for compound_id, name, targets, mechanisms in cursor.fetchall():
                conditions = self._parse_condition_list(targets)
                self.compound_metadata[(compound_id or "").lower()] = {
                    "name": name or compound_id,
                    "conditions": conditions,
                    "mechanisms": mechanisms,
                }
        except sqlite3.Error:
            self.compound_metadata = {}

    # ------------------------------------------------------------------
    # Cannabinoid analysis
    # ------------------------------------------------------------------
    def _analyze_cannabinoids(
        self, cannabinoid_profile: Dict[str, float], target_condition: str
    ) -> Dict[str, Any]:
        normalized = self._normalize(cannabinoid_profile)
        normalized_lower = {k.lower(): v for k, v in normalized.items()}
        score = self._score_condition(normalized_lower, target_condition.lower())

        supporting = []
        for compound, proportion in sorted(
            normalized_lower.items(), key=lambda item: item[1], reverse=True
        )[:5]:
            supporting.append(
                {
                    "compound": compound,
                    "proportion": round(proportion, 3),
                    "conditions": list(self.compound_metadata.get(compound, {}).get("conditions", [])),
                }
            )

        completeness = min(len(cannabinoid_profile) / 6.0, 1.0)
        return {
            "score": score,
            "normalized_profile": normalized_lower,
            "dominant_compounds": supporting,
            "profile_completeness": completeness,
        }

    def _score_condition(
        self, normalized_profile: Dict[str, float], condition: str
    ) -> float:
        rule = self.condition_targets.get(condition, {})
        if not rule:
            return 0.5

        ratio_score = self._ratio_score(
            normalized_profile, rule.get("preferred_ratio"), rule.get("max_thc")
        )
        ratio_floor = rule.get("ratio_floor")
        if ratio_floor is not None:
            ratio_score = max(ratio_score, ratio_floor)
        boost_score = 0.0
        for compound, weight in rule.get("boost", {}).items():
            boost_score += weight * normalized_profile.get(compound, 0.0)

        metadata_boost = 0.0
        for compound, value in normalized_profile.items():
            metadata = self.compound_metadata.get(compound)
            if metadata and condition in metadata.get("conditions", set()):
                metadata_boost += 0.1 * value

        adaptive_bonus = 0.0
        thc_bonus = rule.get("thc_dominant_bonus")
        if thc_bonus:
            thc_value = normalized_profile.get("thc", 0.0)
            threshold = thc_bonus.get("threshold", 0.6)
            weight = thc_bonus.get("weight", 0.2)
            if thc_value >= threshold:
                adaptive_bonus += weight * thc_value

        score = (0.45 * ratio_score) + (0.4 * boost_score) + metadata_boost + adaptive_bonus
        return max(0.0, min(score, 1.0))

    def _ratio_score(
        self,
        normalized_profile: Dict[str, float],
        ratio_rule: Optional[Tuple[str, str, float]],
        max_thc: Optional[float],
    ) -> float:
        if not ratio_rule:
            return 0.5
        a_key, b_key, target = ratio_rule
        a_val = normalized_profile.get(a_key, 0.0)
        b_val = normalized_profile.get(b_key, 0.0)
        if b_val == 0 and a_val == 0:
            return 0.0
        if b_val == 0:
            ratio = target
        else:
            ratio = a_val / max(b_val, 1e-4)
        deviation = abs(ratio - target) / (target + 1e-3)
        ratio_score = max(0.0, 1.0 - deviation)
        if max_thc is not None:
            thc_value = normalized_profile.get("thc", 0.0)
            if thc_value > max_thc:
                ratio_score *= max(0.2, 1 - (thc_value - max_thc))
        return ratio_score

    # ------------------------------------------------------------------
    # Database helpers
    # ------------------------------------------------------------------
    def _get_entourage_effects(
        self,
        cannabinoid_profile: Dict[str, float],
        terpene_profile: Dict[str, float],
        condition: str,
    ) -> List[Dict[str, Any]]:
        if not self.db:
            return []

        compounds = {
            name.lower()
            for name, value in cannabinoid_profile.items()
            if value is not None and value > 0.0
        }
        compounds.update(
            name.lower()
            for name, value in terpene_profile.items()
            if value is not None and value > 0.0
        )
        if not compounds:
            return []

        placeholders = ",".join("?" for _ in compounds)
        try:
            cursor = self.db.cursor()
            cursor.execute(
                f"""
                SELECT compound_a_id, compound_b_id, synergy_score, mechanism,
                       therapeutic_context, confidence_score
                FROM neurobotanica_synergy_predictions
                WHERE compound_a_id IN ({placeholders})
                   OR compound_b_id IN ({placeholders})
                """,
                tuple(compounds) + tuple(compounds),
            )
        except sqlite3.Error:
            return []

        matches = []
        seen = set()
        for row in cursor.fetchall():
            a_id, b_id, score, mechanism, context, confidence = row
            a_name = (a_id or "").lower()
            b_name = (b_id or "").lower()
            if not a_name or not b_name or a_name not in compounds or b_name not in compounds:
                continue
            pair_key = tuple(sorted([a_name, b_name]))
            if pair_key in seen:
                continue
            seen.add(pair_key)
            supports_condition = condition.lower() in self._parse_condition_list(context)
            matches.append(
                {
                    "pair": f"{a_id}+{b_id}",
                    "compound_a": a_id,
                    "compound_b": b_id,
                    "synergy_score": float(score or 0.0),
                    "mechanism": mechanism or "entourage",
                    "supports_condition": supports_condition,
                    "confidence": float(confidence or 0.6),
                }
            )

        matches.sort(key=lambda item: item["synergy_score"], reverse=True)
        return matches[:6]

    def _get_clinical_evidence(
        self, condition: str, cannabinoid_profile: Dict[str, float]
    ) -> List[Dict[str, Any]]:
        if not self.db:
            return []
        try:
            cursor = self.db.cursor()
            cursor.execute(
                "SELECT name FROM sqlite_master WHERE type='table' AND name='neurobotanica_clinical_studies'"
            )
            if not cursor.fetchone():
                return []
            cursor.execute(
                """
                SELECT study_id, study_type, intervention, key_findings, citation,
                       confidence_score, sample_size, publication_year
                FROM neurobotanica_clinical_studies
                WHERE LOWER(condition) = ?
                ORDER BY confidence_score DESC, sample_size DESC
                LIMIT 6
                """,
                (condition.lower(),),
            )
        except sqlite3.Error:
            return []

        evidence = []
        profile_names = {name.lower() for name in cannabinoid_profile.keys()}
        for row in cursor.fetchall():
            study_id, study_type, intervention, findings, citation, confidence, sample_size, year = row
            intervention_text = (intervention or "").lower()
            relevant = any(name in intervention_text for name in profile_names)
            if not relevant and profile_names:
                continue
            evidence.append(
                {
                    "study": study_id,
                    "study_id": study_id,
                    "type": study_type,
                    "citation": citation,
                    "summary": findings,
                    "finding": findings,
                    "confidence_score": float(confidence or 0.6),
                    "sample_size": sample_size,
                    "year": year,
                }
            )
            # Fallback: if no filtered studies, return top N by confidence_score
            if not evidence:
                cursor.execute(
                    """
                    SELECT study_id, study_type, intervention, key_findings, citation,
                           confidence_score, sample_size, publication_year
                    FROM neurobotanica_clinical_studies
                    WHERE LOWER(condition) = ?
                    ORDER BY confidence_score DESC
                    LIMIT ?
                    """,
                    (condition.lower(), top_n),
                )
                evidence = cursor.fetchall()
            return evidence
        return evidence

    # ------------------------------------------------------------------
    # Score calculations
    # ------------------------------------------------------------------
    @staticmethod
    def _calculate_overall_synergy(
        cannabinoid_score: float,
        terpene_score: float,
        entourage_data: Sequence[Dict[str, Any]],
    ) -> float:
        if entourage_data:
            entourage_score = sum(item["synergy_score"] for item in entourage_data) / len(entourage_data)
        else:
            entourage_score = 0.0
        combined = (0.4 * cannabinoid_score) + (0.3 * terpene_score) + (0.3 * entourage_score)
        return max(0.0, min(round(combined, 3), 1.0))

    @staticmethod
    def _calculate_confidence(
        terpene_confidence: float,
        evidence_count: int,
        entourage_data: Sequence[Dict[str, Any]],
        profile_completeness: float,
    ) -> float:
        evidence_strength = min(evidence_count / 12.0, 1.0)
        synergy_quality = (
            sum(item.get("confidence", 0.6) for item in entourage_data) / len(entourage_data)
            if entourage_data
            else 0.6
        )
        composite = (0.4 * evidence_strength) + (0.3 * profile_completeness) + (0.3 * synergy_quality)
        confidence = 0.5 + 0.4 * ((0.5 * terpene_confidence) + (0.5 * composite))
        return max(0.5, min(round(confidence, 3), 0.9))

    def _generate_recommendations(
        self,
        profile: CompoundProfile,
        normalized_cannabinoids: Dict[str, float],
        target_terpene_score: float,
    ) -> List[str]:
        scores: List[Tuple[str, float]] = []
        for condition in self.condition_targets.keys():
            cannabinoid_score = self._score_condition(normalized_cannabinoids, condition)
            terpene_score, _ = self.terpene_scorer.calculate_synergy_score(
                profile.terpene_profile, condition
            )
            combined = (0.4 * cannabinoid_score) + (0.3 * terpene_score) + (0.3 * target_terpene_score)
            scores.append((condition, combined))

        scores.sort(key=lambda item: item[1], reverse=True)
        recommendations = [condition for condition, score in scores if score >= 0.5][:6]
        if profile.target_condition not in recommendations:
            recommendations = [profile.target_condition] + [cond for cond in recommendations if cond != profile.target_condition]
        if len(recommendations) < 2:
            for condition, _ in scores:
                if condition not in recommendations:
                    recommendations.append(condition)
                if len(recommendations) >= 2:
                    break
        seen = set()
        ordered: List[str] = []
        for item in recommendations:
            if item not in seen:
                ordered.append(item)
                seen.add(item)
        if len(ordered) < 2:
            for condition, _ in scores:
                if condition not in ordered:
                    ordered.append(condition)
                if len(ordered) >= 2:
                    break
        return ordered

    # ------------------------------------------------------------------
    # Utility helpers
    # ------------------------------------------------------------------
    @staticmethod
    def _normalize(profile: Dict[str, float]) -> Dict[str, float]:
        total = sum(max(value, 0.0) for value in profile.values()) or 1.0
        return {name: max(value, 0.0) / total for name, value in profile.items()}

    @staticmethod
    def _parse_condition_list(raw_value: Optional[str]) -> List[str]:
        if not raw_value:
            return []
        try:
            parsed = json.loads(raw_value)
            if isinstance(parsed, list):
                return [str(item).lower() for item in parsed]
        except (json.JSONDecodeError, TypeError):
            pass
        return [part.strip().lower() for part in str(raw_value).split(",") if part.strip()]

    @staticmethod
    def _build_condition_targets() -> Dict[str, Dict[str, Any]]:
        return {
            "anxiety": {
                "preferred_ratio": ("cbd", "thc", 2.0),
                "boost": {"cbd": 0.6, "cbg": 0.2, "cbn": 0.1},
                "max_thc": 0.35,
            },
            "chronic_pain": {
                "preferred_ratio": ("thc", "cbd", 1.0),
                "boost": {"cbg": 0.2, "cbc": 0.15, "thcv": 0.15},
                "max_thc": 0.7,
                "ratio_floor": 0.2,
                "thc_dominant_bonus": {"threshold": 0.55, "weight": 0.35},
            },
            "sleep": {
                "preferred_ratio": ("cbn", "thc", 0.5),
                "boost": {"cbn": 0.4, "myrcene": 0.1, "linalool": 0.1},
                "max_thc": 0.6,
            },
            "ptsd": {
                "preferred_ratio": ("cbd", "thc", 1.5),
                "boost": {"cbd": 0.4, "thcv": 0.2, "cbg": 0.1},
                "max_thc": 0.5,
            },
            "inflammation": {
                "preferred_ratio": ("cbd", "thc", 1.2),
                "boost": {"cbc": 0.25, "cbg": 0.25, "beta_caryophyllene": 0.2},
                "max_thc": 0.6,
                "ratio_floor": 0.15,
            },
            "fibromyalgia": {
                "preferred_ratio": ("cbd", "thc", 1.8),
                "boost": {"cbd": 0.4, "cbg": 0.15, "cbn": 0.15},
                "max_thc": 0.5,
                "ratio_floor": 0.15,
            },
            "depression": {
                "preferred_ratio": ("thc", "cbd", 0.8),
                "boost": {"thcv": 0.2, "cbd": 0.2, "limonene": 0.1},
                "max_thc": 0.65,
                "ratio_floor": 0.10,
            },
            "migraine": {
                "preferred_ratio": ("thc", "cbd", 1.2),
                "boost": {"thc": 0.3, "cbg": 0.2, "beta_caryophyllene": 0.2},
                "max_thc": 0.7,
                "ratio_floor": 0.20,
                "thc_dominant_bonus": {"threshold": 0.50, "weight": 0.25},
            },
            "arthritis": {
                "preferred_ratio": ("cbd", "thc", 1.5),
                "boost": {"cbc": 0.2, "cbg": 0.2},
                "max_thc": 0.6,
                "ratio_floor": 0.15,
            },
            "neuropathy": {
                "preferred_ratio": ("thc", "cbd", 1.0),
                "boost": {"cbn": 0.2, "cbg": 0.2},
                "max_thc": 0.6,
                "ratio_floor": 0.15,
                "thc_dominant_bonus": {"threshold": 0.55, "weight": 0.25},
            },
            "epilepsy": {
                "preferred_ratio": ("cbd", "thc", 4.0),
                "boost": {"cbd": 0.7, "cbg": 0.1},
                "max_thc": 0.3,
            },
            "ibs": {
                "preferred_ratio": ("cbd", "thc", 1.5),
                "boost": {"cbd": 0.3, "cbg": 0.2},
                "max_thc": 0.5,
            },
            "crohns": {
                "preferred_ratio": ("thc", "cbd", 1.1),
                "boost": {"thc": 0.25, "cbd": 0.25},
                "max_thc": 0.65,
                "ratio_floor": 0.15,
            },
            "autism": {
                "preferred_ratio": ("cbd", "thc", 3.0),
                "boost": {"cbd": 0.5, "cbg": 0.15},
                "max_thc": 0.35,
            },
            "adhd": {
                "preferred_ratio": ("thc", "cbd", 0.9),
                "boost": {"thcv": 0.3, "cbd": 0.2},
                "max_thc": 0.6,
                "ratio_floor": 0.10,
            },
            "menstrual_pain": {
                "preferred_ratio": ("thc", "cbd", 1.2),
                "boost": {"thc": 0.3, "cbd": 0.2},
                "max_thc": 0.7,
                "ratio_floor": 0.20,
                "thc_dominant_bonus": {"threshold": 0.60, "weight": 0.30},
            },
            "cancer_pain": {
                "preferred_ratio": ("thc", "cbd", 1.5),
                "boost": {"thc": 0.35, "cbg": 0.2},
                "max_thc": 0.75,
                "ratio_floor": 0.25,
                "thc_dominant_bonus": {"threshold": 0.60, "weight": 0.30},
            },
            "parkinsons": {
                "preferred_ratio": ("cbd", "thc", 1.3),
                "boost": {"cbd": 0.35, "thcv": 0.2},
                "max_thc": 0.5,
            },
            "multiple_sclerosis": {
                "preferred_ratio": ("thc", "cbd", 1.0),
                "boost": {"thc": 0.3, "cbd": 0.3},
                "max_thc": 0.7,
                "ratio_floor": 0.15,
                "thc_dominant_bonus": {"threshold": 0.55, "weight": 0.25},
            },
            "alzheimer": {
                "preferred_ratio": ("cbd", "thc", 1.7),
                "boost": {"cbd": 0.4, "cbg": 0.2},
                "max_thc": 0.5,
            },
            "glaucoma": {
                "preferred_ratio": ("thc", "cbd", 1.4),
                "boost": {"thc": 0.35, "cbg": 0.15},
                "max_thc": 0.75,
                "ratio_floor": 0.20,
                "thc_dominant_bonus": {"threshold": 0.60, "weight": 0.30},
            },
            "appetite_loss": {
                "preferred_ratio": ("thc", "cbd", 1.8),
                "boost": {"thc": 0.4, "cbn": 0.2},
                "max_thc": 0.8,
                "ratio_floor": 0.25,
                "thc_dominant_bonus": {"threshold": 0.60, "weight": 0.30},
            },
        }
