from typing import Any, Dict, List, Optional


class _StubInflammatorySynergyEngine:
    """TS-PS-001 development stub with heuristic scoring and personalization."""

    MARKER_THRESHOLDS: Dict[str, Dict[str, float]] = {
        'tnf_alpha': {'moderate': 5.0, 'severe': 20.0},
        'il6': {'moderate': 2.0, 'severe': 10.0},
        'crp': {'moderate': 1.0, 'severe': 10.0},
        'il1b': {'moderate': 0.8, 'severe': 5.0},
    }
    MARKER_WEIGHTS: Dict[str, float] = {
        'tnf_alpha': 0.4,
        'il6': 0.25,
        'crp': 0.25,
        'il1b': 0.1,
    }
    EXPERIENCE_ADJUSTMENTS: Dict[str, float] = {
        'naive': -0.15,
        'beginner': -0.08,
        'intermediate': 0.0,
        'regular': 0.05,
        'experienced': 0.08,
    }
    KINGDOM_COMPOUNDS: Dict[str, List[str]] = {
        'cannabis': ['CBD', 'CBG', 'β-caryophyllene'],
        'fungal': ["Lion's Mane β-glucan", 'Reishi extract'],
        'marine': ['Fucoidan', 'Astaxanthin'],
        'plant': ['Curcumin', 'Quercetin'],
    }

    def predict_inflammatory_synergy(self, *args, **kwargs) -> Dict[str, Any]:
        if args:
            biomarkers = args[0] if len(args) > 0 else kwargs.get('biomarkers', {})
            condition_profile = args[1] if len(args) > 1 else kwargs.get('condition_profile', {})
            kingdoms = args[2] if len(args) > 2 else kwargs.get('available_kingdoms')
        else:
            biomarkers = kwargs.get('biomarkers', {})
            condition_profile = kwargs.get('condition_profile', {})
            kingdoms = kwargs.get('available_kingdoms')

        available_kingdoms = self._normalize_kingdoms(kingdoms)
        biomarker_scores = self._calculate_biomarker_scores(biomarkers)
        biomarker_score = self._weighted_average(biomarker_scores)
        condition_score = self._condition_score(condition_profile)
        goal_alignment = self._goal_alignment(condition_profile)
        experience_level = (condition_profile.get('experience_level') or 'beginner').lower()
        experience_adjustment = self.EXPERIENCE_ADJUSTMENTS.get(experience_level, -0.05)

        synergy_base = (0.55 * biomarker_score) + (0.3 * condition_score) + (0.15 * goal_alignment)
        synergy_score = self._clamp(synergy_base + experience_adjustment, 0.05, 0.99)

        kingdom_scores = self._kingdom_scores(biomarker_scores, condition_profile, available_kingdoms)
        primary_kingdom = self._select_primary_kingdom(kingdom_scores, available_kingdoms)
        secondary_kingdoms = self._select_secondary_kingdoms(kingdom_scores, primary_kingdom)
        recommended_compounds = self.KINGDOM_COMPOUNDS.get(primary_kingdom, ['Standard Botanical Extract'])
        dosing_guidance = self._dosing_guidance(experience_level, condition_profile, biomarkers)
        expected_reduction = self._expected_reduction(biomarkers, synergy_score)

        warning = None
        if not any((value or 0) for value in biomarkers.values()):
            warning = "Biomarkers appear empty — results are a heuristic fallback."
        elif biomarker_score < 0.15 and condition_score < 0.2:
            warning = "Limited inflammatory signals detected; leaning on patient history."

        confidence_level = self._clamp(0.55 + (synergy_score * 0.4), 0.55, 0.98)

        return {
            'primary_kingdom': primary_kingdom,
            'secondary_kingdoms': secondary_kingdoms,
            'synergy_score': round(synergy_score, 4),
            'confidence_level': round(confidence_level, 4),
            'recommended_compounds': recommended_compounds,
            'dosing_guidance': dosing_guidance,
            'expected_reduction': expected_reduction,
            'warning': warning,
        }

    # ------------------------------------------------------------------
    # Helper methods

    def _normalize_kingdoms(self, kingdoms: Optional[List[str]]) -> List[str]:
        if not kingdoms:
            return ['cannabis', 'fungal', 'marine', 'plant']
        seen: List[str] = []
        for k in kingdoms:
            if not isinstance(k, str):
                continue
            normalized = k.lower().strip()
            if normalized and normalized not in seen:
                seen.append(normalized)
        return seen or ['cannabis']

    def _calculate_biomarker_scores(self, biomarkers: Any) -> Dict[str, float]:
        scores: Dict[str, float] = {}
        data = biomarkers if isinstance(biomarkers, dict) else {}
        for marker, thresholds in self.MARKER_THRESHOLDS.items():
            value = data.get(marker) or 0.0
            moderate = thresholds['moderate']
            severe = thresholds['severe']
            if severe <= 0:
                scores[marker] = 0.0
                continue
            if value <= 0:
                norm = 0.0
            elif value <= moderate:
                norm = (value / max(moderate, 1e-6)) * 0.5
            elif value >= severe:
                norm = 1.0
            else:
                norm = 0.5 + ((value - moderate) / max(severe - moderate, 1e-6)) * 0.5
            scores[marker] = self._clamp(norm, 0.0, 1.0)
        return scores

    def _weighted_average(self, scores: Dict[str, float]) -> float:
        total_weight = sum(self.MARKER_WEIGHTS.values())
        if total_weight == 0:
            return 0.0
        aggregate = 0.0
        for marker, weight in self.MARKER_WEIGHTS.items():
            aggregate += weight * scores.get(marker, 0.0)
        return aggregate / total_weight

    def _condition_score(self, profile: Dict[str, Any]) -> float:
        conditions = profile.get('conditions') if isinstance(profile, dict) else None
        if not isinstance(conditions, list) or not conditions:
            return 0.0
        capped = []
        primary_bonus = 0.0
        for condition in conditions:
            severity = min(10, max(0, int(condition.get('severity', 0)))) if isinstance(condition, dict) else 0
            capped.append(severity)
            if condition.get('is_primary') and severity >= 7:
                primary_bonus = max(primary_bonus, 0.12)
        base_score = (sum(capped) / (10 * len(capped))) if capped else 0.0
        return self._clamp(base_score + primary_bonus, 0.0, 1.0)

    def _goal_alignment(self, profile: Dict[str, Any]) -> float:
        goal = (profile.get('primary_goal') or '').lower() if isinstance(profile, dict) else ''
        if not goal:
            return 0.35
        if 'inflamm' in goal:
            return 1.0
        if 'pain' in goal or 'autoimmune' in goal:
            return 0.75
        if 'sleep' in goal or 'anxiety' in goal:
            return 0.5
        return 0.4

    def _kingdom_scores(self, biomarker_scores: Dict[str, float], profile: Dict[str, Any], kingdoms: List[str]) -> Dict[str, float]:
        scores = {k: 0.0 for k in kingdoms}
        keywords = ' '.join(
            condition.get('name', '')
            for condition in (profile.get('conditions') or [])
            if isinstance(condition, dict)
        ).lower()

        for kingdom in scores:
            if kingdom == 'cannabis':
                scores[kingdom] += (biomarker_scores.get('tnf_alpha', 0.0) * 0.6)
                scores[kingdom] += (biomarker_scores.get('il6', 0.0) * 0.3)
                if 'pain' in keywords:
                    scores[kingdom] += 0.1
            elif kingdom == 'fungal':
                scores[kingdom] += (biomarker_scores.get('il1b', 0.0) * 0.5)
                scores[kingdom] += (self._condition_score(profile) * 0.3)
                if 'immune' in keywords or 'auto' in keywords:
                    scores[kingdom] += 0.15
            elif kingdom == 'marine':
                scores[kingdom] += (biomarker_scores.get('crp', 0.0) * 0.7)
                if 'cardio' in keywords or 'vascular' in keywords:
                    scores[kingdom] += 0.1
            elif kingdom == 'plant':
                scores[kingdom] += (biomarker_scores.get('il6', 0.0) * 0.4)
                if 'gut' in keywords or 'digest' in keywords:
                    scores[kingdom] += 0.1
        return scores

    def _select_primary_kingdom(self, scores: Dict[str, float], available: List[str]) -> str:
        if not scores:
            return available[0] if available else 'cannabis'
        ordered = sorted(scores.items(), key=lambda item: item[1], reverse=True)
        top_score = ordered[0][1]
        tied = [k for k, v in ordered if v == top_score]
        for candidate in ordered:
            if candidate[1] == top_score:
                # prefer cannabis when tied for regulatory familiarity
                if candidate[0] == 'cannabis' or len(tied) == 1:
                    return candidate[0]
        return ordered[0][0]

    def _select_secondary_kingdoms(self, scores: Dict[str, float], primary: str) -> List[str]:
        ordered = sorted(
            ((k, v) for k, v in scores.items() if k != primary),
            key=lambda item: item[1],
            reverse=True,
        )
        return [k for k, v in ordered if v > 0][:2]

    def _dosing_guidance(self, experience_level: str, profile: Dict[str, Any], biomarkers: Dict[str, Any]) -> Dict[str, Any]:
        preferences = profile.get('administration_preferences') if isinstance(profile, dict) else None
        preferred_route = (preferences[0] if preferences else 'oral_tincture') or 'oral_tincture'
        preferred_route = preferred_route.replace(' ', '_')

        low_experience = {'naive', 'beginner'}
        primary_plan = 'start_low_titrate' if experience_level in low_experience else 'standard_titrate'
        frequency = 'as_needed_microdose' if 'inhale' in preferred_route or 'vape' in preferred_route else 'twice_daily_with_food'

        notes: List[str] = []
        if experience_level in low_experience:
            notes.append('Limited cannabinoid experience noted; begin with micro-dosing titration.')
        if any((biomarkers.get(marker) or 0) > self.MARKER_THRESHOLDS[marker]['severe'] for marker in self.MARKER_THRESHOLDS):
            notes.append('Severe biomarker elevation detected; schedule follow-up labs in 14 days.')

        history_summary = self._history_summary(profile)
        guidance: Dict[str, Any] = {
            'primary': primary_plan,
            'frequency': frequency,
            'preferred_route': preferred_route,
            'notes': notes,
        }
        if history_summary:
            guidance['history_summary'] = history_summary
        return guidance

    def _history_summary(self, profile: Dict[str, Any]) -> Optional[str]:
        conditions = profile.get('conditions') if isinstance(profile, dict) else None
        if not isinstance(conditions, list) or not conditions:
            return None
        ordered = sorted(
            [c for c in conditions if isinstance(c, dict)],
            key=lambda c: c.get('severity', 0),
            reverse=True,
        )
        top = ordered[:2]
        if not top:
            return None
        parts = [f"{c.get('name', 'condition')} (sev {c.get('severity', 0)})" for c in top]
        return ', '.join(parts)

    def _expected_reduction(self, biomarkers: Dict[str, Any], synergy_score: float) -> Dict[str, float]:
        reduction: Dict[str, float] = {}
        multiplier = 0.25 + (synergy_score * 0.5)
        for marker, value in biomarkers.items():
            if not value:
                continue
            reduction[marker] = round(value * multiplier, 2)
        if not reduction:
            reduction['tnf_alpha'] = round(10 * synergy_score, 2)
            reduction['crp'] = round(5 * synergy_score, 2)
        return reduction

    @staticmethod
    def _clamp(value: float, lower: float, upper: float) -> float:
        return max(lower, min(upper, value))


def get_inflammatory_synergy_engine() -> _StubInflammatorySynergyEngine:
    return _StubInflammatorySynergyEngine()
