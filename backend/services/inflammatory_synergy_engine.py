from typing import Any, Dict, List


class _StubInflammatorySynergyEngine:
    """Minimal stub of TS-PS-001 engine for tests.

    This provides a deterministic, low-dependency implementation used by
    tests and local development. The real engine is proprietary and should
    be injected in production deployments.
    """

    def predict_inflammatory_synergy(self, *, biomarkers: Dict[str, float], condition_profile: Dict[str, Any], available_kingdoms: List[str]) -> Dict[str, Any]:
        # Simple heuristic: pick 'cannabis' as primary if TNF-α or CRP elevated,
        # otherwise choose the first available kingdom. Compute a mock score.
        tnf = biomarkers.get('tnf_alpha', 0) or 0
        crp = biomarkers.get('crp', 0) or 0
        il6 = biomarkers.get('il6', 0) or 0

        primary = 'cannabis' if (tnf > 5 or crp > 3 or il6 > 2) and 'cannabis' in available_kingdoms else (available_kingdoms[0] if available_kingdoms else 'cannabis')
        secondary = [k for k in available_kingdoms if k != primary][:2]

        # Mock scoring
        raw_score = min(1.0, (tnf * 0.02) + (crp * 0.03) + (il6 * 0.01))
        confidence = 0.6 + (min(0.4, raw_score * 0.4))

        recommended = []
        if primary == 'cannabis':
            recommended = ['CBD', 'CBG']
        elif primary == 'fungal':
            recommended = ['Beta-Glucan Extract']
        else:
            recommended = ['Standard Botanical Extract']

        dosing = {
            'primary': 'start_low_titrate',
            'frequency': 'twice_daily'
        }

        expected_reduction = {'tnf_alpha': raw_score * 10, 'crp': raw_score * 5}

        return {
            'primary_kingdom': primary,
            'secondary_kingdoms': secondary,
            'synergy_score': float(raw_score),
            'confidence_level': float(confidence),
            'recommended_compounds': recommended,
            'dosing_guidance': dosing,
            'expected_reduction': expected_reduction,
            'warning': None,
        }


def get_inflammatory_synergy_engine() -> _StubInflammatorySynergyEngine:
    return _StubInflammatorySynergyEngine()
