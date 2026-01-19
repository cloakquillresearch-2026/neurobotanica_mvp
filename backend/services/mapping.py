import json
import os
from typing import Dict

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
MAP_PATH = os.path.join(ROOT, 'backend', 'data', 'condition_cannabinoid_map.json')

_MAP: Dict[str, Dict[str, float]] = {}


def load_mapping() -> Dict[str, Dict[str, float]]:
    global _MAP
    if _MAP:
        return _MAP
    try:
        with open(MAP_PATH, 'r', encoding='utf-8') as f:
            _MAP = json.load(f)
    except Exception:
        _MAP = {}
    return _MAP


def get_cannabinoid_scores_for_condition(condition: str) -> Dict[str, float]:
    mapping = load_mapping()
    return mapping.get(condition, {})
