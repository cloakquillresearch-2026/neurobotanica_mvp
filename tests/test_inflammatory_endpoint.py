import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from fastapi.testclient import TestClient
from backend.main import app

client = TestClient(app)

def test_inflammatory_endpoint_returns_fallback():
    payload = {
        "biomarkers": {"tnf_alpha": 20, "il6": 10, "crp": 5, "il1b": 2},
        "condition_profile": {"conditions": [{"name": "inflammation", "severity": 7, "is_primary": True}]},
        "available_kingdoms": ["cannabis", "fungal", "plant"]
    }
    resp = client.post('/api/dispensary/inflammatory-synergy', json=payload)
    assert resp.status_code == 200
    json_resp = resp.json()
    assert 'primary_kingdom' in json_resp
    assert 'confidence_level' in json_resp


def test_inflammatory_endpoint_empty_biomarkers_returns_fallback():
    payload = {
        "biomarkers": {},
        "condition_profile": {"conditions": [{"name": "inflammation", "severity": 1}]},
        "available_kingdoms": ["cannabis", "fungal"]
    }
    resp = client.post('/api/dispensary/inflammatory-synergy', json=payload)
    assert resp.status_code == 200
    json_resp = resp.json()
    # Fallback response should include warning about limited accuracy
    assert json_resp.get('warning') is not None


def test_inflammatory_endpoint_cache_consistency():
    payload = {
        "biomarkers": {"tnf_alpha": 10, "il6": 5, "crp": 3, "il1b": 1},
        "condition_profile": {"conditions": [{"name": "inflammation", "severity": 4}]},
        "available_kingdoms": ["cannabis", "fungal", "plant"]
    }
    r1 = client.post('/api/dispensary/inflammatory-synergy', json=payload)
    r2 = client.post('/api/dispensary/inflammatory-synergy', json=payload)
    assert r1.status_code == 200 and r2.status_code == 200
    assert r1.json() == r2.json()  # identical payloads should yield identical cached responses


def test_inflammatory_endpoint_respects_patient_history_preferences():
    payload = {
        "biomarkers": {"tnf_alpha": 8, "il6": 3, "crp": 2},
        "condition_profile": {
            "conditions": [
                {"name": "neuroinflammation", "severity": 8, "is_primary": True}
            ],
            "experience_level": "naive",
            "administration_preferences": ["oral"],
            "primary_goal": "inflammation"
        },
        "available_kingdoms": ["cannabis", "fungal", "plant"]
    }
    resp = client.post('/api/dispensary/inflammatory-synergy', json=payload)
    assert resp.status_code == 200
    data = resp.json()
    dosing = data.get('dosing_guidance', {})
    assert dosing.get('preferred_route') == 'oral'
    notes = ' '.join(dosing.get('notes', []))
    assert 'limited cannabinoid experience' in notes.lower()


def test_inflammatory_endpoint_handles_single_kingdom_high_crp():
    payload = {
        "biomarkers": {"crp": 12, "tnf_alpha": 2},
        "condition_profile": {
            "conditions": [
                {"name": "cardio_inflammation", "severity": 6, "is_primary": True}
            ]
        },
        "available_kingdoms": ["marine"]
    }
    resp = client.post('/api/dispensary/inflammatory-synergy', json=payload)
    assert resp.status_code == 200
    data = resp.json()
    assert data['primary_kingdom'] == 'marine'
    assert data['secondary_kingdoms'] == []
