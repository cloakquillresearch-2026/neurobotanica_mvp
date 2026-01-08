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
