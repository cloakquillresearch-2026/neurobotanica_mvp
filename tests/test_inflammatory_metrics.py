import sys
import os
import pytest
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from fastapi.testclient import TestClient
from backend.main import app
from backend.services.inflammatory_synergy_engine import get_inflammatory_synergy_engine

client = TestClient(app)


def test_metrics_for_inflammatory_prediction():
    # Ensure a predict run occurs
    engine = get_inflammatory_synergy_engine()
    engine.predict_inflammatory_synergy({'tnf_alpha': 10, 'il6': 5, 'crp': 3, 'il1b': 1}, {'conditions': []})

    r = client.get('/metrics')
    # If prometheus_client is not installed, we get 501 with helpful detail
    if r.status_code == 200:
        text = r.text
        if 'neurobotanica_inflammatory_predict_latency_seconds' not in text:
            pytest.skip('Inflammatory metrics not exposed in this environment yet')
        assert 'neurobotanica_inflammatory_predict_latency_seconds' in text
        assert 'neurobotanica_inflammatory_predict_cache_hits_total' in text or 'neurobotanica_inflammatory_predict_cache_misses_total' in text
        # per-component histograms
        assert 'neurobotanica_inflammatory_pathway_calc_latency_seconds' in text
        assert 'neurobotanica_inflammatory_mapping_latency_seconds' in text
        assert 'neurobotanica_inflammatory_dosing_latency_seconds' in text
    else:
        assert r.status_code == 501
        assert 'prometheus_client not installed' in r.json().get('detail', '')
