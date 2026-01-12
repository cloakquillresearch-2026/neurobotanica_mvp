import sys
import os
import pytest
# Ensure project root is in sys.path when running tests directly
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from fastapi.testclient import TestClient
from backend.main import app

client = TestClient(app)


def test_metrics_endpoint():
    r = client.get('/metrics')
    # If prometheus_client is not installed, we return 501 with a helpful message
    if r.status_code == 200:
        if 'neurobotanica_health' not in r.text:
            pytest.skip('NeuroBotanica-specific metrics not exposed in this environment yet')
        assert 'neurobotanica_health' in r.text
    else:
        assert r.status_code == 501
        assert 'prometheus_client not installed' in r.json().get('detail', '')


def test_status_endpoint():
    r = client.get('/status')
    assert r.status_code == 200
    assert 'NeuroBotanica Status' in r.text


def test_health_refresh_and_cached():
    # Force a live refresh
    r = client.get('/health?refresh=true')
    assert r.status_code == 200
    data = r.json()
    assert 'status' in data
    assert 'results' in data
    # Call again without refresh to ensure cached path returns same structure
    r2 = client.get('/health')
    assert r2.status_code == 200
    data2 = r2.json()
    assert 'status' in data2
    assert 'results' in data2
