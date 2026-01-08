from fastapi.testclient import TestClient
from backend.main import app

client = TestClient(app)

payload = {
    "biomarkers": {"tnf_alpha": 20, "il6": 10, "crp": 5, "il1b": 2},
    "condition_profile": {"conditions": [{"name": "inflammation", "severity": 7, "is_primary": True}]},
    "available_kingdoms": ["cannabis", "fungal", "plant"]
}

resp = client.post('/api/dispensary/inflammatory-synergy', json=payload)
print('Status', resp.status_code)
print(resp.json())
