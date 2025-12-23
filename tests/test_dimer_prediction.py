import pytest

from fastapi import FastAPI
from fastapi.testclient import TestClient

import sys
from unittest.mock import MagicMock

# Patch backend.database.get_db before importing router
sys.modules['backend.database'] = MagicMock()
from backend.api import dimers

app = FastAPI()
app.include_router(dimers.router)
client = TestClient(app)

def test_predict_homodimer():
    payload = {
        "monomer": {
            "name": "THC",
            "smiles": "CCCCCc1cc(O)c2C3C=C(C)CCC3C(C)(C)Oc2c1"
        },
        "linkage_type": "METHYLENE"
    }
    response = client.post("/api/dimers/predict/homodimer", json=payload)
    assert response.status_code == 200
    data = response.json()
    assert "dimer_name" in data
    assert "synergy_score" in data
    assert "predicted_properties" in data
    assert isinstance(data["synergy_score"], float)

def test_predict_heterodimer():
    payload = {
        "monomer_a": {
            "name": "THC",
            "smiles": "CCCCCc1cc(O)c2C3C=C(C)CCC3C(C)(C)Oc2c1"
        },
        "monomer_b": {
            "name": "CBD",
            "smiles": "CCCCCc1cc(O)c(C2C=C(C)CCC2C(C)=C)c(O)c1"
        },
        "linkage_type": "METHYLENE"
    }
    response = client.post("/api/dimers/predict/heterodimer", json=payload)
    assert response.status_code == 200
    data = response.json()
    assert "dimer_name" in data
    assert "synergy_score" in data
    assert "predicted_properties" in data
    assert isinstance(data["synergy_score"], float)
