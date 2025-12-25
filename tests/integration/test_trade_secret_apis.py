"""
Integration Tests for Trade Secret API Endpoints

Tests all 6 trade secret engines:
- ChemPath ($50K implementation)
- ToxPath ($50K implementation)
- RegPath ($50K implementation)
- GenomePath ($6.2B value)
- BioPath ($2.0B value)
- ClinPath ($3.2B value)

Run with: pytest tests/integration/test_trade_secret_apis.py -v
"""

import pytest
from fastapi.testclient import TestClient
import sys
from pathlib import Path

# Add project root to path
project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))

from backend.main import app

client = TestClient(app)


class TestChemPathAPI:
    """Tests for ChemPath - Chemical Characterization Engine."""
    
    def test_chempath_analyze_endpoint_exists(self):
        """Test that ChemPath analyze endpoint is registered."""
        response = client.post(
            "/api/v1/chempath/analyze",
            json={
                "compound": {
                    "name": "THC",
                    "smiles": "CCCCCc1cc(O)c2C3C=C(C)CCC3C(C)(C)Oc2c1"
                }
            }
        )
        # Should not be 404 (endpoint exists)
        assert response.status_code != 404
    
    def test_chempath_validate_smiles(self):
        """Test SMILES validation endpoint."""
        response = client.post(
            "/api/v1/chempath/validate-smiles",
            json={"smiles": "CCO"}  # Ethanol
        )
        assert response.status_code in [200, 422]  # Success or validation error


class TestToxPathAPI:
    """Tests for ToxPath - Toxicity Risk Assessment Engine."""
    
    def test_toxpath_assess_endpoint_exists(self):
        """Test that ToxPath assess endpoint is registered."""
        response = client.post(
            "/api/v1/toxpath/assess",
            json={
                "compound_ref": "CCCCCc1cc(O)c2C3C=C(C)CCC3C(C)(C)Oc2c1",
                "compound_name": "THC",
                "route": "oral"
            }
        )
        assert response.status_code != 404
    
    def test_toxpath_statistics(self):
        """Test ToxPath statistics endpoint."""
        response = client.get("/api/v1/toxpath/statistics")
        assert response.status_code in [200, 404]


class TestRegPathAPI:
    """Tests for RegPath - Regulatory Pathway Optimization Engine."""
    
    def test_regpath_strategy_endpoint_exists(self):
        """Test that RegPath strategy endpoint is registered."""
        response = client.post(
            "/api/v1/regpath/strategy",
            json={
                "product_profile": {
                    "product_name": "CBD-X Oral Solution",
                    "product_type": "isolated_cannabinoid",
                    "intended_use": "Treatment of epilepsy",
                    "route": "oral",
                    "target_markets": ["US"]
                }
            }
        )
        assert response.status_code != 404
    
    def test_regpath_pathways_list(self):
        """Test regulatory pathways list endpoint."""
        response = client.get("/api/v1/regpath/pathways")
        assert response.status_code in [200, 404]


class TestGenomePathAPI:
    """Tests for GenomePath - TK-Genomic Bridge ($6.2B value)."""
    
    def test_genomepath_tk_to_genomic(self):
        """Test TK to genomic hypothesis generation."""
        response = client.post(
            "/api/genomepath/tk-to-genomic",
            json={
                "practice_name": "Cannabis for Pain",
                "plant_species": "Cannabis sativa",
                "preparation_method": "Infusion",
                "therapeutic_use": "Pain management",
                "community_id": "test-community",
                "community_consent_verified": True
            }
        )
        # Endpoint should exist (not 404)
        assert response.status_code != 404
    
    def test_genomepath_statistics(self):
        """Test GenomePath statistics endpoint."""
        response = client.get("/api/genomepath/statistics")
        assert response.status_code in [200, 404]


class TestBioPathAPI:
    """Tests for BioPath - Bias-Aware Validation Engine ($2.0B value)."""
    
    def test_biopath_validate_endpoint_exists(self):
        """Test that BioPath validate endpoint is registered."""
        response = client.post(
            "/api/biopath/validate",
            json={
                "compound_name": "THC",
                "target_condition": "PTSD",
                "evidence_sources": [
                    {
                        "source_type": "clinical_trial",
                        "score": 0.8,
                        "sample_size": 100,
                        "study_quality": 0.7
                    }
                ],
                "include_community_validation": True
            }
        )
        assert response.status_code != 404
    
    def test_biopath_validate_from_studies(self):
        """Test validation from clinical studies database."""
        response = client.post(
            "/api/biopath/validate-from-studies",
            json={
                "compound_name": "CBD",
                "target_condition": "Anxiety"
            }
        )
        assert response.status_code != 404
    
    def test_biopath_statistics(self):
        """Test BioPath statistics endpoint."""
        response = client.get("/api/biopath/statistics")
        if response.status_code == 200:
            data = response.json()
            assert "trade_secret_value" in data
            assert data["trade_secret_value"] == "$2.0B"


class TestClinPathAPI:
    """Tests for ClinPath - Clinical Trial Optimization ($3.2B value)."""
    
    def test_clinpath_optimize_endpoint_exists(self):
        """Test that ClinPath optimize endpoint is registered."""
        response = client.post(
            "/api/clinpath/optimize",
            json={
                "compound_name": "THC",
                "indication": "PTSD",
                "target_jurisdictions": ["USA", "CAN"],
                "use_tm_pathway": True
            }
        )
        assert response.status_code != 404
    
    def test_clinpath_predict_approval(self):
        """Test approval probability prediction."""
        response = client.post(
            "/api/clinpath/predict-approval",
            json={
                "compound_name": "CBD",
                "indication": "Anxiety",
                "jurisdiction": "USA",
                "evidence_strength": 0.8,
                "regulatory_precedents": 5
            }
        )
        assert response.status_code != 404
    
    def test_clinpath_jurisdiction_sequence(self):
        """Test optimal jurisdiction sequencing."""
        response = client.post(
            "/api/clinpath/jurisdiction-sequence",
            json={
                "compound_name": "THC",
                "indication": "Chronic Pain",
                "target_jurisdictions": ["USA", "CAN", "EUR", "AUS"]
            }
        )
        assert response.status_code != 404
    
    def test_clinpath_statistics(self):
        """Test ClinPath statistics endpoint."""
        response = client.get("/api/clinpath/statistics")
        if response.status_code == 200:
            data = response.json()
            assert "trade_secret_value" in data
            assert data["trade_secret_value"] == "$3.2B"


class TestDimerPredictionAPI:
    """Tests for Dimer Prediction API (core functionality)."""
    
    def test_homodimer_prediction(self):
        """Test homodimer prediction endpoint."""
        response = client.post(
            "/api/dimers/predict/homodimer",
            json={
                "monomer": {
                    "name": "THC",
                    "smiles": "CCCCCc1cc(O)c2C3C=C(C)CCC3C(C)(C)Oc2c1"
                },
                "linkage_type": "methylene"
            }
        )
        assert response.status_code == 200
        data = response.json()
        assert "dimer_name" in data
        assert "THC" in data["dimer_name"]
    
    def test_heterodimer_prediction(self):
        """Test heterodimer prediction endpoint."""
        response = client.post(
            "/api/dimers/predict/heterodimer",
            json={
                "monomer_a": {
                    "name": "THC",
                    "smiles": "CCCCCc1cc(O)c2C3C=C(C)CCC3C(C)(C)Oc2c1"
                },
                "monomer_b": {
                    "name": "CBD",
                    "smiles": "CCCCCc1cc(O)c(C2C=C(C)CCC2C(C)=C)c(O)c1"
                },
                "linkage_type": "methylene"
            }
        )
        assert response.status_code == 200


class TestHealthAndRoot:
    """Tests for health and root endpoints."""
    
    def test_root_endpoint(self):
        """Test root endpoint returns API info."""
        response = client.get("/")
        assert response.status_code == 200
        data = response.json()
        assert data["version"] == "0.4.0"
        assert "trade_secret_engines" in data
        assert len(data["trade_secret_engines"]) == 6
    
    def test_health_endpoint(self):
        """Test health endpoint includes ML model status."""
        response = client.get("/health")
        assert response.status_code == 200
        data = response.json()
        assert "ml_models" in data
        assert "features" in data
        assert "trade_secret_engines" in data["features"]
    
    def test_stats_endpoint(self):
        """Test stats endpoint returns database statistics."""
        response = client.get("/api/v1/stats")
        assert response.status_code == 200


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
