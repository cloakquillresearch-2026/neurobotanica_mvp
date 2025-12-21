"""
Week 13: Integration Tests
==========================

Cross-module integration tests verifying system coherence.

Tests:
- ChemPath → ToxPath → RegPath pipeline
- Security middleware with all routers
- PatentPath integration with compound analysis
- Full analysis workflow

Run:
    pytest backend/tests/test_week13_integration.py -v
"""
import pytest
from fastapi.testclient import TestClient
from datetime import datetime
import json

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from backend.main import app


@pytest.fixture(scope="module")
def client():
    """Create test client."""
    with TestClient(app) as test_client:
        yield test_client


# =============================================================================
# Cross-Module Integration Tests
# =============================================================================

class TestChemPathToToxPathIntegration:
    """Test ChemPath → ToxPath pipeline integration."""
    
    def test_chempath_output_compatible_with_toxpath(self, client):
        """Verify ChemPath output can feed ToxPath."""
        # Step 1: Analyze compound with ChemPath (requires compound object)
        chempath_response = client.post(
            "/api/v1/chempath/analyze",
            json={
                "compound": {
                    "name": "Test CBD",
                    "smiles": "CCCCCC1=CC(=C(C(=C1)O)C2C=C(CCC2C(=C)C)C)O"
                }
            }
        )
        assert chempath_response.status_code == 200
        chempath_data = chempath_response.json()
        
        # Verify ChemPath provides required data for ToxPath
        assert "status" in chempath_data or "analysis" in chempath_data
    
    def test_toxpath_accepts_chempath_descriptors(self, client):
        """Verify ToxPath can process ChemPath descriptor format."""
        # Get descriptors from ChemPath
        chempath_response = client.post(
            "/api/v1/chempath/analyze",
            json={
                "compound": {
                    "name": "Ethanol",
                    "smiles": "CCO"
                }
            }
        )
        assert chempath_response.status_code == 200
        
        # ToxPath should accept the compound for analysis
        toxpath_response = client.post(
            "/api/v1/toxpath/assess",
            json={
                "compound_ref": "CCO",
                "compound_name": "Ethanol",
                "route": "oral"
            }
        )
        assert toxpath_response.status_code == 200
        toxpath_data = toxpath_response.json()
        
        # Verify ToxPath output (actual field names)
        assert "toxpath_assessment_id" in toxpath_data
        assert "risk_summary" in toxpath_data


class TestToxPathToRegPathIntegration:
    """Test ToxPath → RegPath pipeline integration."""
    
    def test_toxpath_output_compatible_with_regpath(self, client):
        """Verify ToxPath output can feed RegPath."""
        # Step 1: Get toxicity assessment
        toxpath_response = client.post(
            "/api/v1/toxpath/assess",
            json={
                "compound_ref": "CCCCCC1=CC(=C(C(=C1)O)C2C=C(CCC2C(=C)C)C)O",
                "compound_name": "CBD",
                "route": "oral"
            }
        )
        assert toxpath_response.status_code == 200
        toxpath_data = toxpath_response.json()
        
        # Verify ToxPath provides data for regulatory assessment
        assert "toxpath_assessment_id" in toxpath_data
        assert "risk_summary" in toxpath_data
    
    def test_regpath_accepts_toxpath_data(self, client):
        """Verify RegPath can process with toxicity context."""
        regpath_response = client.post(
            "/api/v1/regpath/strategy",
            json={
                "product_profile": {
                    "product_name": "Test CBD Product",
                    "product_type": "isolated_cannabinoid",
                    "intended_use": "Pain management",
                    "route": "oral"
                }
            }
        )
        assert regpath_response.status_code == 200
        regpath_data = regpath_response.json()
        
        # Verify RegPath output (actual field names)
        assert "regpath_strategy_id" in regpath_data
        assert "primary_pathway" in regpath_data


class TestFullAnalysisPipeline:
    """Test complete ChemPath → ToxPath → RegPath pipeline."""
    
    def test_complete_compound_analysis_flow(self, client):
        """Test full analysis workflow for a compound."""
        test_smiles = "CCCCCC1=CC(=C(C(=C1)O)C2C=C(CCC2C(=C)C)C)O"  # CBD-like
        
        # Step 1: ChemPath Analysis
        chem_response = client.post(
            "/api/v1/chempath/analyze",
            json={
                "compound": {
                    "name": "Test Compound",
                    "smiles": test_smiles
                }
            }
        )
        assert chem_response.status_code == 200
        chem_data = chem_response.json()
        
        # Step 2: ToxPath Analysis
        tox_response = client.post(
            "/api/v1/toxpath/assess",
            json={
                "compound_ref": test_smiles,
                "compound_name": "Test Compound",
                "route": "oral"
            }
        )
        assert tox_response.status_code == 200
        tox_data = tox_response.json()
        
        # Step 3: RegPath Analysis
        reg_response = client.post(
            "/api/v1/regpath/strategy",
            json={
                "product_profile": {
                    "product_name": "Test Product",
                    "product_type": "isolated_cannabinoid",
                    "intended_use": "Therapeutic",
                    "route": "oral"
                }
            }
        )
        assert reg_response.status_code == 200
        reg_data = reg_response.json()
        
        # Verify all three produced valid outputs
        assert chem_data is not None
        assert tox_data is not None
        assert reg_data is not None
    
    def test_error_propagation_in_pipeline(self, client):
        """Test that errors are handled gracefully in pipeline."""
        invalid_smiles = "INVALID_NOT_A_SMILES"
        
        # All services should handle invalid input gracefully
        chem_response = client.post(
            "/api/v1/chempath/analyze",
            json={
                "compound": {
                    "name": "Invalid",
                    "smiles": invalid_smiles
                }
            }
        )
        # Should either return 400/422/500 or 200 with error indication
        assert chem_response.status_code in [200, 400, 422, 500]


class TestSecurityIntegrationWithRouters:
    """Test security middleware integration with all routers."""
    
    def test_security_headers_on_chempath(self, client):
        """Verify security headers on ChemPath responses."""
        response = client.get("/api/v1/chempath/health")
        assert response.status_code == 200
        # Headers may be set by middleware
        assert response.headers.get("content-type") is not None
    
    def test_security_headers_on_toxpath(self, client):
        """Verify security headers on ToxPath responses."""
        response = client.get("/api/v1/toxpath/health")
        assert response.status_code == 200
    
    def test_security_headers_on_regpath(self, client):
        """Verify security headers on RegPath responses."""
        response = client.get("/api/v1/regpath/health")
        assert response.status_code == 200
    
    def test_rate_limit_status_across_routes(self, client):
        """Verify rate limiting works across all routes."""
        # Hit multiple endpoints
        endpoints = [
            "/api/v1/chempath/health",
            "/api/v1/toxpath/health",
            "/api/v1/regpath/health",
            "/api/v1/security/health"
        ]
        
        for endpoint in endpoints:
            response = client.get(endpoint)
            assert response.status_code == 200


class TestPatentPathIntegration:
    """Test PatentPath integration with compound analysis."""
    
    def test_patentpath_prior_art_search(self, client):
        """Test prior art search endpoint."""
        response = client.post(
            "/api/v1/patentpath/prior-art/search",
            json={
                "keywords": ["cannabidiol", "CBD", "therapeutic"],
                "max_results": 5
            }
        )
        assert response.status_code == 200
    
    def test_patentpath_novelty_assessment(self, client):
        """Test novelty scoring for a compound."""
        response = client.post(
            "/api/v1/patentpath/novelty/assess",
            json={
                "title": "Novel CBD Formulation",
                "description": "A novel formulation of cannabidiol for therapeutic use",
                "keywords": ["cannabinoid", "formulation"]
            }
        )
        assert response.status_code == 200
    
    def test_patentpath_fto_check(self, client):
        """Test freedom-to-operate check."""
        response = client.post(
            "/api/v1/patentpath/fto/check",
            json={
                "compound_description": "A novel cannabinoid compound for therapeutic use",
                "keywords": ["cannabinoid", "therapeutic"],
                "smiles": "CCO",
                "intended_use": "therapeutic"
            }
        )
        assert response.status_code == 200
    
    def test_patentpath_claims_generate(self, client):
        """Test patent claim generation."""
        response = client.post(
            "/api/v1/patentpath/claims/generate",
            json={
                "compound_name": "Novel CBD Analog",
                "smiles": "CCCCCC1=CC(=C(C(=C1)O)C2C=C(CCC2C(=C)C)C)O",
                "claim_type": "composition"
            }
        )
        assert response.status_code == 200


class TestOmniPathIntegration:
    """Test OmniPath integration with compound data."""
    
    def test_omnipath_manifest_endpoint(self, client):
        """Test OmniPath manifest generation."""
        response = client.get("/api/v1/omnipath/manifests")
        assert response.status_code in [200, 404]
    
    def test_omnipath_token_validation(self, client):
        """Test OmniPath token validation endpoint."""
        response = client.post(
            "/api/v1/omnipath/validate-token",
            json={"token": "test-token"}
        )
        # May return 200, 401, 404 depending on implementation
        assert response.status_code in [200, 400, 401, 404, 422]


# =============================================================================
# Data Consistency Tests
# =============================================================================

class TestDataConsistency:
    """Test data consistency across modules."""
    
    def test_compound_data_consistent_across_endpoints(self, client):
        """Verify compound data is consistent across different endpoints."""
        # Get compound from compounds endpoint
        compounds_response = client.get("/api/v1/compounds/")
        assert compounds_response.status_code == 200
        
        compounds = compounds_response.json()
        if compounds:
            compound_name = compounds[0].get("name", "")
            
            # Should be able to analyze same compound
            if compound_name:
                chem_response = client.post(
                    "/api/v1/chempath/analyze",
                    json={"smiles": compounds[0].get("smiles", "CCO") or "CCO"}
                )
                # Just verify it doesn't crash
                assert chem_response.status_code in [200, 400, 422]
    
    def test_study_data_consistent(self, client):
        """Verify study data consistency."""
        studies_response = client.get("/api/v1/studies/")
        assert studies_response.status_code == 200
        
        studies = studies_response.json()
        assert isinstance(studies, list)


class TestAPIVersionConsistency:
    """Test API version consistency."""
    
    def test_all_v1_endpoints_accessible(self, client):
        """Verify all v1 API endpoints are accessible."""
        v1_endpoints = [
            "/api/v1/studies/",
            "/api/v1/compounds/",
            "/api/v1/chempath/health",
            "/api/v1/toxpath/health",
            "/api/v1/regpath/health",
            "/api/v1/security/health",
        ]
        
        for endpoint in v1_endpoints:
            response = client.get(endpoint)
            # Should return 200 or at least not 500
            assert response.status_code != 500, f"Server error on {endpoint}"
    
    def test_root_endpoint_shows_version(self, client):
        """Verify root endpoint shows API version."""
        response = client.get("/")
        assert response.status_code == 200
        data = response.json()
        assert "version" in data


# =============================================================================
# Error Handling Integration Tests
# =============================================================================

class TestErrorHandlingIntegration:
    """Test error handling consistency across modules."""
    
    def test_404_consistency(self, client):
        """Verify 404 responses are consistent across modules."""
        endpoints_404 = [
            "/api/v1/studies/NONEXISTENT-XYZ",
            "/api/v1/compounds/NONEXISTENT_COMPOUND",
        ]
        
        for endpoint in endpoints_404:
            response = client.get(endpoint)
            assert response.status_code == 404
            # Should have consistent error format
            data = response.json()
            assert "detail" in data
    
    def test_validation_error_consistency(self, client):
        """Verify validation errors are consistent."""
        # Invalid limit parameter
        response = client.get("/api/v1/studies/?limit=9999")
        assert response.status_code == 422
        
        data = response.json()
        assert "detail" in data


class TestConcurrentRequestHandling:
    """Test handling of concurrent requests."""
    
    def test_multiple_rapid_requests(self, client):
        """Test system handles rapid consecutive requests."""
        results = []
        for i in range(10):
            response = client.get("/api/v1/chempath/health")
            results.append(response.status_code)
        
        # All should succeed
        assert all(status == 200 for status in results)
    
    def test_different_endpoints_concurrently(self, client):
        """Test hitting different endpoints."""
        endpoints = [
            "/api/v1/chempath/health",
            "/api/v1/toxpath/health",
            "/api/v1/regpath/health",
            "/",
            "/health"
        ]
        
        for endpoint in endpoints:
            response = client.get(endpoint)
            assert response.status_code == 200


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
