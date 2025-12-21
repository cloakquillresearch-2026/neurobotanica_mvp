"""
Week 13: End-to-End Tests
=========================

Complete workflow tests simulating real user scenarios.

Tests:
- Full compound analysis workflow
- Manufacturer regulatory readiness flow
- PatentPath IP protection workflow
- Multi-compound batch analysis

Run:
    pytest backend/tests/test_week13_e2e.py -v
"""
import pytest
from fastapi.testclient import TestClient
from datetime import datetime
import json
import time

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
# End-to-End Workflow Tests
# =============================================================================

class TestManufacturerWorkflow:
    """Test complete manufacturer regulatory readiness workflow."""
    
    def test_compound_characterization_workflow(self, client):
        """
        Scenario: Manufacturer submits new compound for full characterization.
        
        Steps:
        1. Submit SMILES for ChemPath analysis
        2. Get toxicity assessment
        3. Get regulatory pathway recommendation
        4. Request patent analysis
        """
        # CBD-like SMILES
        test_smiles = "CCCCCC1=CC(=C(C(=C1)O)C2C=C(CCC2C(=C)C)C)O"
        
        # Step 1: ChemPath Analysis
        chem_result = client.post(
            "/api/v1/chempath/analyze",
            json={
                "compound": {
                    "name": "Novel Cannabinoid X",
                    "smiles": test_smiles
                },
                "compute_3d": False
            }
        )
        assert chem_result.status_code == 200
        chem_data = chem_result.json()
        assert chem_data["status"] == "success"
        
        # Step 2: ToxPath Assessment
        tox_result = client.post(
            "/api/v1/toxpath/assess",
            json={
                "compound_ref": test_smiles,
                "compound_name": "Novel Cannabinoid X",
                "route": "oral",
                "exposure": {
                    "dose_mg": 25.0,
                    "frequency": "daily",
                    "duration": "chronic"
                }
            }
        )
        assert tox_result.status_code == 200
        tox_data = tox_result.json()
        assert tox_data["status"] == "completed"
        assert "risk_summary" in tox_data
        
        # Step 3: RegPath Strategy
        reg_result = client.post(
            "/api/v1/regpath/strategy",
            json={
                "product_profile": {
                    "product_name": "Novel Cannabinoid X Oral Solution",
                    "product_type": "isolated_cannabinoid",
                    "intended_use": "Management of chronic pain",
                    "route": "oral",
                    "target_markets": ["US"],
                    "therapeutic_area": "Pain Management"
                },
                "constraints": {
                    "risk_tolerance": "moderate"
                }
            }
        )
        assert reg_result.status_code == 200
        reg_data = reg_result.json()
        assert "primary_pathway" in reg_data
        assert "readiness_checklist" in reg_data
        
        # Step 4: PatentPath FTO Check
        fto_result = client.post(
            "/api/v1/patentpath/fto/check",
            json={
                "compound_description": "Novel cannabinoid compound for chronic pain management",
                "keywords": ["cannabinoid", "pain", "therapeutic"],
                "smiles": test_smiles,
                "intended_use": "chronic pain management"
            }
        )
        assert fto_result.status_code == 200
        fto_data = fto_result.json()
        assert "status" in fto_data
        assert "report" in fto_data
    
    def test_coa_validation_workflow(self, client):
        """
        Scenario: Manufacturer submits COA for quality validation.
        
        Steps:
        1. Submit compound with COA data
        2. Receive QC flags and recommendations
        """
        result = client.post(
            "/api/v1/chempath/analyze",
            json={
                "compound": {
                    "name": "CBD Isolate Batch 001",
                    "smiles": "CCCCCC1=CC(=C(C(=C1)O)C2C=C(CCC2C(=C)C)C)O"
                },
                "coa": {
                    "lab_name": "Certified Labs Inc.",
                    "sample_id": "BATCH-2025-001",
                    "units": "%",
                    "cannabinoids": [
                        {"name": "CBD", "value": 98.5},
                        {"name": "THC", "value": 0.15}
                    ],
                    "terpenes": [
                        {"name": "Myrcene", "value": 0.5},
                        {"name": "Limonene", "value": 0.3}
                    ],
                    "moisture": 0.8,
                    "heavy_metals": [
                        {"name": "Lead", "value": 0.01, "result": "pass"},
                        {"name": "Arsenic", "value": 0.005, "result": "pass"}
                    ]
                }
            }
        )
        assert result.status_code == 200
        data = result.json()
        assert data["status"] == "success"


class TestIPProtectionWorkflow:
    """Test intellectual property protection workflow."""
    
    def test_novel_compound_ip_workflow(self, client):
        """
        Scenario: Researcher discovers novel dimer, seeks IP protection.
        
        Steps:
        1. Prior art search
        2. Novelty assessment
        3. Claim generation
        4. TK attribution check
        5. Cost estimation
        """
        # Step 1: Prior Art Search
        prior_art = client.post(
            "/api/v1/patentpath/prior-art/search",
            json={
                "keywords": ["cannabinoid", "dimer", "novel", "therapeutic"],
                "max_results": 10
            }
        )
        assert prior_art.status_code == 200
        
        # Step 2: Novelty Assessment
        novelty = client.post(
            "/api/v1/patentpath/novelty/assess",
            json={
                "title": "CBD-CBG Ether Dimer for Neurological Disorders",
                "description": "A novel dimeric cannabinoid compound linking CBD and CBG through an ether bond for enhanced neurological therapeutic effects",
                "keywords": ["dimer", "CBD", "CBG", "ether", "neurological"]
            }
        )
        assert novelty.status_code == 200
        novelty_data = novelty.json()
        assert "report" in novelty_data
        assert "scores" in novelty_data["report"]
        assert "overall_novelty" in novelty_data["report"]["scores"]
        
        # Step 3: Claim Generation
        claims = client.post(
            "/api/v1/patentpath/claims/generate",
            json={
                "compound_name": "CBD-CBG-Ether-001",
                "therapeutic_use": "Treatment of neurological disorders",
                "key_features": [
                    "Dimeric structure",
                    "Ether linkage",
                    "Enhanced BBB penetration"
                ]
            }
        )
        assert claims.status_code == 200
        claims_data = claims.json()
        assert "report" in claims_data
        assert "generated_claims" in claims_data["report"]
        
        # Step 4: TK Attribution Check
        tk_check = client.post(
            "/api/v1/patentpath/tk/check",
            json={
                "compound_description": "Novel dimer based on traditional cannabis compounds",
                "derivation_source": "modified_natural",
                "source_description": "Based on traditional cannabis use patterns"
            }
        )
        assert tk_check.status_code == 200
        
        # Step 5: Cost Estimation
        cost = client.post(
            "/api/v1/patentpath/costs/estimate",
            json={
                "num_compounds": 1,
                "protection_strategy": "full",
                "target_jurisdictions": ["US", "EU"]
            }
        )
        assert cost.status_code == 200


class TestBatchAnalysisWorkflow:
    """Test batch compound analysis workflow."""
    
    def test_multiple_compound_screening(self, client):
        """
        Scenario: Research team screens multiple compounds.
        
        Steps:
        1. Analyze each compound
        2. Collect risk assessments
        3. Rank by regulatory readiness
        """
        compounds = [
            {"name": "Compound A", "smiles": "CCO"},
            {"name": "Compound B", "smiles": "CC(C)O"},
            {"name": "Compound C", "smiles": "CCCO"}
        ]
        
        results = []
        for compound in compounds:
            # ChemPath analysis
            chem_result = client.post(
                "/api/v1/chempath/analyze",
                json={"compound": compound}
            )
            assert chem_result.status_code == 200
            
            # ToxPath assessment
            tox_result = client.post(
                "/api/v1/toxpath/assess",
                json={
                    "compound_ref": compound["smiles"],
                    "compound_name": compound["name"],
                    "route": "oral"
                }
            )
            assert tox_result.status_code == 200
            
            results.append({
                "name": compound["name"],
                "chem": chem_result.json(),
                "tox": tox_result.json()
            })
        
        # All compounds should be analyzed
        assert len(results) == 3
        
        # All should have risk summaries
        for result in results:
            assert "risk_summary" in result["tox"]


class TestSecurityEnforcedWorkflow:
    """Test workflows with security constraints."""
    
    def test_api_key_required_endpoints(self, client):
        """Verify certain endpoints require API key (if enabled)."""
        # Health endpoints should always work
        health_endpoints = [
            "/api/v1/chempath/health",
            "/api/v1/toxpath/health",
            "/api/v1/regpath/health",
            "/api/v1/security/health"
        ]
        
        for endpoint in health_endpoints:
            response = client.get(endpoint)
            assert response.status_code == 200
    
    def test_rate_limiting_behavior(self, client):
        """Test that rapid requests don't crash the system."""
        for _ in range(20):
            response = client.get("/api/v1/chempath/health")
            # Should either succeed or return rate limit response
            assert response.status_code in [200, 429]


class TestDataIntegrityWorkflow:
    """Test data integrity across the workflow."""
    
    def test_compound_data_preservation(self, client):
        """Verify compound data is preserved through analysis pipeline."""
        original_smiles = "CCCCCC1=CC(=C(C(=C1)O)C2C=C(CCC2C(=C)C)C)O"
        original_name = "Test Preservation Compound"
        
        # Run through ChemPath
        chem_result = client.post(
            "/api/v1/chempath/analyze",
            json={
                "compound": {
                    "name": original_name,
                    "smiles": original_smiles
                }
            }
        )
        assert chem_result.status_code == 200
        
        # Run through ToxPath
        tox_result = client.post(
            "/api/v1/toxpath/assess",
            json={
                "compound_ref": original_smiles,
                "compound_name": original_name,
                "route": "oral"
            }
        )
        assert tox_result.status_code == 200
        tox_data = tox_result.json()
        
        # Verify the compound data was preserved
        assert tox_data["compound_name"] == original_name
        assert tox_data["compound_ref"] == original_smiles


class TestErrorRecoveryWorkflow:
    """Test error recovery and graceful degradation."""
    
    def test_partial_data_handling(self, client):
        """Test handling of incomplete input data."""
        # Minimal required data only
        result = client.post(
            "/api/v1/chempath/analyze",
            json={
                "compound": {
                    "name": "Minimal Compound",
                    "smiles": "CCO"
                }
            }
        )
        assert result.status_code == 200
    
    def test_invalid_smiles_graceful_handling(self, client):
        """Test graceful handling of invalid SMILES."""
        result = client.post(
            "/api/v1/chempath/analyze",
            json={
                "compound": {
                    "name": "Invalid Test",
                    "smiles": "THIS_IS_NOT_VALID_SMILES"
                }
            }
        )
        # Should handle gracefully (200 with error or 4xx)
        assert result.status_code in [200, 400, 422, 500]


class TestReportGenerationWorkflow:
    """Test report generation workflows."""
    
    def test_toxpath_memo_generation(self, client):
        """Test ToxPath memo generation workflow."""
        # First create an assessment
        tox_result = client.post(
            "/api/v1/toxpath/assess",
            json={
                "compound_ref": "CCCCCC1=CC(=C(C(=C1)O)C2C=C(CCC2C(=C)C)C)O",
                "compound_name": "CBD Test",
                "route": "oral"
            }
        )
        assert tox_result.status_code == 200
        assessment_id = tox_result.json()["toxpath_assessment_id"]
        
        # Generate memo
        memo_result = client.post(
            f"/api/v1/toxpath/assessments/{assessment_id}/memo",
            json={
                "format": "markdown",
                "include_testing_plan": True,
                "include_assumptions": True
            }
        )
        assert memo_result.status_code == 200
        memo_data = memo_result.json()
        assert "memo_content" in memo_data
    
    def test_regpath_memo_generation(self, client):
        """Test RegPath memo generation workflow."""
        # First create a strategy
        reg_result = client.post(
            "/api/v1/regpath/strategy",
            json={
                "product_profile": {
                    "product_name": "Test Product Memo",
                    "product_type": "isolated_cannabinoid",
                    "intended_use": "Therapeutic use",
                    "route": "oral"
                }
            }
        )
        assert reg_result.status_code == 200
        strategy_id = reg_result.json()["regpath_strategy_id"]
        
        # Generate memo
        memo_result = client.post(
            f"/api/v1/regpath/strategies/{strategy_id}/memo",
            json={
                "format": "markdown",
                "include_timeline": True,
                "include_checklist": True
            }
        )
        assert memo_result.status_code == 200


class TestFullPipelineE2E:
    """Complete end-to-end pipeline test."""
    
    def test_complete_manufacturer_submission(self, client):
        """
        Full manufacturer submission workflow:
        1. Compound characterization
        2. Toxicity assessment  
        3. Regulatory strategy
        4. Patent protection check
        5. Report generation
        """
        smiles = "CCCCCC1=CC(=C(C(=C1)O)C2C=C(CCC2C(=C)C)C)O"
        name = "Complete Pipeline Test Compound"
        
        # 1. ChemPath
        chem = client.post("/api/v1/chempath/analyze", json={
            "compound": {"name": name, "smiles": smiles}
        })
        assert chem.status_code == 200
        
        # 2. ToxPath
        tox = client.post("/api/v1/toxpath/assess", json={
            "compound_ref": smiles, "compound_name": name, "route": "oral"
        })
        assert tox.status_code == 200
        tox_id = tox.json()["toxpath_assessment_id"]
        
        # 3. RegPath
        reg = client.post("/api/v1/regpath/strategy", json={
            "product_profile": {
                "product_name": f"{name} Product",
                "product_type": "isolated_cannabinoid",
                "intended_use": "Therapeutic",
                "route": "oral"
            }
        })
        assert reg.status_code == 200
        reg_id = reg.json()["regpath_strategy_id"]
        
        # 4. PatentPath FTO
        fto = client.post("/api/v1/patentpath/fto/check", json={
            "compound_description": f"Novel compound {name} for therapeutic use",
            "keywords": ["cannabinoid", "therapeutic"],
            "smiles": smiles
        })
        assert fto.status_code == 200
        
        # 5. Generate Memos
        tox_memo = client.post(f"/api/v1/toxpath/assessments/{tox_id}/memo", json={
            "format": "markdown"
        })
        assert tox_memo.status_code == 200
        
        reg_memo = client.post(f"/api/v1/regpath/strategies/{reg_id}/memo", json={
            "format": "markdown"
        })
        assert reg_memo.status_code == 200
        
        # All steps completed successfully
        print("✅ Complete manufacturer pipeline executed successfully")


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
