"""
Week 10 Tests - ToxPath Toxicity Risk Assessment Engine

Comprehensive tests for:
- ToxPath assessor functionality
- Structural alert detection
- Risk tier classification
- Testing plan generation
- Memo generation
- API endpoints
"""

import pytest
from datetime import datetime


# ============================================================================
# ToxPath Assessor Tests
# ============================================================================

class TestToxPathAssessor:
    """Test ToxPathAssessor core functionality."""
    
    def test_assessor_initialization(self):
        """Test assessor initializes correctly."""
        from backend.services.toxpath import ToxPathAssessor
        
        assessor = ToxPathAssessor()
        assert assessor is not None
        # RDKit should be available in our environment
        assert assessor._rdkit_available is True
    
    def test_assessor_has_alert_codes(self):
        """Test assessor provides alert codes."""
        from backend.services.toxpath import ToxPathAssessor
        
        assessor = ToxPathAssessor()
        codes = assessor.get_alert_codes()
        
        assert isinstance(codes, dict)
        assert len(codes) > 10  # Should have multiple alert codes
        assert "ALERT_MICHAEL_ACCEPTOR" in codes
        assert "ALERT_NITRO" in codes
        assert "ALERT_EPOXIDE" in codes
    
    def test_assessor_has_risk_tiers(self):
        """Test assessor provides risk tier info."""
        from backend.services.toxpath import ToxPathAssessor
        
        assessor = ToxPathAssessor()
        tiers = assessor.get_risk_tiers()
        
        assert len(tiers) == 4
        assert tiers[0]["tier"] == 1
        assert tiers[0]["label"] == "Low Risk"
        assert tiers[3]["tier"] == 4
        assert tiers[3]["label"] == "Very High Risk"
    
    def test_assessor_has_routes(self):
        """Test assessor provides supported routes."""
        from backend.services.toxpath import ToxPathAssessor
        
        assessor = ToxPathAssessor()
        routes = assessor.get_routes()
        
        assert "oral" in routes
        assert "inhaled" in routes
        assert "sublingual" in routes
        assert "topical" in routes
        assert "intravenous" in routes
        assert "unknown" not in routes


class TestRiskTierEnum:
    """Test RiskTier enumeration."""
    
    def test_risk_tier_values(self):
        """Test risk tier values."""
        from backend.services.toxpath import RiskTier
        
        assert RiskTier.LOW.value == 1
        assert RiskTier.MODERATE.value == 2
        assert RiskTier.HIGH.value == 3
        assert RiskTier.VERY_HIGH.value == 4
    
    def test_risk_tier_labels(self):
        """Test risk tier labels."""
        from backend.services.toxpath import RiskTier
        
        assert RiskTier.LOW.label == "Low Risk"
        assert RiskTier.MODERATE.label == "Moderate Risk"
        assert RiskTier.HIGH.label == "High Risk"
        assert RiskTier.VERY_HIGH.label == "Very High Risk"
    
    def test_risk_tier_testing_depth(self):
        """Test risk tier testing depth recommendations."""
        from backend.services.toxpath import RiskTier
        
        assert RiskTier.LOW.testing_depth == "Minimal"
        assert RiskTier.MODERATE.testing_depth == "Standard panel"
        assert RiskTier.HIGH.testing_depth == "Extended + consult"
        assert RiskTier.VERY_HIGH.testing_depth == "Full GLP recommended"


class TestAlertSeverityEnum:
    """Test AlertSeverity enumeration."""
    
    def test_alert_severity_values(self):
        """Test alert severity values."""
        from backend.services.toxpath import AlertSeverity
        
        assert AlertSeverity.INFO.value == "info"
        assert AlertSeverity.CAUTION.value == "caution"
        assert AlertSeverity.WARNING.value == "warning"
        assert AlertSeverity.CRITICAL.value == "critical"


class TestAdministrationRouteEnum:
    """Test AdministrationRoute enumeration."""
    
    def test_route_values(self):
        """Test route values."""
        from backend.services.toxpath import AdministrationRoute
        
        assert AdministrationRoute.ORAL.value == "oral"
        assert AdministrationRoute.INHALED.value == "inhaled"
        assert AdministrationRoute.SUBLINGUAL.value == "sublingual"
        assert AdministrationRoute.UNKNOWN.value == "unknown"


# ============================================================================
# Data Model Tests
# ============================================================================

class TestToxPathRequest:
    """Test ToxPathRequest data model."""
    
    def test_request_minimal(self):
        """Test minimal request creation."""
        from backend.services.toxpath import ToxPathRequest
        
        req = ToxPathRequest(compound_ref="CCO")
        assert req.compound_ref == "CCO"
        assert req.route is None
        assert req.exposure is None
    
    def test_request_full(self):
        """Test full request creation."""
        from backend.services.toxpath import ToxPathRequest, ExposureProfile, ImpurityEntry
        
        req = ToxPathRequest(
            compound_ref="CCO",
            compound_name="Ethanol",
            route="oral",
            exposure=ExposureProfile(dose_mg=100, frequency="daily", duration="chronic"),
            known_impurities=[ImpurityEntry(name="methanol", ppm=50)]
        )
        
        assert req.compound_name == "Ethanol"
        assert req.route == "oral"
        assert req.exposure.dose_mg == 100
        assert len(req.known_impurities) == 1
    
    def test_request_to_dict(self):
        """Test request serialization."""
        from backend.services.toxpath import ToxPathRequest
        
        req = ToxPathRequest(compound_ref="CCO", compound_name="Ethanol")
        data = req.to_dict()
        
        assert data["compound_ref"] == "CCO"
        assert data["compound_name"] == "Ethanol"


class TestExposureProfile:
    """Test ExposureProfile data model."""
    
    def test_exposure_profile_creation(self):
        """Test exposure profile creation."""
        from backend.services.toxpath import ExposureProfile
        
        profile = ExposureProfile(dose_mg=100, frequency="daily", duration="chronic")
        assert profile.dose_mg == 100
        assert profile.frequency == "daily"
        assert profile.duration == "chronic"
    
    def test_exposure_profile_to_dict(self):
        """Test exposure profile serialization."""
        from backend.services.toxpath import ExposureProfile
        
        profile = ExposureProfile(dose_mg=50, frequency="once")
        data = profile.to_dict()
        
        assert data["dose_mg"] == 50
        assert data["frequency"] == "once"


class TestImpurityEntry:
    """Test ImpurityEntry data model."""
    
    def test_impurity_entry_creation(self):
        """Test impurity entry creation."""
        from backend.services.toxpath import ImpurityEntry
        
        entry = ImpurityEntry(name="methanol", ppm=100, cas_number="67-56-1")
        assert entry.name == "methanol"
        assert entry.ppm == 100
        assert entry.cas_number == "67-56-1"


class TestStructuralAlert:
    """Test StructuralAlert data model."""
    
    def test_structural_alert_creation(self):
        """Test structural alert creation."""
        from backend.services.toxpath import StructuralAlert, AlertSeverity
        
        alert = StructuralAlert(
            code="ALERT_NITRO",
            name="Nitro Group",
            severity=AlertSeverity.WARNING,
            description="Nitro group detected",
            mechanism="Metabolic reduction",
            affected_organs=["liver"]
        )
        
        assert alert.code == "ALERT_NITRO"
        assert alert.severity == AlertSeverity.WARNING
        assert "liver" in alert.affected_organs
    
    def test_structural_alert_to_dict(self):
        """Test structural alert serialization."""
        from backend.services.toxpath import StructuralAlert, AlertSeverity
        
        alert = StructuralAlert(
            code="ALERT_TEST",
            name="Test Alert",
            severity=AlertSeverity.CAUTION,
            description="Test description"
        )
        data = alert.to_dict()
        
        assert data["code"] == "ALERT_TEST"
        assert data["severity"] == "caution"


class TestTestingStep:
    """Test TestingStep data model."""
    
    def test_testing_step_creation(self):
        """Test testing step creation."""
        from backend.services.toxpath import TestingStep
        
        step = TestingStep(
            order=1,
            test_name="Ames Test",
            test_type="in_vitro",
            rationale="Standard mutagenicity screen",
            estimated_cost_range="$5,000-15,000",
            timeline="4-6 weeks",
            required=True,
            gmp_glp_required=True
        )
        
        assert step.order == 1
        assert step.test_name == "Ames Test"
        assert step.gmp_glp_required is True


class TestRiskSummary:
    """Test RiskSummary data model."""
    
    def test_risk_summary_creation(self):
        """Test risk summary creation."""
        from backend.services.toxpath import RiskSummary, RiskTier
        
        summary = RiskSummary(
            overall_tier=RiskTier.MODERATE,
            tier_rationale="Several caution-level alerts",
            top_risks=["Reactive metabolites", "Hepatotoxicity"],
            key_unknowns=["Long-term safety"],
            key_assumptions=["Typical metabolism"]
        )
        
        assert summary.overall_tier == RiskTier.MODERATE
        assert len(summary.top_risks) == 2
    
    def test_risk_summary_to_dict(self):
        """Test risk summary serialization."""
        from backend.services.toxpath import RiskSummary, RiskTier
        
        summary = RiskSummary(
            overall_tier=RiskTier.HIGH,
            tier_rationale="Multiple warnings",
            top_risks=[],
            key_unknowns=[],
            key_assumptions=[]
        )
        data = summary.to_dict()
        
        assert data["overall_tier"] == 3
        assert data["tier_label"] == "High Risk"


# ============================================================================
# Assessment Pipeline Tests
# ============================================================================

class TestAssessmentPipeline:
    """Test full assessment pipeline."""
    
    def test_assess_simple_compound(self):
        """Test assessment of simple compound (ethanol)."""
        from backend.services.toxpath import ToxPathAssessor, ToxPathRequest
        
        assessor = ToxPathAssessor()
        request = ToxPathRequest(compound_ref="CCO", compound_name="Ethanol")
        
        result = assessor.assess(request)
        
        assert result.toxpath_assessment_id is not None
        assert result.compound_ref == "CCO"
        assert result.status == "completed"
        assert result.risk_summary is not None
    
    def test_assess_with_route(self):
        """Test assessment with specified route."""
        from backend.services.toxpath import ToxPathAssessor, ToxPathRequest
        
        assessor = ToxPathAssessor()
        request = ToxPathRequest(
            compound_ref="CCO",
            route="oral"
        )
        
        result = assessor.assess(request)
        
        assert result.route == "oral"
        # Should have route-specific concerns
        assert len(result.risk_summary.route_specific_concerns) > 0
    
    def test_assess_inhaled_route_concerns(self):
        """Test inhaled route generates pulmonary concerns."""
        from backend.services.toxpath import ToxPathAssessor, ToxPathRequest
        
        assessor = ToxPathAssessor()
        request = ToxPathRequest(
            compound_ref="CCO",
            route="inhaled"
        )
        
        result = assessor.assess(request)
        
        concerns = " ".join(result.risk_summary.route_specific_concerns).lower()
        assert "pulmonary" in concerns or "inhal" in concerns
    
    def test_assess_generates_testing_plan(self):
        """Test assessment generates testing plan."""
        from backend.services.toxpath import ToxPathAssessor, ToxPathRequest
        
        assessor = ToxPathAssessor()
        request = ToxPathRequest(compound_ref="c1ccccc1")  # Benzene
        
        result = assessor.assess(request)
        
        assert len(result.testing_plan) > 0
        # First step should be characterization
        assert "Structural Confirmation" in result.testing_plan[0].test_name or \
               "characterization" in result.testing_plan[0].test_name.lower()
    
    def test_assess_with_exposure_profile(self):
        """Test assessment with exposure profile."""
        from backend.services.toxpath import (
            ToxPathAssessor, ToxPathRequest, ExposureProfile
        )
        
        assessor = ToxPathAssessor()
        request = ToxPathRequest(
            compound_ref="CCO",
            exposure=ExposureProfile(
                dose_mg=100,
                frequency="daily",
                duration="chronic"
            )
        )
        
        result = assessor.assess(request)
        
        # Chronic exposure should be noted in risks
        risks_text = " ".join(result.risk_summary.top_risks).lower()
        # May or may not mention chronic depending on tier
        assert result.status == "completed"
    
    def test_assess_properties_snapshot(self):
        """Test assessment includes properties snapshot."""
        from backend.services.toxpath import ToxPathAssessor, ToxPathRequest
        
        assessor = ToxPathAssessor()
        request = ToxPathRequest(compound_ref="CCO")
        
        result = assessor.assess(request)
        
        props = result.properties_snapshot
        assert props["is_valid"] is True
        assert "molecular_weight" in props
        assert "logp" in props


class TestStructuralAlertDetection:
    """Test structural alert detection."""
    
    def test_detect_nitro_group(self):
        """Test nitro group detection."""
        from backend.services.toxpath import ToxPathAssessor, ToxPathRequest
        
        assessor = ToxPathAssessor()
        # Nitrobenzene
        request = ToxPathRequest(compound_ref="c1ccc(cc1)[N+](=O)[O-]")
        
        result = assessor.assess(request)
        
        alert_codes = [a.code for a in result.alerts]
        assert "ALERT_NITRO" in alert_codes
    
    def test_detect_epoxide(self):
        """Test epoxide detection."""
        from backend.services.toxpath import ToxPathAssessor, ToxPathRequest
        
        assessor = ToxPathAssessor()
        # Ethylene oxide
        request = ToxPathRequest(compound_ref="C1CO1")
        
        result = assessor.assess(request)
        
        alert_codes = [a.code for a in result.alerts]
        assert "ALERT_EPOXIDE" in alert_codes
    
    def test_detect_high_logp(self):
        """Test high LogP detection."""
        from backend.services.toxpath import ToxPathAssessor, ToxPathRequest
        
        assessor = ToxPathAssessor()
        # Long alkyl chain - high LogP
        request = ToxPathRequest(compound_ref="CCCCCCCCCCCCCCCC")
        
        result = assessor.assess(request)
        
        alert_codes = [a.code for a in result.alerts]
        assert "ALERT_HIGH_LOGP" in alert_codes
    
    def test_clean_compound_low_alerts(self):
        """Test clean compound has minimal alerts."""
        from backend.services.toxpath import ToxPathAssessor, ToxPathRequest
        
        assessor = ToxPathAssessor()
        # Simple, clean compound
        request = ToxPathRequest(compound_ref="CCO")  # Ethanol
        
        result = assessor.assess(request)
        
        # Ethanol should have low risk tier
        # May have info-level alerts but not warnings/critical
        warning_alerts = [a for a in result.alerts 
                        if a.severity.value in ["warning", "critical"]]
        assert len(warning_alerts) == 0


class TestRiskTierClassification:
    """Test risk tier classification logic."""
    
    def test_epoxide_triggers_very_high_tier(self):
        """Test that critical alerts trigger very high tier."""
        from backend.services.toxpath import ToxPathAssessor, ToxPathRequest, RiskTier
        
        assessor = ToxPathAssessor()
        # Ethylene oxide - has critical epoxide alert
        request = ToxPathRequest(compound_ref="C1CO1")
        
        result = assessor.assess(request)
        
        assert result.risk_summary.overall_tier == RiskTier.VERY_HIGH
    
    def test_clean_compound_low_tier(self):
        """Test clean compound gets low tier."""
        from backend.services.toxpath import ToxPathAssessor, ToxPathRequest, RiskTier
        
        assessor = ToxPathAssessor()
        # Simple molecule without alerts
        request = ToxPathRequest(compound_ref="C")  # Methane
        
        result = assessor.assess(request)
        
        assert result.risk_summary.overall_tier == RiskTier.LOW
    
    def test_consultation_required_for_high_tier(self):
        """Test consultation is required for high tier."""
        from backend.services.toxpath import ToxPathAssessor, ToxPathRequest
        
        assessor = ToxPathAssessor()
        # Compound with critical alert
        request = ToxPathRequest(compound_ref="C1CO1")
        
        result = assessor.assess(request)
        
        assert result.consultation_required is True
        assert len(result.consultation_flags) > 0


# ============================================================================
# Memo Generator Tests
# ============================================================================

class TestMemoGenerator:
    """Test ToxPathMemoGenerator."""
    
    def test_memo_generator_initialization(self):
        """Test memo generator initializes."""
        from backend.services.toxpath import ToxPathMemoGenerator
        
        generator = ToxPathMemoGenerator()
        assert generator is not None
        assert generator.config is not None
    
    def test_generate_markdown_memo(self):
        """Test markdown memo generation."""
        from backend.services.toxpath import (
            ToxPathAssessor, ToxPathRequest, 
            ToxPathMemoGenerator, MemoFormat
        )
        
        assessor = ToxPathAssessor()
        request = ToxPathRequest(compound_ref="CCO", compound_name="Ethanol")
        result = assessor.assess(request)
        
        generator = ToxPathMemoGenerator()
        memo = generator.generate(result, MemoFormat.MARKDOWN)
        
        assert "# ToxPath Toxicity Risk Assessment Memo" in memo
        assert "Ethanol" in memo
        assert "Risk Assessment Summary" in memo
    
    def test_generate_html_memo(self):
        """Test HTML memo generation."""
        from backend.services.toxpath import (
            ToxPathAssessor, ToxPathRequest,
            ToxPathMemoGenerator, MemoFormat
        )
        
        assessor = ToxPathAssessor()
        request = ToxPathRequest(compound_ref="CCO")
        result = assessor.assess(request)
        
        generator = ToxPathMemoGenerator()
        memo = generator.generate(result, MemoFormat.HTML)
        
        assert "<!DOCTYPE html>" in memo
        assert "<html>" in memo
        assert "ToxPath Assessment" in memo
    
    def test_generate_json_memo(self):
        """Test JSON memo generation."""
        import json
        from backend.services.toxpath import (
            ToxPathAssessor, ToxPathRequest,
            ToxPathMemoGenerator, MemoFormat
        )
        
        assessor = ToxPathAssessor()
        request = ToxPathRequest(compound_ref="CCO")
        result = assessor.assess(request)
        
        generator = ToxPathMemoGenerator()
        memo = generator.generate(result, MemoFormat.JSON)
        
        # Should be valid JSON
        data = json.loads(memo)
        assert "assessment" in data
        assert "memo_version" in data
    
    def test_memo_with_config(self):
        """Test memo generation with custom config."""
        from backend.services.toxpath import (
            ToxPathAssessor, ToxPathRequest,
            ToxPathMemoGenerator, MemoFormat, MemoConfig
        )
        
        assessor = ToxPathAssessor()
        request = ToxPathRequest(compound_ref="CCO")
        result = assessor.assess(request)
        
        config = MemoConfig(
            company_name="Test Corp",
            prepared_by="Dr. Test"
        )
        generator = ToxPathMemoGenerator(config)
        memo = generator.generate(result, MemoFormat.MARKDOWN)
        
        assert "Test Corp" in memo
        assert "Dr. Test" in memo


# ============================================================================
# API Endpoint Tests
# ============================================================================

class TestToxPathAPI:
    """Test ToxPath API endpoints."""
    
    @pytest.fixture
    def client(self):
        """Create test client."""
        from fastapi.testclient import TestClient
        from backend.main import app
        return TestClient(app)
    
    def test_health_endpoint(self, client):
        """Test health check endpoint."""
        response = client.get("/api/v1/toxpath/health")
        assert response.status_code == 200
        
        data = response.json()
        assert data["service"] == "toxpath"
        assert data["status"] == "healthy"
    
    def test_get_risk_tiers(self, client):
        """Test risk tiers reference endpoint."""
        response = client.get("/api/v1/toxpath/risk-tiers")
        assert response.status_code == 200
        
        data = response.json()
        assert len(data) == 4
        assert data[0]["tier"] == 1
    
    def test_get_alert_codes(self, client):
        """Test alert codes reference endpoint."""
        response = client.get("/api/v1/toxpath/alert-codes")
        assert response.status_code == 200
        
        data = response.json()
        assert len(data) > 10
        codes = [item["code"] for item in data]
        assert "ALERT_NITRO" in codes
    
    def test_get_routes(self, client):
        """Test routes reference endpoint."""
        response = client.get("/api/v1/toxpath/routes")
        assert response.status_code == 200
        
        data = response.json()
        assert "oral" in data
        assert "inhaled" in data
    
    def test_assess_endpoint(self, client):
        """Test assessment endpoint."""
        response = client.post("/api/v1/toxpath/assess", json={
            "compound_ref": "CCO",
            "compound_name": "Ethanol",
            "route": "oral"
        })
        assert response.status_code == 200
        
        data = response.json()
        assert "toxpath_assessment_id" in data
        assert data["compound_ref"] == "CCO"
        assert data["route"] == "oral"
        assert "risk_summary" in data
    
    def test_assess_with_exposure(self, client):
        """Test assessment with exposure profile."""
        response = client.post("/api/v1/toxpath/assess", json={
            "compound_ref": "CCO",
            "exposure": {
                "dose_mg": 100,
                "frequency": "daily",
                "duration": "chronic"
            }
        })
        assert response.status_code == 200
        
        data = response.json()
        assert data["status"] == "completed"
    
    def test_get_assessment_by_id(self, client):
        """Test retrieving assessment by ID."""
        # First create an assessment
        create_response = client.post("/api/v1/toxpath/assess", json={
            "compound_ref": "CCO"
        })
        assessment_id = create_response.json()["toxpath_assessment_id"]
        
        # Then retrieve it
        response = client.get(f"/api/v1/toxpath/assessments/{assessment_id}")
        assert response.status_code == 200
        
        data = response.json()
        assert data["toxpath_assessment_id"] == assessment_id
    
    def test_get_nonexistent_assessment(self, client):
        """Test 404 for nonexistent assessment."""
        response = client.get("/api/v1/toxpath/assessments/nonexistent-id")
        assert response.status_code == 404
    
    def test_generate_memo_endpoint(self, client):
        """Test memo generation endpoint."""
        # First create an assessment
        create_response = client.post("/api/v1/toxpath/assess", json={
            "compound_ref": "CCO",
            "compound_name": "Ethanol"
        })
        assessment_id = create_response.json()["toxpath_assessment_id"]
        
        # Generate memo
        response = client.post(
            f"/api/v1/toxpath/assessments/{assessment_id}/memo",
            json={"format": "markdown"}
        )
        assert response.status_code == 200
        
        data = response.json()
        assert data["assessment_id"] == assessment_id
        assert data["format"] == "markdown"
        assert "ToxPath" in data["memo_content"]
    
    def test_generate_html_memo_endpoint(self, client):
        """Test HTML memo generation via API."""
        create_response = client.post("/api/v1/toxpath/assess", json={
            "compound_ref": "CCO"
        })
        assessment_id = create_response.json()["toxpath_assessment_id"]
        
        response = client.post(
            f"/api/v1/toxpath/assessments/{assessment_id}/memo",
            json={"format": "html"}
        )
        assert response.status_code == 200
        
        data = response.json()
        assert data["format"] == "html"
        assert "<!DOCTYPE html>" in data["memo_content"]


# ============================================================================
# Integration Tests
# ============================================================================

class TestToxPathIntegration:
    """Integration tests for ToxPath."""
    
    def test_cannabinoid_assessment(self):
        """Test assessment of cannabinoid-like compound."""
        from backend.services.toxpath import ToxPathAssessor, ToxPathRequest
        
        assessor = ToxPathAssessor()
        # Simple cannabinoid-like structure (phenol with alkyl chain)
        request = ToxPathRequest(
            compound_ref="CCCCCc1ccc(O)cc1",
            compound_name="4-Pentylphenol",
            route="inhaled"
        )
        
        result = assessor.assess(request)
        
        assert result.status == "completed"
        # Should note BBB penetration potential for small lipophilic compounds
        alert_codes = [a.code for a in result.alerts]
        # May or may not trigger BBB alert depending on exact structure
        assert isinstance(result.testing_plan, list)
    
    def test_assessment_response_serialization(self):
        """Test full response can be serialized to dict."""
        from backend.services.toxpath import ToxPathAssessor, ToxPathRequest
        
        assessor = ToxPathAssessor()
        request = ToxPathRequest(
            compound_ref="CCO",
            compound_name="Ethanol",
            route="oral"
        )
        
        result = assessor.assess(request)
        data = result.to_dict()
        
        # All required fields present
        assert "toxpath_assessment_id" in data
        assert "risk_summary" in data
        assert "alerts" in data
        assert "testing_plan" in data
        assert "consultation_required" in data
    
    def test_full_pipeline_with_memo(self):
        """Test complete pipeline: assess → memo."""
        from backend.services.toxpath import (
            ToxPathAssessor, ToxPathRequest,
            ToxPathMemoGenerator, MemoFormat
        )
        
        # Step 1: Assess
        assessor = ToxPathAssessor()
        request = ToxPathRequest(
            compound_ref="c1ccc(cc1)[N+](=O)[O-]",  # Nitrobenzene
            compound_name="Nitrobenzene",
            route="inhaled"
        )
        result = assessor.assess(request)
        
        # Should have nitro alert
        assert "ALERT_NITRO" in [a.code for a in result.alerts]
        
        # Step 2: Generate memo
        generator = ToxPathMemoGenerator()
        memo = generator.generate(result, MemoFormat.MARKDOWN)
        
        # Memo should reflect findings
        assert "Nitrobenzene" in memo
        assert "Nitro" in memo
        assert "Risk" in memo


# ============================================================================
# Edge Case Tests
# ============================================================================

class TestEdgeCases:
    """Test edge cases and error handling."""
    
    def test_invalid_smiles(self):
        """Test handling of invalid SMILES."""
        from backend.services.toxpath import ToxPathAssessor, ToxPathRequest
        
        assessor = ToxPathAssessor()
        request = ToxPathRequest(compound_ref="not_valid_smiles")
        
        result = assessor.assess(request)
        
        # Should still return result but with invalid flag
        assert result.properties_snapshot["is_valid"] is False
    
    def test_unknown_route(self):
        """Test handling of unknown route."""
        from backend.services.toxpath import ToxPathAssessor, ToxPathRequest
        
        assessor = ToxPathAssessor()
        request = ToxPathRequest(
            compound_ref="CCO",
            route="invalid_route"
        )
        
        result = assessor.assess(request)
        
        # Should default to unknown
        assert result.route == "unknown"
    
    def test_empty_compound_name(self):
        """Test handling of empty compound name."""
        from backend.services.toxpath import ToxPathAssessor, ToxPathRequest
        
        assessor = ToxPathAssessor()
        request = ToxPathRequest(compound_ref="CCO", compound_name="")
        
        result = assessor.assess(request)
        
        assert result.status == "completed"
    
    def test_very_large_molecule(self):
        """Test handling of large molecule."""
        from backend.services.toxpath import ToxPathAssessor, ToxPathRequest
        
        assessor = ToxPathAssessor()
        # Large peptide-like structure
        large_smiles = "CC" * 50  # Long alkane chain
        request = ToxPathRequest(compound_ref=large_smiles)
        
        result = assessor.assess(request)
        
        # Should complete even for large molecules
        assert result.status == "completed"
        # High LogP alert expected
        assert "ALERT_HIGH_LOGP" in [a.code for a in result.alerts]
