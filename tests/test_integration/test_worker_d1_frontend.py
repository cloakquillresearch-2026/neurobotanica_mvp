"""
Integration tests for Cloudflare Worker + D1 + Frontend flow
Tests the complete data flow from frontend through worker to D1 database

These tests validate:
1. Worker endpoints respond correctly
2. D1 database queries execute properly
3. Frontend can fetch and display results
4. Error handling works across the stack
5. Performance meets <10 second target
"""

import pytest
import requests
import time
import json
from typing import Optional
from unittest.mock import patch, MagicMock

# API endpoints for testing
NEUROBOTANICA_API_BASE = "https://neurobotanica-api.contessapetrini.workers.dev"
BUDTENDER_API_BASE = "https://budtender.neuro-botanica.com"


class TestWorkerHealthCheck:
    """Test worker health check endpoint"""

    @pytest.mark.integration
    def test_health_endpoint_returns_status(self):
        """Health check should return status and engine list"""
        response = requests.get(
            f"{NEUROBOTANICA_API_BASE}/api/neurobotanica/health",
            timeout=10
        )

        # Should return 200 or 503 (degraded)
        assert response.status_code in [200, 503]

        data = response.json()
        assert "status" in data
        assert data["status"] in ["healthy", "degraded"]
        assert "engines" in data
        assert isinstance(data["engines"], list)

    @pytest.mark.integration
    def test_health_includes_database_status(self):
        """Health check should report database status"""
        response = requests.get(
            f"{NEUROBOTANICA_API_BASE}/api/neurobotanica/health",
            timeout=10
        )

        if response.status_code == 200:
            data = response.json()
            # May include database status
            if "database" in data:
                assert "healthy" in data["database"]
                assert "latency_ms" in data["database"]


class TestAnalyzeEndpoint:
    """Test the main analyze endpoint"""

    @pytest.mark.integration
    def test_analyze_requires_compound_ids(self):
        """Analyze should reject requests without compound_ids"""
        response = requests.post(
            f"{NEUROBOTANICA_API_BASE}/api/neurobotanica/analyze",
            json={"demographics": {"age": 30}},
            timeout=10
        )

        assert response.status_code == 400
        data = response.json()
        assert "error" in data

    @pytest.mark.integration
    def test_analyze_with_valid_request(self):
        """Analyze should return results for valid request"""
        response = requests.post(
            f"{NEUROBOTANICA_API_BASE}/api/neurobotanica/analyze",
            json={
                "compound_ids": ["cbd", "thc"],
                "demographics": {
                    "age": 35,
                    "gender": "female",
                    "weight": 140,
                    "condition": "anxiety"
                },
                "customer_tier": "computational_only"
            },
            timeout=15
        )

        assert response.status_code == 200
        data = response.json()

        # Check required fields
        assert "interactions" in data
        assert "bias_correction" in data
        assert "synergy" in data
        assert "polysaccharide_effects" in data
        assert "processing_time_ms" in data

        # Check bias correction structure
        bias = data["bias_correction"]
        assert "adjusted_dose_mg" in bias
        assert "factors_applied" in bias
        assert "evidence" in bias

    @pytest.mark.integration
    def test_analyze_performance_under_10_seconds(self):
        """Analysis should complete in under 10 seconds"""
        start_time = time.time()

        response = requests.post(
            f"{NEUROBOTANICA_API_BASE}/api/neurobotanica/analyze",
            json={
                "compound_ids": ["cbd"],
                "demographics": {"age": 30},
                "customer_tier": "computational_only"
            },
            timeout=15
        )

        elapsed = time.time() - start_time

        assert response.status_code == 200
        assert elapsed < 10, f"Analysis took {elapsed:.2f}s, exceeds 10s target"

        # Also check reported processing time
        data = response.json()
        assert data["processing_time_ms"] < 10000


class TestBiasCorrectionEndpoint:
    """Test standalone bias correction endpoint"""

    @pytest.mark.integration
    def test_bias_correction_with_demographics(self):
        """Bias correction should adjust dose based on demographics"""
        response = requests.post(
            f"{NEUROBOTANICA_API_BASE}/api/neurobotanica/bias-correction",
            json={
                "base_dose": 10.0,
                "compound_id": "cbd",
                "demographics": {
                    "age": 70,
                    "gender": "female",
                    "weight": 120
                }
            },
            timeout=10
        )

        assert response.status_code == 200
        data = response.json()

        # Dose should be reduced for elderly female with low weight
        assert data["adjusted_dose_mg"] < 10.0
        assert "factors_applied" in data
        assert len(data["factors_applied"]["demographics_considered"]) > 0

    @pytest.mark.integration
    def test_bias_correction_respects_bounds(self):
        """Bias correction should cap adjustments at safe bounds"""
        response = requests.post(
            f"{NEUROBOTANICA_API_BASE}/api/neurobotanica/bias-correction",
            json={
                "base_dose": 10.0,
                "compound_id": "cbd",
                "demographics": {
                    "age": 85,  # Very elderly
                    "gender": "female",
                    "weight": 90,  # Very light
                    "metabolism_type": "slow"
                }
            },
            timeout=10
        )

        assert response.status_code == 200
        data = response.json()

        # Adjusted dose should be at least 5mg (0.5x cap)
        assert data["adjusted_dose_mg"] >= 5.0

        # Check for warning about capping
        if "warnings" in data and len(data["warnings"]) > 0:
            assert any("cap" in w.lower() for w in data["warnings"])


class TestD1RecommendationsEndpoint:
    """Test D1-backed recommendations endpoint"""

    @pytest.mark.integration
    def test_recommendations_for_valid_condition(self):
        """Should return recommendations for valid condition"""
        response = requests.post(
            f"{BUDTENDER_API_BASE}/api/recommendations",
            json={
                "condition": "ANXIETY",
                "severity": "moderate"
            },
            timeout=10
        )

        assert response.status_code == 200
        data = response.json()

        assert "condition" in data
        assert "evidence_summary" in data
        assert "recommended_cannabinoids" in data
        assert "dosing_guidance" in data

    @pytest.mark.integration
    def test_recommendations_for_invalid_condition(self):
        """Should return 404 with available conditions for invalid condition"""
        response = requests.post(
            f"{BUDTENDER_API_BASE}/api/recommendations",
            json={
                "condition": "INVALID_CONDITION_XYZ",
                "severity": "moderate"
            },
            timeout=10
        )

        assert response.status_code == 404
        data = response.json()

        assert "error" in data
        assert "available_conditions" in data or "suggestion" in data


class TestInflammatorySynergyEndpoint:
    """Test TS-PS-001 inflammatory synergy endpoint"""

    @pytest.mark.integration
    def test_synergy_with_biomarkers(self):
        """Should calculate synergy based on biomarkers"""
        response = requests.post(
            f"{BUDTENDER_API_BASE}/api/dispensary/inflammatory-synergy",
            json={
                "biomarkers": {
                    "tnf_alpha": 15.0,
                    "il6": 8.0,
                    "crp": 5.0,
                    "il1b": 3.0
                },
                "condition_profile": {
                    "conditions": [{"name": "inflammation", "severity": 7, "is_primary": True}],
                    "experience_level": "intermediate"
                },
                "available_kingdoms": ["cannabis", "fungal", "plant", "marine"]
            },
            timeout=10
        )

        assert response.status_code == 200
        data = response.json()

        assert "primary_kingdom" in data
        assert "synergy_score" in data
        assert "recommended_compounds" in data
        assert "expected_reduction" in data

    @pytest.mark.integration
    def test_synergy_without_biomarkers_uses_conditions(self):
        """Should fall back to condition-based logic without biomarkers"""
        response = requests.post(
            f"{BUDTENDER_API_BASE}/api/dispensary/inflammatory-synergy",
            json={
                "biomarkers": {},
                "condition_profile": {
                    "conditions": [{"name": "anxiety", "severity": 6, "is_primary": True}],
                    "experience_level": "beginner"
                },
                "available_kingdoms": ["cannabis", "plant"]
            },
            timeout=10
        )

        assert response.status_code == 200
        data = response.json()

        assert "primary_kingdom" in data
        assert "synergy_score" in data


class TestCORSHandling:
    """Test CORS headers are properly set"""

    @pytest.mark.integration
    def test_cors_preflight_response(self):
        """OPTIONS request should return CORS headers"""
        response = requests.options(
            f"{NEUROBOTANICA_API_BASE}/api/neurobotanica/analyze",
            headers={
                "Origin": "https://neurobotanica.pages.dev",
                "Access-Control-Request-Method": "POST"
            },
            timeout=10
        )

        assert response.status_code == 204
        assert "Access-Control-Allow-Origin" in response.headers
        assert "Access-Control-Allow-Methods" in response.headers

    @pytest.mark.integration
    def test_cors_headers_in_response(self):
        """POST response should include CORS headers"""
        response = requests.post(
            f"{NEUROBOTANICA_API_BASE}/api/neurobotanica/analyze",
            json={"compound_ids": ["cbd"]},
            headers={"Origin": "https://neurobotanica.pages.dev"},
            timeout=10
        )

        assert "Access-Control-Allow-Origin" in response.headers


class TestErrorHandling:
    """Test error handling across the stack"""

    @pytest.mark.integration
    def test_invalid_json_returns_400(self):
        """Invalid JSON should return 400 error"""
        response = requests.post(
            f"{NEUROBOTANICA_API_BASE}/api/neurobotanica/analyze",
            data="not valid json",
            headers={"Content-Type": "application/json"},
            timeout=10
        )

        assert response.status_code == 400

    @pytest.mark.integration
    def test_404_for_unknown_endpoint(self):
        """Unknown endpoint should return 404"""
        response = requests.get(
            f"{NEUROBOTANICA_API_BASE}/api/unknown/endpoint",
            timeout=10
        )

        assert response.status_code == 404
        data = response.json()
        assert "error" in data


class TestFullIntegrationFlow:
    """Test complete frontend -> worker -> D1 flow"""

    @pytest.mark.integration
    @pytest.mark.slow
    def test_complete_budtender_flow(self):
        """
        Simulate complete budtender flow:
        1. Get recommendations for condition
        2. Run NeuroBotanica analysis
        3. Get inflammatory synergy insights
        """
        start_time = time.time()

        # Step 1: Get recommendations
        rec_response = requests.post(
            f"{BUDTENDER_API_BASE}/api/recommendations",
            json={"condition": "CHRONIC_PAIN", "severity": "moderate"},
            timeout=10
        )

        rec_time = time.time() - start_time
        assert rec_response.status_code == 200, f"Recommendations failed: {rec_response.text}"

        # Step 2: Run analysis
        analysis_response = requests.post(
            f"{NEUROBOTANICA_API_BASE}/api/neurobotanica/analyze",
            json={
                "compound_ids": ["cbd", "thc"],
                "demographics": {
                    "age": 45,
                    "gender": "male",
                    "weight": 180,
                    "condition": "chronic_pain"
                },
                "customer_tier": "computational_only"
            },
            timeout=15
        )

        analysis_time = time.time() - start_time
        assert analysis_response.status_code == 200, f"Analysis failed: {analysis_response.text}"

        # Step 3: Get synergy insights
        synergy_response = requests.post(
            f"{BUDTENDER_API_BASE}/api/dispensary/inflammatory-synergy",
            json={
                "biomarkers": {"tnf_alpha": 10.0, "il6": 5.0},
                "condition_profile": {
                    "conditions": [{"name": "chronic_pain", "severity": 7, "is_primary": True}],
                    "experience_level": "regular"
                },
                "available_kingdoms": ["cannabis", "plant"]
            },
            timeout=10
        )

        total_time = time.time() - start_time

        assert synergy_response.status_code == 200, f"Synergy failed: {synergy_response.text}"
        assert total_time < 10, f"Complete flow took {total_time:.2f}s, exceeds 10s target"

        # Validate all responses have required data
        rec_data = rec_response.json()
        analysis_data = analysis_response.json()
        synergy_data = synergy_response.json()

        assert "evidence_summary" in rec_data
        assert "bias_correction" in analysis_data
        assert "synergy_score" in synergy_data

        print(f"\nIntegration flow completed in {total_time:.2f}s")
        print(f"  - Recommendations: {rec_time:.2f}s")
        print(f"  - Analysis: {analysis_time - rec_time:.2f}s")
        print(f"  - Synergy: {total_time - analysis_time:.2f}s")


# Unit tests that can run without network
class TestBiasCorrectionCalculations:
    """Unit tests for bias correction calculations (no network required)"""

    def test_age_factor_young(self):
        """Young adults (<25) should have reduced dose"""
        age = 22
        factor = self._calculate_age_factor(age)
        assert factor == 0.85

    def test_age_factor_adult(self):
        """Adults (25-50) should have standard dose"""
        for age in [25, 30, 40, 50]:
            factor = self._calculate_age_factor(age)
            assert factor == 1.0

    def test_age_factor_senior(self):
        """Seniors (65+) should have reduced dose"""
        factor = self._calculate_age_factor(70)
        assert factor == 0.85

    def test_weight_factor_ranges(self):
        """Weight factors should follow expected ranges"""
        assert self._calculate_weight_factor(100) == 0.85  # Light
        assert self._calculate_weight_factor(140) == 0.9   # Below average
        assert self._calculate_weight_factor(175) == 1.0   # Average
        assert self._calculate_weight_factor(220) == 1.05  # Above average
        assert self._calculate_weight_factor(260) == 1.1   # Heavy

    def test_combined_factors(self):
        """Combined factors should multiply correctly"""
        base_dose = 10.0
        age_factor = 0.85  # 70 years
        weight_factor = 0.9  # 120 lbs
        gender_factor = 0.92  # female

        combined = age_factor * weight_factor * gender_factor
        adjusted = base_dose * combined

        assert abs(combined - 0.7038) < 0.001
        assert abs(adjusted - 7.038) < 0.01

    def test_factor_capping(self):
        """Factors should be capped at safe bounds"""
        assert self._cap_factor(0.3) == 0.5
        assert self._cap_factor(0.8) == 0.8
        assert self._cap_factor(2.5) == 2.0

    @staticmethod
    def _calculate_age_factor(age: int) -> float:
        if age < 25:
            return 0.85
        if age > 75:
            return 0.8
        if age > 65:
            return 0.85
        if age > 50:
            return 0.95
        return 1.0

    @staticmethod
    def _calculate_weight_factor(weight_lbs: int) -> float:
        if weight_lbs < 120:
            return 0.85
        if weight_lbs < 150:
            return 0.9
        if weight_lbs > 280:
            return 1.15
        if weight_lbs > 250:
            return 1.1
        if weight_lbs > 200:
            return 1.05
        return 1.0

    @staticmethod
    def _cap_factor(factor: float) -> float:
        if factor < 0.5:
            return 0.5
        if factor > 2.0:
            return 2.0
        return factor


if __name__ == "__main__":
    pytest.main([__file__, "-v", "--tb=short"])
