"""
Integration Tests for Dispensary API (Use Case 4)

Tests the Nevada pilot dispensary recommendation system.

Run with: pytest tests/integration/test_dispensary_api.py -v
"""

import pytest
from fastapi.testclient import TestClient
import sys
from pathlib import Path

project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))

from backend.main import app

client = TestClient(app)


class TestDispensaryRecommendationAPI:
    """Tests for dispensary personalized recommendation system."""
    
    def test_recommend_endpoint_exists(self):
        """Test that recommend endpoint is registered."""
        response = client.post(
            "/api/dispensary/recommend",
            json={
                "customer_profile": {
                    "age": 42,
                    "weight_kg": 70,
                    "conditions": [
                        {"name": "chronic_pain", "severity": 7, "is_primary": True}
                    ],
                    "experience_level": "regular"
                },
                "available_inventory": [
                    {
                        "product_id": "SKU_001",
                        "product_name": "ACDC Vape",
                        "product_type": "vape",
                        "cbd_percent": 20,
                        "thc_percent": 1
                    }
                ],
                "max_recommendations": 3
            }
        )
        assert response.status_code == 200
        data = response.json()
        assert "recommendation_id" in data
        assert "recommendations" in data
    
    def test_recommend_returns_match_score(self):
        """Test that recommendations include match scores."""
        response = client.post(
            "/api/dispensary/recommend",
            json={
                "customer_profile": {
                    "age": 35,
                    "weight_kg": 75,
                    "conditions": [
                        {"name": "anxiety", "severity": 6, "is_primary": True}
                    ],
                    "experience_level": "occasional"
                },
                "available_inventory": [
                    {
                        "product_id": "SKU_002",
                        "product_name": "High CBD Tincture",
                        "product_type": "tincture",
                        "cbd_percent": 25,
                        "thc_percent": 0
                    }
                ]
            }
        )
        assert response.status_code == 200
        data = response.json()
        assert len(data["recommendations"]) > 0
        rec = data["recommendations"][0]
        assert "match_score" in rec
        assert 0 <= rec["match_score"] <= 1
    
    def test_profile_creation(self):
        """Test customer profile creation."""
        response = client.post(
            "/api/dispensary/profile",
            json={
                "age": 45,
                "weight_kg": 80,
                "conditions": [
                    {"name": "insomnia", "severity": 8, "is_primary": True}
                ],
                "experience_level": "regular"
            }
        )
        assert response.status_code == 200
        data = response.json()
        assert "profile_id" in data
        assert data["profile_id"].startswith("prof_")
        assert "completeness_score" in data
    
    def test_statistics_endpoint(self):
        """Test dispensary statistics endpoint."""
        response = client.get("/api/dispensary/statistics")
        assert response.status_code == 200
        data = response.json()
        assert "clinical_studies_used" in data
        assert data["clinical_studies_used"] == 398
        assert data["conditions_covered"] == 22


class TestAdjuvantOptimizationAPI:
    """Tests for adjuvant optimization endpoint."""
    
    def test_adjuvant_optimize_endpoint_exists(self):
        """Test that adjuvant optimization endpoint is registered."""
        response = client.post(
            "/api/dispensary/adjuvants/optimize",
            json={
                "primary_compound": "CBD",
                "therapeutic_target": "insomnia"
            }
        )
        assert response.status_code == 200
        data = response.json()
        assert "recommended_adjuvant" in data
        assert "dosage_mg" in data
        assert "timing_offset_minutes" in data
    
    def test_adjuvant_for_chronic_pain(self):
        """Test adjuvant recommendation for chronic pain."""
        response = client.post(
            "/api/dispensary/adjuvants/optimize",
            json={
                "primary_compound": "THC",
                "therapeutic_target": "chronic_pain"
            }
        )
        assert response.status_code == 200
        data = response.json()
        assert data["expected_enhancement_percent"] > 0
        assert data["evidence_tier"] in [1, 2, 3, 4, 5]
    
    def test_adjuvant_with_patient_profile(self):
        """Test adjuvant with patient profile for personalization."""
        response = client.post(
            "/api/dispensary/adjuvants/optimize",
            json={
                "primary_compound": "CBD",
                "therapeutic_target": "anxiety",
                "patient_profile": {
                    "age": 55,
                    "weight_kg": 65,
                    "conditions": [
                        {"name": "anxiety", "severity": 6, "is_primary": True}
                    ],
                    "supplements": ["magnesium"],
                    "experience_level": "regular"
                }
            }
        )
        assert response.status_code == 200
        data = response.json()
        assert "mechanism" in data
        assert len(data["citations"]) > 0


class TestProductAvoidance:
    """Tests for product avoidance logic."""
    
    def test_avoids_anxiety_triggering_products(self):
        """Test that high-anxiety-risk products are flagged to avoid."""
        response = client.post(
            "/api/dispensary/recommend",
            json={
                "customer_profile": {
                    "age": 30,
                    "weight_kg": 60,
                    "conditions": [
                        {"name": "anxiety", "severity": 8, "is_primary": True}
                    ],
                    "negative_experiences": [
                        {"strain": "Sour Diesel", "issue": "increased_anxiety"}
                    ],
                    "experience_level": "occasional"
                },
                "available_inventory": [
                    {
                        "product_id": "SKU_GOOD",
                        "product_name": "High CBD Calm",
                        "product_type": "tincture",
                        "cbd_percent": 20,
                        "thc_percent": 1
                    },
                    {
                        "product_id": "SKU_BAD",
                        "product_name": "Sour Diesel Vape",
                        "product_type": "vape",
                        "strain_name": "Sour Diesel",
                        "cbd_percent": 0,
                        "thc_percent": 26
                    }
                ]
            }
        )
        assert response.status_code == 200
        data = response.json()
        
        # Check that Sour Diesel is in avoid list
        avoid_ids = [p["product_id"] for p in data["products_to_avoid"]]
        assert "SKU_BAD" in avoid_ids


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
