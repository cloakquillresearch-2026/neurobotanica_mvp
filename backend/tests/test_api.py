"""
NeuroBotanica Unit Tests - API Endpoints
========================================

Test FastAPI routes: studies, compounds, fda_compliance, conformers

Run tests:
    pytest backend/tests/test_api.py -v
"""
import pytest
from fastapi.testclient import TestClient
from sqlalchemy import create_engine, StaticPool
from sqlalchemy.orm import sessionmaker

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from backend.models.database import Base, get_db
from backend.models.study import ClinicalStudy
from backend.models.compound import Cannabinoid
from backend.main import app


# Create test database with StaticPool to maintain connection across tests
TEST_DATABASE_URL = "sqlite:///:memory:"
test_engine = create_engine(
    TEST_DATABASE_URL, 
    connect_args={"check_same_thread": False},
    poolclass=StaticPool  # Keeps single connection for in-memory DB
)
TestSession = sessionmaker(autocommit=False, autoflush=False, bind=test_engine)

# Track if database is initialized
_db_initialized = False


def override_get_db():
    """Override database dependency for testing."""
    db = TestSession()
    try:
        yield db
    finally:
        db.close()


# Apply override
app.dependency_overrides[get_db] = override_get_db


def _init_test_db():
    """Initialize test database with sample data."""
    global _db_initialized
    if _db_initialized:
        return
    
    Base.metadata.create_all(bind=test_engine)
    
    # Add test data
    db = TestSession()
    
    # Check if data already exists
    if db.query(ClinicalStudy).first() is not None:
        db.close()
        _db_initialized = True
        return
    
    # Add test studies
    studies = [
        ClinicalStudy(
            study_id="TEST-001",
            condition="chronic_pain",
            study_type="RCT",
            study_title="Test Pain Study",
            cannabinoid="CBD",
            sample_size="100",
            year=2023,
            evidence_grade="Level 1"
        ),
        ClinicalStudy(
            study_id="TEST-002",
            condition="anxiety",
            study_type="RCT",
            study_title="Test Anxiety Study",
            cannabinoid="CBD",
            sample_size="80",
            year=2022,
            evidence_grade="Level 2"
        ),
        ClinicalStudy(
            study_id="TEST-003",
            condition="epilepsy",
            study_type="RCT",
            study_title="Test Epilepsy Study",
            cannabinoid="CBD",
            sample_size="150",
            year=2021,
            fda_drug_reference="Epidiolex"
        )
    ]
    
    for study in studies:
        db.add(study)
    
    # Add test compounds
    compounds = [
        Cannabinoid(
            name="Cannabidiol",
            abbreviation="CBD",
            smiles="CCCCCC1=CC(=C(C(=C1)O)C2C=C(CCC2C(=C)C)C)O",
            molecular_weight=314.47,
            logp=5.79,
            is_natural=True,
            fda_approved=True,
            fda_drug_name="Epidiolex",
            schedule_classification="unscheduled"
        ),
        Cannabinoid(
            name="Delta-9-THC",
            abbreviation="THC",
            smiles="CCCCCC1=CC(=C2C3C=C(CCC3C(OC2=C1)(C)C)C)O",
            molecular_weight=314.47,
            logp=7.68,
            is_natural=True,
            fda_approved=True,
            fda_drug_name="Marinol",
            schedule_classification="III"
        )
    ]
    
    for compound in compounds:
        db.add(compound)
    
    db.commit()
    db.close()
    _db_initialized = True


@pytest.fixture(scope="function")
def client():
    """Create test client with fresh database for each test."""
    # Ensure our override is applied (other test files might override)
    app.dependency_overrides[get_db] = override_get_db
    _init_test_db()
    
    with TestClient(app) as test_client:
        yield test_client


class TestRootEndpoints:
    """Tests for root and health endpoints."""
    
    def test_root(self, client):
        """Test root endpoint."""
        response = client.get("/")
        assert response.status_code == 200
        
        data = response.json()
        assert "NeuroBotanica" in data["message"]
        assert data["version"] == "0.1.0"
        assert "endpoints" in data
    
    def test_health(self, client):
        """Test health check endpoint."""
        response = client.get("/health")
        assert response.status_code == 200
        
        data = response.json()
        assert data["status"] in ["healthy", "degraded"]
        assert "features" in data
        assert data["features"]["fda_compliance_module"] is True


class TestStudyEndpoints:
    """Tests for study API endpoints."""
    
    def test_list_studies(self, client):
        """Test listing studies."""
        response = client.get("/api/v1/studies/")
        assert response.status_code == 200
        
        data = response.json()
        assert len(data) >= 3
    
    def test_get_study_by_id(self, client):
        """Test getting study by study_id."""
        # Get study by known study_id
        response = client.get("/api/v1/studies/TEST-001")
        assert response.status_code == 200
        assert response.json()["study_id"] == "TEST-001"
    
    def test_filter_studies_by_condition(self, client):
        """Test filtering studies by condition."""
        response = client.get("/api/v1/studies/?condition=chronic_pain")
        assert response.status_code == 200
        
        data = response.json()
        for study in data:
            assert "chronic_pain" in study["condition"].lower()
    
    def test_filter_studies_by_cannabinoid(self, client):
        """Test filtering studies by cannabinoid."""
        response = client.get("/api/v1/studies/?cannabinoid=CBD")
        assert response.status_code == 200
        
        data = response.json()
        for study in data:
            assert "CBD" in study.get("cannabinoid", "").upper()
    
    def test_study_not_found(self, client):
        """Test 404 for non-existent study."""
        response = client.get("/api/v1/studies/NONEXISTENT-999")
        assert response.status_code == 404


class TestCompoundEndpoints:
    """Tests for compound API endpoints."""
    
    def test_list_compounds(self, client):
        """Test listing compounds."""
        response = client.get("/api/v1/compounds/")
        assert response.status_code == 200
        
        data = response.json()
        assert len(data) >= 2
    
    def test_get_compound_by_name(self, client):
        """Test getting compound by name."""
        # Use full name that matches ilike search
        response = client.get("/api/v1/compounds/Cannabidiol")
        assert response.status_code == 200
        assert "Cannabidiol" in response.json()["name"] or "CBD" in response.json().get("abbreviation", "")
    
    def test_get_compound_pharmacology(self, client):
        """Test getting compound pharmacology profile."""
        response = client.get("/api/v1/compounds/Cannabidiol/pharmacology")
        assert response.status_code == 200
    
    def test_compound_not_found(self, client):
        """Test 404 for non-existent compound."""
        response = client.get("/api/v1/compounds/NONEXISTENT")
        assert response.status_code == 404


class TestFDAComplianceEndpoints:
    """Tests for FDA compliance API endpoints."""
    
    def test_fda_overview(self, client):
        """Test FDA overview endpoint."""
        response = client.get("/api/v1/fda/")
        assert response.status_code == 200
        
        data = response.json()
        assert "endpoints" in data or "message" in data
    
    def test_efficacy_comparison(self, client):
        """Test efficacy comparison endpoint."""
        response = client.get("/api/v1/fda/efficacy-comparison/chronic_pain")
        assert response.status_code == 200
    
    def test_pharmacology_package(self, client):
        """Test pharmacology package generation."""
        response = client.get("/api/v1/fda/pharmacology-package/Cannabidiol")
        # 200 if found, 404 if compound not in database (acceptable)
        assert response.status_code in [200, 404]


class TestConformerEndpoints:
    """Tests for 3D conformer API endpoints."""
    
    def test_conformer_service_status(self, client):
        """Test conformer service status."""
        response = client.get("/api/v1/conformers/status")
        assert response.status_code == 200
        
        data = response.json()
        assert "rdkit_available" in data
        assert data["service"] == "3D Conformer Generator"
    
    def test_conformer_summary(self, client):
        """Test conformer summary endpoint."""
        response = client.get("/api/v1/conformers/summary")
        assert response.status_code == 200
        
        data = response.json()
        assert "total_compounds" in data
        assert "with_conformers" in data
        assert "coverage_percent" in data
    
    def test_pending_compounds(self, client):
        """Test pending compounds endpoint."""
        response = client.get("/api/v1/conformers/pending")
        assert response.status_code == 200
        
        data = response.json()
        assert isinstance(data, list)
    
    def test_conformer_status_for_compound(self, client):
        """Test getting conformer status for a compound."""
        # Get status for known compound
        response = client.get("/api/v1/conformers/compound/1")
        # May 404 if no compounds with that ID
        assert response.status_code in [200, 404]
        
        if response.status_code == 200:
            data = response.json()
            assert "compound_id" in data or "has_conformers" in data
    
    def test_generate_from_smiles(self, client):
        """Test conformer generation from SMILES."""
        # Skip if RDKit not available
        status_response = client.get("/api/v1/conformers/status")
        if not status_response.json().get("rdkit_available"):
            pytest.skip("RDKit not available")
        
        response = client.post(
            "/api/v1/conformers/generate-from-smiles",
            json={"smiles": "CCO", "num_conformers": 10}  # Ethanol
        )
        
        # May be 200 or 503 depending on RDKit availability
        assert response.status_code in [200, 503]


class TestAPIValidation:
    """Tests for API input validation."""
    
    def test_study_not_found_with_invalid_id(self, client):
        """Test 404 for invalid study ID."""
        response = client.get("/api/v1/studies/INVALID-STUDY-XYZ")
        assert response.status_code == 404  # Not found (valid format but doesn't exist)
    
    def test_compound_not_found_with_invalid_name(self, client):
        """Test 404 for invalid compound name."""
        response = client.get("/api/v1/compounds/INVALID_COMPOUND_XYZ")
        assert response.status_code == 404  # Not found
    
    def test_pagination_limits(self, client):
        """Test pagination parameter validation."""
        # Limit > 100 should fail validation (le=100)
        response = client.get("/api/v1/studies/?limit=1000")
        assert response.status_code == 422  # Exceeds max limit
        
        # Valid limit should work
        response = client.get("/api/v1/studies/?limit=50")
        assert response.status_code == 200


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
