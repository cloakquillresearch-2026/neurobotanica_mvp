"""
Week 11 Tests - RegPath Regulatory Pathway Optimization Engine

Comprehensive tests for:
- RegPath strategist functionality
- Pathway decision matrix
- Timeline generation
- Checklist generation
- Memo generation
- API endpoints
"""

import pytest
from datetime import datetime


# ============================================================================
# RegPath Strategist Tests
# ============================================================================

class TestRegPathStrategist:
    """Test RegPathStrategist core functionality."""
    
    def test_strategist_initialization(self):
        """Test strategist initializes correctly."""
        from backend.services.regpath import RegPathStrategist
        
        strategist = RegPathStrategist()
        assert strategist is not None
    
    def test_strategist_has_pathways(self):
        """Test strategist provides pathway info."""
        from backend.services.regpath import RegPathStrategist
        
        strategist = RegPathStrategist()
        pathways = strategist.get_pathways()
        
        assert len(pathways) >= 5
        pathway_names = [p["pathway"] for p in pathways]
        assert "IND" in pathway_names
        assert "NDA" in pathway_names
        assert "505(b)(2)" in pathway_names
    
    def test_strategist_has_product_types(self):
        """Test strategist provides product types."""
        from backend.services.regpath import RegPathStrategist
        
        strategist = RegPathStrategist()
        types = strategist.get_product_types()
        
        assert len(types) >= 4
        type_values = [t["type"] for t in types]
        assert "isolated_cannabinoid" in type_values
        assert "novel_dimer" in type_values
        assert "formulation" in type_values
    
    def test_strategist_has_checklist_categories(self):
        """Test strategist provides checklist categories."""
        from backend.services.regpath import RegPathStrategist
        
        strategist = RegPathStrategist()
        categories = strategist.get_checklist_categories()
        
        assert "CMC" in categories
        assert "Nonclinical" in categories
        assert "Clinical" in categories
        assert "Documentation" in categories


class TestProductTypeEnum:
    """Test ProductType enumeration."""
    
    def test_product_type_values(self):
        """Test product type values."""
        from backend.services.regpath import ProductType
        
        assert ProductType.ISOLATED_CANNABINOID.value == "isolated_cannabinoid"
        assert ProductType.NOVEL_DIMER.value == "novel_dimer"
        assert ProductType.FORMULATION.value == "formulation"
    
    def test_product_type_display_names(self):
        """Test product type display names."""
        from backend.services.regpath import ProductType
        
        assert "Isolated Cannabinoid" in ProductType.ISOLATED_CANNABINOID.display_name
        assert "Novel" in ProductType.NOVEL_DIMER.display_name


class TestRegulatoryPathwayEnum:
    """Test RegulatoryPathway enumeration."""
    
    def test_pathway_values(self):
        """Test pathway values."""
        from backend.services.regpath import RegulatoryPathway
        
        assert RegulatoryPathway.IND.value == "IND"
        assert RegulatoryPathway.NDA.value == "NDA"
        assert RegulatoryPathway.NDA_505B2.value == "505(b)(2)"
        assert RegulatoryPathway.ANDA.value == "ANDA"
    
    def test_pathway_full_names(self):
        """Test pathway full names."""
        from backend.services.regpath import RegulatoryPathway
        
        assert "Investigational New Drug" in RegulatoryPathway.IND.full_name
        assert "New Drug Application" in RegulatoryPathway.NDA.full_name
    
    def test_pathway_timelines(self):
        """Test pathway typical timelines."""
        from backend.services.regpath import RegulatoryPathway
        
        ind_timeline = RegulatoryPathway.IND.typical_timeline_months
        assert ind_timeline[0] >= 12
        assert ind_timeline[1] <= 24
        
        nda_timeline = RegulatoryPathway.NDA.typical_timeline_months
        assert nda_timeline[0] >= 36


class TestNoveltyLevelEnum:
    """Test NoveltyLevel enumeration."""
    
    def test_novelty_values(self):
        """Test novelty level values."""
        from backend.services.regpath import NoveltyLevel
        
        assert NoveltyLevel.LOW.value == "low"
        assert NoveltyLevel.MODERATE.value == "moderate"
        assert NoveltyLevel.HIGH.value == "high"
    
    def test_novelty_descriptions(self):
        """Test novelty level descriptions."""
        from backend.services.regpath import NoveltyLevel
        
        assert "existing safety data" in NoveltyLevel.LOW.description.lower()
        assert "novel" in NoveltyLevel.HIGH.description.lower()


# ============================================================================
# Data Model Tests
# ============================================================================

class TestProductProfile:
    """Test ProductProfile data model."""
    
    def test_profile_creation(self):
        """Test product profile creation."""
        from backend.services.regpath import ProductProfile
        
        profile = ProductProfile(
            product_name="CBD-X",
            product_type="isolated_cannabinoid",
            intended_use="Treatment of epilepsy",
            route="oral"
        )
        
        assert profile.product_name == "CBD-X"
        assert profile.route == "oral"
        assert profile.target_markets == ["US"]
    
    def test_profile_to_dict(self):
        """Test profile serialization."""
        from backend.services.regpath import ProductProfile
        
        profile = ProductProfile(
            product_name="Test Product",
            product_type="formulation",
            intended_use="Pain relief",
            route="topical",
            therapeutic_area="Pain Management"
        )
        data = profile.to_dict()
        
        assert data["product_name"] == "Test Product"
        assert data["therapeutic_area"] == "Pain Management"


class TestEvidenceInputs:
    """Test EvidenceInputs data model."""
    
    def test_evidence_creation(self):
        """Test evidence inputs creation."""
        from backend.services.regpath import EvidenceInputs
        
        evidence = EvidenceInputs(
            chempath_job_id="job-123",
            toxpath_assessment_id="tox-456",
            predicate_drug="Epidiolex"
        )
        
        assert evidence.chempath_job_id == "job-123"
        assert evidence.predicate_drug == "Epidiolex"
    
    def test_evidence_to_dict(self):
        """Test evidence serialization."""
        from backend.services.regpath import EvidenceInputs
        
        evidence = EvidenceInputs(
            existing_literature_refs=["ref1", "ref2"]
        )
        data = evidence.to_dict()
        
        assert len(data["existing_literature_refs"]) == 2


class TestStrategyConstraints:
    """Test StrategyConstraints data model."""
    
    def test_constraints_creation(self):
        """Test constraints creation."""
        from backend.services.regpath import StrategyConstraints
        
        constraints = StrategyConstraints(
            budget_range="10M-50M",
            risk_tolerance="low",
            excluded_pathways=["ANDA"]
        )
        
        assert constraints.budget_range == "10M-50M"
        assert constraints.risk_tolerance == "low"
        assert "ANDA" in constraints.excluded_pathways


class TestPathwayRecommendation:
    """Test PathwayRecommendation data model."""
    
    def test_recommendation_creation(self):
        """Test pathway recommendation creation."""
        from backend.services.regpath import PathwayRecommendation, RegulatoryPathway
        
        rec = PathwayRecommendation(
            pathway=RegulatoryPathway.NDA_505B2,
            confidence="high",
            rationale=["Strong predicate drug", "Established safety"],
            key_requirements=["BA/BE study"],
            estimated_cost_range="$15M - $75M",
            estimated_timeline_months=(24, 42)
        )
        
        assert rec.pathway == RegulatoryPathway.NDA_505B2
        assert rec.confidence == "high"
    
    def test_recommendation_to_dict(self):
        """Test recommendation serialization."""
        from backend.services.regpath import PathwayRecommendation, RegulatoryPathway
        
        rec = PathwayRecommendation(
            pathway=RegulatoryPathway.IND,
            confidence="moderate",
            rationale=["Novel compound"],
            key_requirements=["GLP tox"],
            estimated_cost_range="$2M - $10M",
            estimated_timeline_months=(12, 24)
        )
        data = rec.to_dict()
        
        assert data["pathway"] == "IND"
        assert data["estimated_timeline_months"]["min"] == 12


class TestChecklistItem:
    """Test ChecklistItem data model."""
    
    def test_checklist_item_creation(self):
        """Test checklist item creation."""
        from backend.services.regpath import ChecklistItem, ChecklistCategory
        
        item = ChecklistItem(
            category=ChecklistCategory.CMC,
            item="Drug Substance Spec",
            description="Establish specifications",
            priority="required"
        )
        
        assert item.category == ChecklistCategory.CMC
        assert item.priority == "required"


class TestTimelineMilestone:
    """Test TimelineMilestone data model."""
    
    def test_milestone_creation(self):
        """Test milestone creation."""
        from backend.services.regpath import TimelineMilestone, MilestonePhase
        
        milestone = TimelineMilestone(
            phase=MilestonePhase.IND_ENABLING,
            milestone="GLP Toxicology",
            description="Complete GLP tox studies",
            duration_months=9,
            start_month=3,
            deliverables=["Tox report"]
        )
        
        assert milestone.phase == MilestonePhase.IND_ENABLING
        assert milestone.duration_months == 9
    
    def test_milestone_to_dict(self):
        """Test milestone serialization."""
        from backend.services.regpath import TimelineMilestone, MilestonePhase
        
        milestone = TimelineMilestone(
            phase=MilestonePhase.PHASE_1,
            milestone="Phase 1 Trial",
            description="First in human",
            duration_months=12,
            start_month=12
        )
        data = milestone.to_dict()
        
        assert data["phase"] == "Phase 1"
        assert data["end_month"] == 24


# ============================================================================
# Strategy Generation Tests
# ============================================================================

class TestStrategyGeneration:
    """Test strategy generation pipeline."""
    
    def test_generate_strategy_isolated_cannabinoid(self):
        """Test strategy for isolated cannabinoid."""
        from backend.services.regpath import (
            RegPathStrategist, RegPathRequest, ProductProfile
        )
        
        strategist = RegPathStrategist()
        request = RegPathRequest(
            product_profile=ProductProfile(
                product_name="CBD Solution",
                product_type="isolated_cannabinoid",
                intended_use="Epilepsy treatment",
                route="oral"
            )
        )
        
        result = strategist.generate_strategy(request)
        
        assert result.regpath_strategy_id is not None
        assert result.status == "completed"
        assert result.primary_pathway is not None
    
    def test_generate_strategy_novel_dimer(self):
        """Test strategy for novel dimer (always IND)."""
        from backend.services.regpath import (
            RegPathStrategist, RegPathRequest, ProductProfile,
            RegulatoryPathway
        )
        
        strategist = RegPathStrategist()
        request = RegPathRequest(
            product_profile=ProductProfile(
                product_name="Dimer-X",
                product_type="novel_dimer",
                intended_use="Pain management",
                route="sublingual"
            )
        )
        
        result = strategist.generate_strategy(request)
        
        # Novel dimers should always get IND as primary
        assert result.primary_pathway.pathway == RegulatoryPathway.IND
        assert result.novelty_assessment == "high"
    
    def test_generate_strategy_with_evidence(self):
        """Test strategy with evidence inputs."""
        from backend.services.regpath import (
            RegPathStrategist, RegPathRequest, ProductProfile, EvidenceInputs
        )
        
        strategist = RegPathStrategist()
        request = RegPathRequest(
            product_profile=ProductProfile(
                product_name="CBD-Plus",
                product_type="formulation",
                intended_use="Anxiety treatment",
                route="oral"
            ),
            evidence_inputs=EvidenceInputs(
                chempath_job_id="chem-123",
                toxpath_assessment_id="tox-456",
                predicate_drug="Epidiolex"
            )
        )
        
        result = strategist.generate_strategy(request)
        
        # With predicate drug, should have stronger evidence
        assert result.evidence_strength in ["moderate", "strong"]
    
    def test_generate_strategy_with_constraints(self):
        """Test strategy with constraints."""
        from backend.services.regpath import (
            RegPathStrategist, RegPathRequest, ProductProfile, StrategyConstraints
        )
        
        strategist = RegPathStrategist()
        request = RegPathRequest(
            product_profile=ProductProfile(
                product_name="Test Product",
                product_type="isolated_cannabinoid",
                intended_use="Test use",
                route="oral"
            ),
            constraints=StrategyConstraints(
                budget_range="5M-20M",
                risk_tolerance="low"
            )
        )
        
        result = strategist.generate_strategy(request)
        
        assert result.status == "completed"
        assert len(result.key_assumptions) > 0


class TestPathwaySelection:
    """Test pathway selection logic."""
    
    def test_505b2_with_predicate(self):
        """Test 505(b)(2) selection with predicate drug."""
        from backend.services.regpath import (
            RegPathStrategist, RegPathRequest, ProductProfile, 
            EvidenceInputs, RegulatoryPathway
        )
        
        strategist = RegPathStrategist()
        request = RegPathRequest(
            product_profile=ProductProfile(
                product_name="CBD Generic",
                product_type="formulation",
                intended_use="Epilepsy",
                route="oral"
            ),
            evidence_inputs=EvidenceInputs(
                predicate_drug="Epidiolex",
                chempath_job_id="chem-1",
                toxpath_assessment_id="tox-1"
            )
        )
        
        result = strategist.generate_strategy(request)
        
        # Strong evidence + low novelty formulation should get 505(b)(2) or ANDA
        assert result.primary_pathway.pathway in [
            RegulatoryPathway.NDA_505B2, 
            RegulatoryPathway.ANDA
        ]
    
    def test_ind_for_weak_evidence(self):
        """Test IND selection for weak evidence."""
        from backend.services.regpath import (
            RegPathStrategist, RegPathRequest, ProductProfile,
            RegulatoryPathway
        )
        
        strategist = RegPathStrategist()
        request = RegPathRequest(
            product_profile=ProductProfile(
                product_name="New Compound",
                product_type="synthetic",
                intended_use="Pain treatment",
                route="oral"
            )
            # No evidence inputs = weak evidence
        )
        
        result = strategist.generate_strategy(request)
        
        # Weak evidence + synthetic = likely IND
        assert result.evidence_strength == "weak"


class TestChecklistGeneration:
    """Test checklist generation."""
    
    def test_checklist_has_cmc_items(self):
        """Test checklist includes CMC items."""
        from backend.services.regpath import (
            RegPathStrategist, RegPathRequest, ProductProfile,
            ChecklistCategory
        )
        
        strategist = RegPathStrategist()
        request = RegPathRequest(
            product_profile=ProductProfile(
                product_name="Test",
                product_type="isolated_cannabinoid",
                intended_use="Test",
                route="oral"
            )
        )
        
        result = strategist.generate_strategy(request)
        
        cmc_items = [i for i in result.readiness_checklist 
                    if i.category == ChecklistCategory.CMC]
        assert len(cmc_items) >= 3
    
    def test_checklist_priorities(self):
        """Test checklist has required items."""
        from backend.services.regpath import (
            RegPathStrategist, RegPathRequest, ProductProfile
        )
        
        strategist = RegPathStrategist()
        request = RegPathRequest(
            product_profile=ProductProfile(
                product_name="Test",
                product_type="formulation",
                intended_use="Test",
                route="oral"
            )
        )
        
        result = strategist.generate_strategy(request)
        
        required_items = [i for i in result.readiness_checklist 
                        if i.priority == "required"]
        assert len(required_items) > 0


class TestTimelineGeneration:
    """Test timeline generation."""
    
    def test_timeline_has_milestones(self):
        """Test timeline includes milestones."""
        from backend.services.regpath import (
            RegPathStrategist, RegPathRequest, ProductProfile
        )
        
        strategist = RegPathStrategist()
        request = RegPathRequest(
            product_profile=ProductProfile(
                product_name="Test",
                product_type="isolated_cannabinoid",
                intended_use="Test",
                route="oral"
            )
        )
        
        result = strategist.generate_strategy(request)
        
        assert len(result.timeline) >= 3
    
    def test_timeline_sequence(self):
        """Test timeline milestones are sequential."""
        from backend.services.regpath import (
            RegPathStrategist, RegPathRequest, ProductProfile
        )
        
        strategist = RegPathStrategist()
        request = RegPathRequest(
            product_profile=ProductProfile(
                product_name="Test",
                product_type="novel_dimer",
                intended_use="Test",
                route="oral"
            )
        )
        
        result = strategist.generate_strategy(request)
        
        # Check milestones don't have negative start months
        for milestone in result.timeline:
            assert milestone.start_month >= 0


# ============================================================================
# Memo Generator Tests
# ============================================================================

class TestMemoGenerator:
    """Test RegPathMemoGenerator."""
    
    def test_memo_generator_initialization(self):
        """Test memo generator initializes."""
        from backend.services.regpath import RegPathMemoGenerator
        
        generator = RegPathMemoGenerator()
        assert generator is not None
        assert generator.config is not None
    
    def test_generate_markdown_memo(self):
        """Test markdown memo generation."""
        from backend.services.regpath import (
            RegPathStrategist, RegPathRequest, ProductProfile,
            RegPathMemoGenerator, MemoFormat
        )
        
        strategist = RegPathStrategist()
        request = RegPathRequest(
            product_profile=ProductProfile(
                product_name="CBD Solution",
                product_type="isolated_cannabinoid",
                intended_use="Epilepsy treatment",
                route="oral"
            )
        )
        result = strategist.generate_strategy(request)
        
        generator = RegPathMemoGenerator()
        memo = generator.generate(result, MemoFormat.MARKDOWN)
        
        assert "# Regulatory Strategy Memo" in memo
        assert "CBD Solution" in memo
        assert "Pathway Recommendation" in memo
    
    def test_generate_html_memo(self):
        """Test HTML memo generation."""
        from backend.services.regpath import (
            RegPathStrategist, RegPathRequest, ProductProfile,
            RegPathMemoGenerator, MemoFormat
        )
        
        strategist = RegPathStrategist()
        request = RegPathRequest(
            product_profile=ProductProfile(
                product_name="Test Product",
                product_type="formulation",
                intended_use="Test",
                route="oral"
            )
        )
        result = strategist.generate_strategy(request)
        
        generator = RegPathMemoGenerator()
        memo = generator.generate(result, MemoFormat.HTML)
        
        assert "<!DOCTYPE html>" in memo
        assert "<html>" in memo
    
    def test_generate_json_memo(self):
        """Test JSON memo generation."""
        import json
        from backend.services.regpath import (
            RegPathStrategist, RegPathRequest, ProductProfile,
            RegPathMemoGenerator, MemoFormat
        )
        
        strategist = RegPathStrategist()
        request = RegPathRequest(
            product_profile=ProductProfile(
                product_name="Test",
                product_type="isolated_cannabinoid",
                intended_use="Test",
                route="oral"
            )
        )
        result = strategist.generate_strategy(request)
        
        generator = RegPathMemoGenerator()
        memo = generator.generate(result, MemoFormat.JSON)
        
        data = json.loads(memo)
        assert "strategy" in data
        assert "memo_version" in data
    
    def test_memo_with_config(self):
        """Test memo with custom config."""
        from backend.services.regpath import (
            RegPathStrategist, RegPathRequest, ProductProfile,
            RegPathMemoGenerator, MemoFormat, MemoConfig
        )
        
        strategist = RegPathStrategist()
        request = RegPathRequest(
            product_profile=ProductProfile(
                product_name="Test",
                product_type="formulation",
                intended_use="Test",
                route="oral"
            )
        )
        result = strategist.generate_strategy(request)
        
        config = MemoConfig(
            company_name="Test Pharma",
            prepared_by="Dr. Test"
        )
        generator = RegPathMemoGenerator(config)
        memo = generator.generate(result, MemoFormat.MARKDOWN)
        
        assert "Test Pharma" in memo
        assert "Dr. Test" in memo


# ============================================================================
# API Endpoint Tests
# ============================================================================

class TestRegPathAPI:
    """Test RegPath API endpoints."""
    
    @pytest.fixture
    def client(self):
        """Create test client."""
        from fastapi.testclient import TestClient
        from backend.main import app
        return TestClient(app)
    
    def test_health_endpoint(self, client):
        """Test health check endpoint."""
        response = client.get("/api/v1/regpath/health")
        assert response.status_code == 200
        
        data = response.json()
        assert data["service"] == "regpath"
        assert data["status"] == "healthy"
    
    def test_get_pathways(self, client):
        """Test pathways reference endpoint."""
        response = client.get("/api/v1/regpath/pathways")
        assert response.status_code == 200
        
        data = response.json()
        assert len(data) >= 5
        pathways = [p["pathway"] for p in data]
        assert "IND" in pathways
    
    def test_get_product_types(self, client):
        """Test product types reference endpoint."""
        response = client.get("/api/v1/regpath/product-types")
        assert response.status_code == 200
        
        data = response.json()
        assert len(data) >= 4
    
    def test_get_checklist_categories(self, client):
        """Test checklist categories endpoint."""
        response = client.get("/api/v1/regpath/checklist-categories")
        assert response.status_code == 200
        
        data = response.json()
        assert "CMC" in data
        assert "Clinical" in data
    
    def test_generate_strategy_endpoint(self, client):
        """Test strategy generation endpoint."""
        response = client.post("/api/v1/regpath/strategy", json={
            "product_profile": {
                "product_name": "CBD Test Solution",
                "product_type": "isolated_cannabinoid",
                "intended_use": "Epilepsy treatment",
                "route": "oral",
                "target_markets": ["US"]
            }
        })
        assert response.status_code == 200
        
        data = response.json()
        assert "regpath_strategy_id" in data
        assert "primary_pathway" in data
        assert "readiness_checklist" in data
        assert "timeline" in data
    
    def test_generate_strategy_with_evidence(self, client):
        """Test strategy with evidence inputs."""
        response = client.post("/api/v1/regpath/strategy", json={
            "product_profile": {
                "product_name": "CBD-Plus",
                "product_type": "formulation",
                "intended_use": "Anxiety",
                "route": "oral"
            },
            "evidence_inputs": {
                "chempath_job_id": "chem-123",
                "predicate_drug": "Epidiolex"
            }
        })
        assert response.status_code == 200
        
        data = response.json()
        assert data["evidence_strength"] in ["moderate", "strong"]
    
    def test_get_strategy_by_id(self, client):
        """Test retrieving strategy by ID."""
        # First create a strategy
        create_response = client.post("/api/v1/regpath/strategy", json={
            "product_profile": {
                "product_name": "Test",
                "product_type": "formulation",
                "intended_use": "Test",
                "route": "oral"
            }
        })
        strategy_id = create_response.json()["regpath_strategy_id"]
        
        # Then retrieve it
        response = client.get(f"/api/v1/regpath/strategies/{strategy_id}")
        assert response.status_code == 200
        
        data = response.json()
        assert data["regpath_strategy_id"] == strategy_id
    
    def test_get_nonexistent_strategy(self, client):
        """Test 404 for nonexistent strategy."""
        response = client.get("/api/v1/regpath/strategies/nonexistent-id")
        assert response.status_code == 404
    
    def test_generate_memo_endpoint(self, client):
        """Test memo generation endpoint."""
        # First create a strategy
        create_response = client.post("/api/v1/regpath/strategy", json={
            "product_profile": {
                "product_name": "CBD Oral",
                "product_type": "isolated_cannabinoid",
                "intended_use": "Epilepsy",
                "route": "oral"
            }
        })
        strategy_id = create_response.json()["regpath_strategy_id"]
        
        # Generate memo
        response = client.post(
            f"/api/v1/regpath/strategies/{strategy_id}/memo",
            json={"format": "markdown"}
        )
        assert response.status_code == 200
        
        data = response.json()
        assert data["strategy_id"] == strategy_id
        assert "Regulatory Strategy Memo" in data["memo_content"]
    
    def test_generate_html_memo_endpoint(self, client):
        """Test HTML memo generation via API."""
        create_response = client.post("/api/v1/regpath/strategy", json={
            "product_profile": {
                "product_name": "Test",
                "product_type": "formulation",
                "intended_use": "Test",
                "route": "oral"
            }
        })
        strategy_id = create_response.json()["regpath_strategy_id"]
        
        response = client.post(
            f"/api/v1/regpath/strategies/{strategy_id}/memo",
            json={"format": "html"}
        )
        assert response.status_code == 200
        
        data = response.json()
        assert data["format"] == "html"
        assert "<!DOCTYPE html>" in data["memo_content"]


# ============================================================================
# Integration Tests
# ============================================================================

class TestRegPathIntegration:
    """Integration tests for RegPath."""
    
    def test_full_strategy_serialization(self):
        """Test full strategy can be serialized to dict."""
        from backend.services.regpath import (
            RegPathStrategist, RegPathRequest, ProductProfile, EvidenceInputs
        )
        
        strategist = RegPathStrategist()
        request = RegPathRequest(
            product_profile=ProductProfile(
                product_name="Full Test Product",
                product_type="isolated_cannabinoid",
                intended_use="Comprehensive test",
                route="oral",
                therapeutic_area="Neurology",
                active_ingredients=["Cannabidiol"]
            ),
            evidence_inputs=EvidenceInputs(
                chempath_job_id="chem-test",
                toxpath_assessment_id="tox-test"
            )
        )
        
        result = strategist.generate_strategy(request)
        data = result.to_dict()
        
        # All required fields present
        assert "regpath_strategy_id" in data
        assert "primary_pathway" in data
        assert "readiness_checklist" in data
        assert "timeline" in data
        assert "next_actions" in data
    
    def test_pipeline_integration(self):
        """Test ChemPath → ToxPath → RegPath conceptual pipeline."""
        from backend.services.regpath import (
            RegPathStrategist, RegPathRequest, ProductProfile, EvidenceInputs
        )
        
        # Simulate having completed ChemPath and ToxPath
        strategist = RegPathStrategist()
        request = RegPathRequest(
            product_profile=ProductProfile(
                product_name="Pipeline Product",
                product_type="novel_dimer",
                intended_use="Pain management",
                route="sublingual"
            ),
            evidence_inputs=EvidenceInputs(
                chempath_job_id="chempath-job-123",
                toxpath_assessment_id="toxpath-assess-456"
            )
        )
        
        result = strategist.generate_strategy(request)
        
        # Should be completed with evidence strength boosted
        assert result.status == "completed"
        assert result.evidence_strength in ["moderate", "strong"]
        
        # Novel dimer should get IND
        from backend.services.regpath import RegulatoryPathway
        assert result.primary_pathway.pathway == RegulatoryPathway.IND
    
    def test_full_memo_generation(self):
        """Test complete strategy → memo pipeline."""
        from backend.services.regpath import (
            RegPathStrategist, RegPathRequest, ProductProfile,
            RegPathMemoGenerator, MemoFormat, MemoConfig
        )
        
        # Step 1: Generate strategy
        strategist = RegPathStrategist()
        request = RegPathRequest(
            product_profile=ProductProfile(
                product_name="CBD-X Oral Solution",
                product_type="isolated_cannabinoid",
                intended_use="Treatment of drug-resistant epilepsy",
                route="oral",
                therapeutic_area="Neurology",
                active_ingredients=["Cannabidiol"]
            )
        )
        result = strategist.generate_strategy(request)
        
        # Step 2: Generate comprehensive memo
        config = MemoConfig(
            company_name="NeuroBotanica Research",
            prepared_by="Regulatory Affairs Team",
            include_cost_estimates=True,
            include_timeline_table=True
        )
        generator = RegPathMemoGenerator(config)
        memo = generator.generate(result, MemoFormat.MARKDOWN)
        
        # Verify memo completeness
        assert "CBD-X Oral Solution" in memo
        assert "Pathway Recommendation" in memo
        assert "Timeline" in memo
        assert "Checklist" in memo
        assert "NeuroBotanica Research" in memo


# ============================================================================
# Edge Case Tests
# ============================================================================

class TestEdgeCases:
    """Test edge cases and error handling."""
    
    def test_unknown_product_type(self):
        """Test handling of unknown product type."""
        from backend.services.regpath import (
            RegPathStrategist, RegPathRequest, ProductProfile
        )
        
        strategist = RegPathStrategist()
        request = RegPathRequest(
            product_profile=ProductProfile(
                product_name="Unknown Type Product",
                product_type="unknown_type",
                intended_use="Test",
                route="oral"
            )
        )
        
        result = strategist.generate_strategy(request)
        
        # Should still complete with default pathway
        assert result.status == "completed"
    
    def test_empty_intended_use(self):
        """Test handling of empty intended use."""
        from backend.services.regpath import (
            RegPathStrategist, RegPathRequest, ProductProfile
        )
        
        strategist = RegPathStrategist()
        request = RegPathRequest(
            product_profile=ProductProfile(
                product_name="Test",
                product_type="formulation",
                intended_use="",
                route="oral"
            )
        )
        
        result = strategist.generate_strategy(request)
        
        assert result.status == "completed"
    
    def test_multiple_target_markets(self):
        """Test handling of multiple target markets."""
        from backend.services.regpath import (
            RegPathStrategist, RegPathRequest, ProductProfile
        )
        
        strategist = RegPathStrategist()
        request = RegPathRequest(
            product_profile=ProductProfile(
                product_name="Global Product",
                product_type="isolated_cannabinoid",
                intended_use="Test",
                route="oral",
                target_markets=["US", "EU", "Canada"]
            )
        )
        
        result = strategist.generate_strategy(request)
        
        assert result.status == "completed"
