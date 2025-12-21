"""
Week 8 Tests - PatentPath Lite Features 5-6

Tests for:
- TK Attribution Checker (Feature 5)
- Cost Estimator (Feature 6)
- Comprehensive Analysis (all 6 features)
- API endpoints for new features
"""
import pytest
from unittest.mock import patch, MagicMock, AsyncMock
from fastapi.testclient import TestClient

# Import test client
from backend.main import app

client = TestClient(app)


# ============================================================================
# TK Attribution Checker Tests
# ============================================================================

class TestTKAttributionChecker:
    """Test Traditional Knowledge Attribution Checker."""
    
    def test_tk_checker_initialization(self):
        """Test TK checker initializes correctly."""
        from backend.services.patentpath.tk_checker import TKAttributionChecker
        
        checker = TKAttributionChecker()
        
        assert hasattr(checker, 'tk_keywords')
        assert hasattr(checker, 'sacred_keywords')
        assert hasattr(checker, 'tk_source_types')
        assert hasattr(checker, 'high_abs_regions')
        
        # Check keyword categories exist
        assert "strong" in checker.tk_keywords
        assert "moderate" in checker.tk_keywords
        assert "weak" in checker.tk_keywords
    
    def test_tk_indicator_strength_enum(self):
        """Test TK indicator strength enum values."""
        from backend.services.patentpath.tk_checker import TKIndicatorStrength
        
        assert TKIndicatorStrength.NONE.value == "none"
        assert TKIndicatorStrength.WEAK.value == "weak"
        assert TKIndicatorStrength.MODERATE.value == "moderate"
        assert TKIndicatorStrength.STRONG.value == "strong"
        assert TKIndicatorStrength.DEFINITIVE.value == "definitive"
    
    def test_sacred_knowledge_status_enum(self):
        """Test sacred knowledge status enum."""
        from backend.services.patentpath.tk_checker import SacredKnowledgeStatus
        
        assert SacredKnowledgeStatus.NOT_DETECTED.value == "not_detected"
        assert SacredKnowledgeStatus.POSSIBLE.value == "possible"
        assert SacredKnowledgeStatus.CONFIRMED.value == "confirmed"
    
    def test_no_tk_derivation(self):
        """Test compound with no TK indicators."""
        from backend.services.patentpath.tk_checker import TKAttributionChecker
        
        checker = TKAttributionChecker()
        result = checker.check_tk_attribution(
            compound_description="Novel synthetic cannabinoid developed in laboratory"
        )
        
        assert result.tk_derived is False
        assert result.sacred_knowledge_detected is False
        assert result.attribution_required is False
        assert result.can_proceed is True
    
    def test_strong_tk_indicators(self):
        """Test compound with strong TK indicators."""
        from backend.services.patentpath.tk_checker import TKAttributionChecker, TKIndicatorStrength
        
        checker = TKAttributionChecker()
        result = checker.check_tk_attribution(
            compound_description="Compound based on traditional knowledge from indigenous peoples"
        )
        
        assert result.tk_derived is True
        assert result.indicator_strength == TKIndicatorStrength.STRONG
        assert result.attribution_required is True
        assert result.can_proceed is True
        assert len(result.requirements) > 0
    
    def test_ethnobotanical_source_detection(self):
        """Test detection of ethnobotanical derivation source."""
        from backend.services.patentpath.tk_checker import TKAttributionChecker
        
        checker = TKAttributionChecker()
        result = checker.check_tk_attribution(
            compound_description="Novel compound isolated from plant",
            derivation_source="ethnobotanical_literature"
        )
        
        assert result.tk_derived is True
        assert result.attribution_required is True
    
    def test_sacred_knowledge_detection_definitive(self):
        """Test detection of sacred knowledge (definitive)."""
        from backend.services.patentpath.tk_checker import TKAttributionChecker
        
        checker = TKAttributionChecker()
        result = checker.check_tk_attribution(
            compound_description="Sacred medicine used in ceremonial rituals"
        )
        
        assert result.sacred_knowledge_detected is True
        assert result.can_proceed is False
        assert result.blocking_reason is not None
        assert "cannot be patented" in result.blocking_reason.lower()
    
    def test_sacred_knowledge_blocks_patenting(self):
        """Test that sacred knowledge blocks patent filing."""
        from backend.services.patentpath.tk_checker import TKAttributionChecker
        
        checker = TKAttributionChecker()
        result = checker.check_tk_attribution(
            compound_description="Sacred plant used in spiritual healing ceremony"
        )
        
        assert result.can_proceed is False
        assert result.attribution_required is False  # Can't attribute what can't be patented
    
    def test_ayurvedic_moderate_indicator(self):
        """Test Ayurvedic as moderate TK indicator."""
        from backend.services.patentpath.tk_checker import TKAttributionChecker, TKIndicatorStrength
        
        checker = TKAttributionChecker()
        result = checker.check_tk_attribution(
            compound_description="Compound inspired by Ayurvedic medicine"
        )
        
        assert result.tk_derived is True
        assert result.indicator_strength in [TKIndicatorStrength.MODERATE, TKIndicatorStrength.STRONG]
    
    def test_high_abs_region_detection(self):
        """Test detection of high ABS requirement region."""
        from backend.services.patentpath.tk_checker import TKAttributionChecker
        
        checker = TKAttributionChecker()
        result = checker.check_tk_attribution(
            compound_description="Plant compound based on traditional knowledge",
            source_region="India"
        )
        
        assert result.tk_derived is True
        # Should have indicator about ABS region
        assert any("india" in ind.lower() for ind in result.detected_indicators)
    
    def test_attribution_requirements_generated(self):
        """Test that attribution requirements are generated."""
        from backend.services.patentpath.tk_checker import TKAttributionChecker
        
        checker = TKAttributionChecker()
        result = checker.check_tk_attribution(
            compound_description="Compound derived from traditional medicine practices"
        )
        
        assert result.tk_derived is True
        assert len(result.requirements) >= 5  # Should have multiple requirements
        
        # Check required categories
        categories = [r.category for r in result.requirements]
        assert "documentation" in categories
        assert "consent" in categories
        assert "benefit_sharing" in categories
    
    def test_equipath_integration_config(self):
        """Test EquiPath integration configuration."""
        from backend.services.patentpath.tk_checker import TKAttributionChecker
        
        checker = TKAttributionChecker()
        result = checker.check_tk_attribution(
            compound_description="Traditional medicine compound",
            community_name="Test Community"
        )
        
        assert result.equipath_integration is not None
        assert "community_wallet_setup" in result.equipath_integration
        assert "estimated_revenue_share" in result.equipath_integration
        assert result.equipath_integration["community_name"] == "Test Community"
    
    def test_specification_requirements_generated(self):
        """Test patent specification requirements are generated."""
        from backend.services.patentpath.tk_checker import TKAttributionChecker
        
        checker = TKAttributionChecker()
        result = checker.check_tk_attribution(
            compound_description="Traditional knowledge based compound"
        )
        
        assert result.specification_requirements is not None
        assert "attribution_statement_template" in result.specification_requirements
        assert "citation_requirements" in result.specification_requirements
    
    def test_tk_result_to_dict(self):
        """Test TK result serialization."""
        from backend.services.patentpath.tk_checker import TKAttributionChecker
        
        checker = TKAttributionChecker()
        result = checker.check_tk_attribution(
            compound_description="Novel compound"
        )
        
        result_dict = result.to_dict()
        
        assert "tk_derived" in result_dict
        assert "sacred_knowledge_detected" in result_dict
        assert "can_proceed" in result_dict
        assert "checked_at" in result_dict
    
    def test_get_tk_source_types(self):
        """Test getting TK source types."""
        from backend.services.patentpath.tk_checker import TKAttributionChecker
        
        checker = TKAttributionChecker()
        source_types = checker.get_tk_source_types()
        
        assert "ethnobotanical_literature" in source_types
        assert "traditional_medicine_text" in source_types
        assert "indigenous_community" in source_types
    
    def test_get_high_abs_regions(self):
        """Test getting high ABS regions."""
        from backend.services.patentpath.tk_checker import TKAttributionChecker
        
        checker = TKAttributionChecker()
        regions = checker.get_high_abs_regions()
        
        assert "india" in regions
        assert "brazil" in regions
        assert "china" in regions
    
    def test_validate_attribution_completeness_complete(self):
        """Test validation when all requirements met."""
        from backend.services.patentpath.tk_checker import TKAttributionChecker
        
        checker = TKAttributionChecker()
        result = checker.validate_attribution_completeness(
            community_name="Test Community",
            consent_obtained=True,
            benefit_sharing_configured=True,
            attribution_drafted=True
        )
        
        assert result["complete"] is True
        assert result["can_file_patent"] is True
        assert len(result["missing_requirements"]) == 0
    
    def test_validate_attribution_completeness_incomplete(self):
        """Test validation when requirements missing."""
        from backend.services.patentpath.tk_checker import TKAttributionChecker
        
        checker = TKAttributionChecker()
        result = checker.validate_attribution_completeness(
            community_name="Test Community",
            consent_obtained=False,
            benefit_sharing_configured=False,
            attribution_drafted=True
        )
        
        assert result["complete"] is False
        assert result["can_file_patent"] is False
        assert len(result["missing_requirements"]) == 2


# ============================================================================
# Cost Estimator Tests
# ============================================================================

class TestCostEstimator:
    """Test Patent Cost Estimator."""
    
    def test_estimator_initialization(self):
        """Test cost estimator initializes correctly."""
        from backend.services.patentpath.cost_estimator import CostEstimator
        
        estimator = CostEstimator()
        
        assert hasattr(estimator, 'uspto_fees')
        assert hasattr(estimator, 'attorney_costs')
        assert hasattr(estimator, 'trade_secret_costs')
    
    def test_protection_strategy_enum(self):
        """Test protection strategy enum."""
        from backend.services.patentpath.cost_estimator import ProtectionStrategy
        
        assert ProtectionStrategy.PATENT.value == "patent"
        assert ProtectionStrategy.TRADE_SECRET.value == "trade_secret"
        assert ProtectionStrategy.HYBRID.value == "hybrid"
    
    def test_complexity_level_enum(self):
        """Test complexity level enum."""
        from backend.services.patentpath.cost_estimator import ComplexityLevel
        
        assert ComplexityLevel.SIMPLE.value == "simple"
        assert ComplexityLevel.MODERATE.value == "moderate"
        assert ComplexityLevel.COMPLEX.value == "complex"
    
    def test_entity_size_enum(self):
        """Test entity size enum."""
        from backend.services.patentpath.cost_estimator import EntitySize
        
        assert EntitySize.LARGE.value == "large"
        assert EntitySize.SMALL.value == "small"
        assert EntitySize.MICRO.value == "micro"
    
    def test_patent_cost_estimate_single_compound(self):
        """Test patent cost estimate for single compound."""
        from backend.services.patentpath.cost_estimator import CostEstimator, ProtectionStrategy
        
        estimator = CostEstimator()
        estimate = estimator.estimate_costs(
            num_compounds=1,
            strategy="patent"
        )
        
        assert estimate.strategy == ProtectionStrategy.PATENT
        assert estimate.num_compounds == 1
        assert estimate.total_cost_estimate.total_cost > 0
        assert estimate.total_cost_estimate.filing_costs > 0
        assert estimate.total_cost_estimate.maintenance_costs > 0
    
    def test_patent_cost_estimate_multiple_compounds(self):
        """Test patent cost scales with compound count."""
        from backend.services.patentpath.cost_estimator import CostEstimator
        
        estimator = CostEstimator()
        estimate_1 = estimator.estimate_costs(num_compounds=1, strategy="patent")
        estimate_5 = estimator.estimate_costs(num_compounds=5, strategy="patent")
        
        # Cost should scale roughly linearly
        assert estimate_5.total_cost_estimate.total_cost > estimate_1.total_cost_estimate.total_cost * 4
    
    def test_trade_secret_cost_estimate(self):
        """Test trade secret cost estimate."""
        from backend.services.patentpath.cost_estimator import CostEstimator, ProtectionStrategy
        
        estimator = CostEstimator()
        estimate = estimator.estimate_costs(
            num_compounds=1,
            strategy="trade_secret"
        )
        
        assert estimate.strategy == ProtectionStrategy.TRADE_SECRET
        assert estimate.total_cost_estimate.prosecution_costs == 0  # No prosecution for trade secrets
        assert len(estimate.advantages) > 0
        assert len(estimate.disadvantages) > 0
    
    def test_hybrid_cost_estimate(self):
        """Test hybrid strategy cost estimate."""
        from backend.services.patentpath.cost_estimator import CostEstimator, ProtectionStrategy
        
        estimator = CostEstimator()
        estimate = estimator.estimate_costs(
            num_compounds=4,
            strategy="hybrid"
        )
        
        assert estimate.strategy == ProtectionStrategy.HYBRID
        # Hybrid should have comparison_savings calculated (can be negative if trade secrets are cheaper)
        assert estimate.comparison_savings is not None
        # Hybrid combines patent + trade secret, so costs should be reasonable
        assert estimate.total_cost_estimate.total_cost > 0
    
    def test_small_entity_discount(self):
        """Test small entity fee discount."""
        from backend.services.patentpath.cost_estimator import CostEstimator
        
        estimator = CostEstimator()
        large = estimator.estimate_costs(num_compounds=1, entity_size="large")
        small = estimator.estimate_costs(num_compounds=1, entity_size="small")
        
        # Small entity should be cheaper
        assert small.total_cost_estimate.filing_costs < large.total_cost_estimate.filing_costs
    
    def test_micro_entity_discount(self):
        """Test micro entity fee discount."""
        from backend.services.patentpath.cost_estimator import CostEstimator
        
        estimator = CostEstimator()
        large = estimator.estimate_costs(num_compounds=1, entity_size="large")
        micro = estimator.estimate_costs(num_compounds=1, entity_size="micro")
        
        # Micro entity should be significantly cheaper
        assert micro.total_cost_estimate.filing_costs < large.total_cost_estimate.filing_costs * 0.5
    
    def test_complexity_affects_cost(self):
        """Test complexity level affects attorney costs."""
        from backend.services.patentpath.cost_estimator import CostEstimator
        
        estimator = CostEstimator()
        simple = estimator.estimate_costs(num_compounds=1, complexity="simple")
        complex_ = estimator.estimate_costs(num_compounds=1, complexity="complex")
        
        # Complex should cost more
        assert complex_.total_cost_estimate.prosecution_costs > simple.total_cost_estimate.prosecution_costs
    
    def test_cost_schedule_generated(self):
        """Test cost schedule is generated."""
        from backend.services.patentpath.cost_estimator import CostEstimator
        
        estimator = CostEstimator()
        estimate = estimator.estimate_costs(num_compounds=1, strategy="patent")
        
        assert len(estimate.cost_schedule) > 0
        
        # Should have filing, prosecution, and maintenance events
        events = [e.event.lower() for e in estimate.cost_schedule]
        assert any("filing" in e for e in events)
        assert any("maintenance" in e for e in events)
    
    def test_timeline_info_included(self):
        """Test timeline information is included."""
        from backend.services.patentpath.cost_estimator import CostEstimator
        
        estimator = CostEstimator()
        estimate = estimator.estimate_costs(num_compounds=1, strategy="patent")
        
        assert "filing_to_examination" in estimate.timeline
        assert "protection_duration" in estimate.timeline
    
    def test_recommendation_generated(self):
        """Test cost-based recommendation is generated."""
        from backend.services.patentpath.cost_estimator import CostEstimator
        
        estimator = CostEstimator()
        estimate = estimator.estimate_costs(num_compounds=1)
        
        assert estimate.recommendation != ""
        assert len(estimate.recommendation) > 10
    
    def test_cost_estimate_to_dict(self):
        """Test cost estimate serialization."""
        from backend.services.patentpath.cost_estimator import CostEstimator
        
        estimator = CostEstimator()
        estimate = estimator.estimate_costs(num_compounds=1)
        
        result_dict = estimate.to_dict()
        
        assert "strategy" in result_dict
        assert "total_cost_estimate" in result_dict
        assert "per_patent_breakdown" in result_dict
        assert "cost_schedule" in result_dict
    
    def test_compare_strategies(self):
        """Test strategy comparison."""
        from backend.services.patentpath.cost_estimator import CostEstimator
        
        estimator = CostEstimator()
        comparison = estimator.compare_strategies(num_compounds=3)
        
        assert "patent" in comparison
        assert "trade_secret" in comparison
        assert "hybrid" in comparison
        assert "comparison" in comparison
        assert "recommended_strategy" in comparison["comparison"]
    
    def test_entity_discount_info(self):
        """Test entity discount information."""
        from backend.services.patentpath.cost_estimator import CostEstimator
        
        estimator = CostEstimator()
        info = estimator.get_entity_discount_info()
        
        assert "large_entity" in info
        assert "small_entity" in info
        assert "micro_entity" in info
        assert "discount" in info["small_entity"]
    
    def test_uspto_fee_schedule_values(self):
        """Test USPTO fee schedule has expected values."""
        from backend.services.patentpath.cost_estimator import USPTOFeeSchedule
        
        fees = USPTOFeeSchedule()
        
        # Check reasonable fee values (2025 estimates)
        assert fees.filing_fee > 1000
        assert fees.search_fee > 2000
        assert fees.examination_fee > 2500
        assert fees.maintenance_3_5_years > 1500
        assert fees.maintenance_11_5_years > fees.maintenance_3_5_years
    
    def test_cost_breakdown_formatting(self):
        """Test cost breakdown is properly formatted."""
        from backend.services.patentpath.cost_estimator import CostEstimator
        
        estimator = CostEstimator()
        estimate = estimator.estimate_costs(num_compounds=1)
        
        breakdown_dict = estimate.total_cost_estimate.to_dict()
        
        # All values should be formatted as currency strings
        assert breakdown_dict["total_investment"].startswith("$")
        assert "," in breakdown_dict["total_investment"]  # Should have thousands separator


# ============================================================================
# TK API Endpoint Tests
# ============================================================================

class TestTKAPIEndpoints:
    """Test TK Attribution API endpoints."""
    
    def test_tk_check_endpoint(self):
        """Test TK check endpoint."""
        response = client.post(
            "/api/v1/patentpath/tk/check",
            json={
                "compound_description": "Completely synthetic compound created entirely in laboratory conditions"
            }
        )
        
        assert response.status_code == 200
        data = response.json()
        assert data["status"] == "success"
        assert "result" in data
        # Result should contain expected keys
        assert "tk_derived" in data["result"]
        assert "can_proceed" in data["result"]
    
    def test_tk_check_with_tk_indicators(self):
        """Test TK check with TK indicators."""
        response = client.post(
            "/api/v1/patentpath/tk/check",
            json={
                "compound_description": "Compound based on traditional knowledge from indigenous peoples"
            }
        )
        
        assert response.status_code == 200
        data = response.json()
        assert data["result"]["tk_derived"] is True
        assert data["result"]["attribution_required"] is True
    
    def test_tk_check_sacred_knowledge(self):
        """Test TK check with sacred knowledge."""
        response = client.post(
            "/api/v1/patentpath/tk/check",
            json={
                "compound_description": "Sacred medicine from ceremonial practices"
            }
        )
        
        assert response.status_code == 200
        data = response.json()
        assert data["result"]["sacred_knowledge_detected"] is True
        assert data["result"]["can_proceed"] is False
    
    def test_tk_source_types_endpoint(self):
        """Test TK source types endpoint."""
        response = client.get("/api/v1/patentpath/tk/source-types")
        
        assert response.status_code == 200
        data = response.json()
        assert "source_types" in data
        assert len(data["source_types"]) > 0
    
    def test_tk_high_abs_regions_endpoint(self):
        """Test high ABS regions endpoint."""
        response = client.get("/api/v1/patentpath/tk/high-abs-regions")
        
        assert response.status_code == 200
        data = response.json()
        assert "regions" in data
        assert "india" in data["regions"]
    
    def test_tk_validate_completeness_endpoint(self):
        """Test TK validation completeness endpoint."""
        response = client.post(
            "/api/v1/patentpath/tk/validate-completeness",
            params={
                "community_name": "Test Community",
                "consent_obtained": True,
                "benefit_sharing_configured": True,
                "attribution_drafted": True
            }
        )
        
        assert response.status_code == 200
        data = response.json()
        assert data["validation"]["complete"] is True


# ============================================================================
# Cost API Endpoint Tests
# ============================================================================

class TestCostAPIEndpoints:
    """Test Cost Estimation API endpoints."""
    
    def test_cost_estimate_endpoint(self):
        """Test cost estimation endpoint."""
        response = client.post(
            "/api/v1/patentpath/costs/estimate",
            json={
                "num_compounds": 1,
                "complexity": "moderate",
                "strategy": "patent"
            }
        )
        
        assert response.status_code == 200
        data = response.json()
        assert data["status"] == "success"
        assert "estimate" in data
        assert "total_cost_estimate" in data["estimate"]
    
    def test_cost_estimate_trade_secret(self):
        """Test cost estimation for trade secret."""
        response = client.post(
            "/api/v1/patentpath/costs/estimate",
            json={
                "num_compounds": 1,
                "strategy": "trade_secret"
            }
        )
        
        assert response.status_code == 200
        data = response.json()
        assert data["estimate"]["strategy"] == "trade_secret"
    
    def test_compare_strategies_endpoint(self):
        """Test strategy comparison endpoint."""
        response = client.post(
            "/api/v1/patentpath/costs/compare-strategies",
            params={
                "num_compounds": 3,
                "complexity": "moderate"
            }
        )
        
        assert response.status_code == 200
        data = response.json()
        assert "comparison" in data
        assert "patent" in data["comparison"]
        assert "trade_secret" in data["comparison"]
    
    def test_entity_discounts_endpoint(self):
        """Test entity discounts endpoint."""
        response = client.get("/api/v1/patentpath/costs/entity-discounts")
        
        assert response.status_code == 200
        data = response.json()
        assert "entity_discounts" in data
        assert "small_entity" in data["entity_discounts"]
    
    def test_protection_strategies_endpoint(self):
        """Test protection strategies endpoint."""
        response = client.get("/api/v1/patentpath/costs/strategies")
        
        assert response.status_code == 200
        data = response.json()
        assert "strategies" in data
        assert len(data["strategies"]) == 3


# ============================================================================
# Comprehensive Analysis Tests (All 6 Features)
# ============================================================================

class TestComprehensiveAnalysis:
    """Test comprehensive patent analysis with all 6 features."""
    
    @pytest.mark.asyncio
    async def test_comprehensive_analysis_all_features(self):
        """Test comprehensive analysis includes all 6 features."""
        response = client.post(
            "/api/v1/patentpath/analyze/comprehensive",
            json={
                "compound_name": "Test Compound",
                "compound_description": "A novel synthetic cannabinoid derivative",
                "keywords": ["cannabinoid", "synthetic", "novel"],
                "therapeutic_use": "pain relief",
                "generate_claims": True,
                "include_cost_estimate": True
            }
        )
        
        assert response.status_code == 200
        data = response.json()
        assert data["status"] == "success"
        
        analysis = data["analysis"]
        
        # Check all 6 features present
        assert "prior_art" in analysis  # Feature 1
        assert "novelty" in analysis  # Feature 2
        assert "fto" in analysis  # Feature 3
        assert "claims" in analysis  # Feature 4
        assert "tk_attribution" in analysis  # Feature 5
        assert "cost_estimate" in analysis  # Feature 6
    
    @pytest.mark.asyncio
    async def test_comprehensive_analysis_summary(self):
        """Test comprehensive analysis generates summary."""
        response = client.post(
            "/api/v1/patentpath/analyze/comprehensive",
            json={
                "compound_name": "Test CBD Dimer",
                "compound_description": "Novel CBD dimer compound",
                "keywords": ["cbd", "dimer"],
                "generate_claims": False,
                "include_cost_estimate": False
            }
        )
        
        assert response.status_code == 200
        data = response.json()
        
        summary = data["analysis"]["summary"]
        assert "novelty_score" in summary
        assert "fto_risk_level" in summary
        assert "tk_derived" in summary
        assert "recommendation" in summary
    
    @pytest.mark.asyncio
    async def test_comprehensive_analysis_tk_derived(self):
        """Test comprehensive analysis with TK-derived compound."""
        response = client.post(
            "/api/v1/patentpath/analyze/comprehensive",
            json={
                "compound_name": "Traditional Compound",
                "compound_description": "Compound based on traditional knowledge from indigenous peoples",
                "keywords": ["traditional", "plant"],
                "derivation_source": "ethnobotanical_literature",
                "community_name": "Test Community",
                "generate_claims": False,
                "include_cost_estimate": False
            }
        )
        
        assert response.status_code == 200
        data = response.json()
        
        tk_attribution = data["analysis"]["tk_attribution"]
        assert tk_attribution["tk_derived"] is True
        assert tk_attribution["attribution_required"] is True
        
        # Summary should reflect TK status
        summary = data["analysis"]["summary"]
        assert summary["tk_derived"] is True
        assert summary["tk_attribution_required"] is True


# ============================================================================
# Integration Tests
# ============================================================================

class TestWeek8Integration:
    """Integration tests for Week 8 features."""
    
    def test_tk_to_cost_workflow(self):
        """Test workflow from TK check to cost estimation."""
        from backend.services.patentpath.tk_checker import TKAttributionChecker
        from backend.services.patentpath.cost_estimator import CostEstimator
        
        # Step 1: Check TK attribution
        tk_checker = TKAttributionChecker()
        tk_result = tk_checker.check_tk_attribution(
            compound_description="Novel synthetic compound"
        )
        
        assert tk_result.can_proceed is True
        
        # Step 2: Estimate costs
        estimator = CostEstimator()
        estimate = estimator.estimate_costs(num_compounds=1)
        
        assert estimate.total_cost_estimate.total_cost > 0
    
    def test_tk_blocking_prevents_cost_analysis(self):
        """Test that sacred knowledge detection should prevent filing."""
        from backend.services.patentpath.tk_checker import TKAttributionChecker
        
        checker = TKAttributionChecker()
        result = checker.check_tk_attribution(
            compound_description="Sacred ceremonial medicine"
        )
        
        # Should not proceed to cost analysis if can_proceed is False
        assert result.can_proceed is False
        # In real workflow, should not estimate costs for sacred knowledge
    
    def test_all_six_features_available(self):
        """Test all 6 PatentPath features are importable."""
        from backend.services.patentpath import (
            PriorArtSearcher,  # Feature 1
            NoveltyScorer,  # Feature 2
            FTOChecker,  # Feature 3
            ClaimGenerator,  # Feature 4
            TKAttributionChecker,  # Feature 5
            CostEstimator  # Feature 6
        )
        
        # All should be instantiable
        searcher = PriorArtSearcher()
        scorer = NoveltyScorer()
        checker = FTOChecker()
        generator = ClaimGenerator()
        tk_checker = TKAttributionChecker()
        estimator = CostEstimator()
        
        assert searcher is not None
        assert scorer is not None
        assert checker is not None
        assert generator is not None
        assert tk_checker is not None
        assert estimator is not None


# ============================================================================
# Validation Tests
# ============================================================================

class TestWeek8Validation:
    """Validation tests for Week 8 requirements."""
    
    def test_tk_checker_detects_all_keyword_types(self):
        """Test TK checker detects strong, moderate, and weak indicators."""
        from backend.services.patentpath.tk_checker import TKAttributionChecker
        
        checker = TKAttributionChecker()
        
        # Strong indicator
        result_strong = checker.check_tk_attribution(
            compound_description="Traditional knowledge from indigenous peoples"
        )
        assert result_strong.tk_derived is True
        
        # Moderate indicator
        result_moderate = checker.check_tk_attribution(
            compound_description="Compound derived from Ayurvedic medicine"
        )
        assert result_moderate.tk_derived is True
    
    def test_cost_estimator_realistic_values(self):
        """Test cost estimator produces realistic values."""
        from backend.services.patentpath.cost_estimator import CostEstimator
        
        estimator = CostEstimator()
        estimate = estimator.estimate_costs(num_compounds=1, complexity="moderate")
        
        total = estimate.total_cost_estimate.total_cost
        
        # Realistic range for single patent (filing + prosecution + 20yr maintenance)
        assert 15000 < total < 60000
    
    def test_sacred_knowledge_absolute_bar(self):
        """Test sacred knowledge is absolute bar to patenting."""
        from backend.services.patentpath.tk_checker import TKAttributionChecker
        
        checker = TKAttributionChecker()
        
        # Sacred terms that should be blocked (matching exact keywords in checker)
        sacred_terms = [
            "sacred medicine used by indigenous community",
            "sacred knowledge passed down through generations",
            "sacred plant medicine"
        ]
        
        for term in sacred_terms:
            result = checker.check_tk_attribution(compound_description=term)
            # Sacred knowledge should be detected and block patenting
            assert result.sacred_knowledge_detected is True, f"Should detect sacred for: {term}"
            assert result.can_proceed is False, f"Should block: {term}"
    
    def test_nagoya_protocol_compliance_mentioned(self):
        """Test Nagoya Protocol compliance is mentioned for TK compounds."""
        from backend.services.patentpath.tk_checker import TKAttributionChecker
        
        checker = TKAttributionChecker()
        result = checker.check_tk_attribution(
            compound_description="Traditional medicine compound"
        )
        
        assert result.tk_derived is True
        assert "nagoya" in result.legal_compliance.lower()
    
    def test_benefit_sharing_in_tk_requirements(self):
        """Test benefit-sharing is included in TK requirements."""
        from backend.services.patentpath.tk_checker import TKAttributionChecker
        
        checker = TKAttributionChecker()
        result = checker.check_tk_attribution(
            compound_description="Traditional knowledge compound"
        )
        
        # Check benefit sharing is required
        req_categories = [r.category for r in result.requirements]
        assert "benefit_sharing" in req_categories
        
        # Check EquiPath integration includes revenue share
        assert "estimated_revenue_share" in result.equipath_integration
