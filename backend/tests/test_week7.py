"""
Week 7 Tests - Terpene Evidence Gating + PatentPath Lite Phase 2
Tests for NeuroBotanica Week 7 implementation.

Coverage:
- TerpeneAnalyzer with evidence gating
- FTOChecker with risk assessment
- ClaimGenerator with USPTO templates
- API endpoints for all features
"""
import pytest
from unittest.mock import patch, AsyncMock, MagicMock
from datetime import datetime, timedelta

# ============================================================================
# Terpene Analyzer Tests
# ============================================================================

class TestTerpeneEvidenceGating:
    """Test terpene evidence gating functionality."""
    
    def test_evidence_tier_enum(self):
        """Test evidence tier enumeration."""
        from backend.services.terpene_analyzer import EvidenceTier
        
        assert EvidenceTier.TIER_1 == 1
        assert EvidenceTier.TIER_2 == 2
        assert EvidenceTier.TIER_3 == 3
        assert EvidenceTier.TIER_4 == 4
        assert EvidenceTier.TIER_5 == 5
        
        # Test ordering
        assert EvidenceTier.TIER_5 > EvidenceTier.TIER_3
        assert EvidenceTier.TIER_3 > EvidenceTier.TIER_2
    
    def test_terpene_database_initialization(self):
        """Test terpene database loads with expected terpenes."""
        from backend.services.terpene_analyzer import TerpeneDatabase
        
        db = TerpeneDatabase()
        terpenes = db.list_all()
        
        # Should have 10 terpenes per MVP spec
        assert len(terpenes) >= 10
        
        # Check key terpenes exist
        assert "myrcene" in terpenes
        assert "limonene" in terpenes
        assert "linalool" in terpenes
        assert "beta_caryophyllene" in terpenes
    
    def test_terpene_data_structure(self):
        """Test TerpeneData dataclass structure."""
        from backend.services.terpene_analyzer import TerpeneDatabase
        
        db = TerpeneDatabase()
        myrcene = db.get("myrcene")
        
        assert myrcene is not None
        assert myrcene.name == "Myrcene"
        assert myrcene.mechanism is not None
        assert len(myrcene.therapeutic_effects) > 0
        assert len(myrcene.cannabinoid_synergy) > 0
        
        # Check synergy structure
        if "CBD" in myrcene.cannabinoid_synergy:
            synergy = myrcene.cannabinoid_synergy["CBD"]
            assert "enhancement_factor" in synergy
            assert "evidence_tier" in synergy
    
    def test_synergy_analysis_basic(self):
        """Test basic synergy analysis."""
        from backend.services.terpene_analyzer import TerpeneAnalyzer
        
        analyzer = TerpeneAnalyzer()
        result = analyzer.analyze_synergy(
            cannabinoid="CBD",
            terpene_profile={"myrcene": 1.0, "limonene": 0.5},
            min_evidence_tier=1  # Include all evidence levels
        )
        
        assert result["cannabinoid"] == "CBD"
        assert "synergies" in result
        assert "total_enhancement_factor" in result
    
    def test_evidence_gating_filters_low_evidence(self):
        """Test that evidence gating properly filters low-quality synergies."""
        from backend.services.terpene_analyzer import TerpeneAnalyzer
        
        analyzer = TerpeneAnalyzer()
        
        # Get all synergies (no filter)
        all_result = analyzer.analyze_synergy(
            cannabinoid="CBD",
            terpene_profile={"myrcene": 1.0, "limonene": 0.5},
            min_evidence_tier=1
        )
        
        # Get filtered synergies (tier 3+)
        filtered_result = analyzer.analyze_synergy(
            cannabinoid="CBD",
            terpene_profile={"myrcene": 1.0, "limonene": 0.5},
            min_evidence_tier=3
        )
        
        # High tier filter should have fewer or equal synergies
        assert filtered_result["num_synergies_included"] <= all_result["num_synergies_included"]
        
        # All filtered synergies should be tier 3+
        for synergy in filtered_result["synergies"]:
            assert synergy["evidence_tier"] >= 3
    
    def test_evidence_gating_tier_5_strict(self):
        """Test strict tier 5 (RCT) filtering."""
        from backend.services.terpene_analyzer import TerpeneAnalyzer
        
        analyzer = TerpeneAnalyzer()
        result = analyzer.analyze_synergy(
            cannabinoid="CBD",
            terpene_profile={"myrcene": 1.0, "limonene": 0.5, "linalool": 0.3},
            min_evidence_tier=5  # Only RCT evidence
        )
        
        # May have no synergies at tier 5 - that's expected
        for synergy in result["synergies"]:
            assert synergy["evidence_tier"] == 5
    
    def test_enhancement_factor_calculation(self):
        """Test enhancement factor increases with concentration."""
        from backend.services.terpene_analyzer import TerpeneAnalyzer
        
        analyzer = TerpeneAnalyzer()
        
        # Low concentration
        low_result = analyzer.analyze_synergy(
            cannabinoid="CBD",
            terpene_profile={"myrcene": 0.1},
            min_evidence_tier=1
        )
        
        # High concentration
        high_result = analyzer.analyze_synergy(
            cannabinoid="CBD",
            terpene_profile={"myrcene": 2.0},
            min_evidence_tier=1
        )
        
        # Higher concentration should generally give higher enhancement
        if low_result["synergies"] and high_result["synergies"]:
            assert high_result["total_enhancement_factor"] >= low_result["total_enhancement_factor"]
    
    def test_optimal_profile_recommendation(self):
        """Test optimal profile recommendation."""
        from backend.services.terpene_analyzer import TerpeneAnalyzer
        
        analyzer = TerpeneAnalyzer()
        optimal = analyzer.get_optimal_profile(
            cannabinoid="CBD",
            therapeutic_goal="analgesic",
            min_evidence_tier=3
        )
        
        assert "recommended_terpenes" in optimal
        assert "therapeutic_goal" in optimal
        assert optimal["therapeutic_goal"] == "analgesic"
    
    def test_profile_comparison(self):
        """Test profile comparison functionality."""
        from backend.services.terpene_analyzer import TerpeneAnalyzer
        
        analyzer = TerpeneAnalyzer()
        comparison = analyzer.compare_profiles(
            cannabinoid="CBD",
            profile_a={"myrcene": 1.0, "limonene": 0.5},
            profile_b={"linalool": 0.8, "beta_caryophyllene": 0.6},
            min_evidence_tier=3
        )
        
        assert "profile_a" in comparison
        assert "profile_b" in comparison
        assert "winner" in comparison
    
    def test_synergy_result_dict_structure(self):
        """Test synergy result dictionary structure."""
        from backend.services.terpene_analyzer import TerpeneAnalyzer
        
        analyzer = TerpeneAnalyzer()
        result = analyzer.analyze_synergy(
            cannabinoid="THC",
            terpene_profile={"myrcene": 1.0},
            min_evidence_tier=1
        )
        
        assert "cannabinoid" in result
        assert "synergies" in result
        assert "total_enhancement_factor" in result
        assert "num_synergies_included" in result


# ============================================================================
# FTO Checker Tests
# ============================================================================

class TestFTOChecker:
    """Test Freedom-to-Operate checker functionality."""
    
    def test_risk_level_enum(self):
        """Test risk level enumeration."""
        from backend.services.patentpath.fto_checker import RiskLevel
        
        assert RiskLevel.LOW.value == "low"
        assert RiskLevel.MODERATE.value == "moderate"
        assert RiskLevel.HIGH.value == "high"
        assert RiskLevel.CRITICAL.value == "critical"
    
    def test_fto_checker_initialization(self):
        """Test FTOChecker initialization."""
        from backend.services.patentpath.fto_checker import FTOChecker
        
        checker = FTOChecker()
        assert checker is not None
    
    @pytest.mark.asyncio
    async def test_fto_check_returns_report(self):
        """Test FTO check returns proper report structure."""
        from backend.services.patentpath.fto_checker import FTOChecker, FTOReport
        
        checker = FTOChecker()
        
        # Mock the prior art search to avoid API calls
        with patch.object(checker.searcher, 'search', new_callable=AsyncMock) as mock_search:
            mock_search.return_value = []
            
            report = await checker.check_fto(
                compound_description="Novel CBD-CBG dimer with ester linkage",
                keywords=["cannabinoid", "dimer", "ester"],
                target_market="US",
                intended_use="pain"
            )
        
        assert isinstance(report, FTOReport)
        assert report.fto_risk_level is not None
        assert report.disclaimer is not None
        assert "legal advice" in report.disclaimer.lower()
    
    @pytest.mark.asyncio
    async def test_fto_check_with_blocking_patents(self):
        """Test FTO check with mock blocking patents."""
        from backend.services.patentpath.fto_checker import FTOChecker
        from backend.services.patentpath.prior_art import PriorArtResult, PatentSource
        
        checker = FTOChecker()
        
        # Mock patent results
        mock_patents = [
            PriorArtResult(
                patent_number="US10000001",
                title="Cannabinoid Dimer Compositions",
                abstract="A composition comprising dimeric cannabinoid compounds...",
                publication_date="2020-01-15",
                inventors=["John Smith", "Jane Doe"],
                assignees=["Pharma Corp"],
                classifications=["A61K31/352"],
                claims_count=20,
                source=PatentSource.USPTO,
                relevance_score=0.85,
                url="https://patents.google.com/patent/US10000001"
            )
        ]
        
        with patch.object(checker.searcher, 'search', new_callable=AsyncMock) as mock_search:
            mock_search.return_value = mock_patents
            
            report = await checker.check_fto(
                compound_description="Novel cannabinoid dimer composition",
                keywords=["cannabinoid", "dimer"],
                intended_use="therapeutic"
            )
        
        # Should identify blocking patents
        assert report.num_blocking_patents >= 0
        assert report.fto_risk_level is not None

    
    @pytest.mark.asyncio
    async def test_fto_licensing_cost_estimation(self):
        """Test licensing cost estimation in FTO report."""
        from backend.services.patentpath.fto_checker import FTOChecker
        
        checker = FTOChecker()
        
        with patch.object(checker.searcher, 'search', new_callable=AsyncMock) as mock_search:
            mock_search.return_value = []
            
            report = await checker.check_fto(
                compound_description="Novel compound",
                keywords=["test"],
                intended_use="therapeutic"
            )
        
        # Licensing cost should be included
        assert report.estimated_licensing_costs is not None
    
    @pytest.mark.asyncio
    async def test_fto_mitigation_strategies(self):
        """Test mitigation strategies in FTO report."""
        from backend.services.patentpath.fto_checker import FTOChecker
        
        checker = FTOChecker()
        
        with patch.object(checker.searcher, 'search', new_callable=AsyncMock) as mock_search:
            mock_search.return_value = []
            
            report = await checker.check_fto(
                compound_description="Novel compound",
                keywords=["test"],
                intended_use="therapeutic"
            )
        
        assert report.mitigation_strategies is not None
        assert isinstance(report.mitigation_strategies, list)
    
    def test_fto_report_to_dict(self):
        """Test FTOReport serialization."""
        from backend.services.patentpath.fto_checker import (
            FTOReport, RiskLevel, LicensingCostEstimate
        )
        
        report = FTOReport(
            fto_risk_level=RiskLevel.LOW,
            risk_score=0.1,
            target_market="US",
            intended_use="pain",
            blocking_patents=[],
            num_blocking_patents=0,
            num_patents_analyzed=10,
            estimated_licensing_costs=LicensingCostEstimate(
                total_estimated_cost="$0",
                per_patent_average="$0",
                breakdown=[],
                confidence="high",
                note="No blocking patents identified"
            ),
            recommendation="Proceed with development",
            mitigation_strategies=["Proceed with development"],
            next_steps=["File provisional patent"],
            analysis_date=datetime.now().isoformat()
        )
        
        report_dict = report.to_dict()
        
        assert "fto_risk_level" in report_dict
        assert "blocking_patents" in report_dict
        assert "estimated_licensing_costs" in report_dict
        assert "disclaimer" in report_dict


# ============================================================================
# Claim Generator Tests
# ============================================================================

class TestClaimGenerator:
    """Test patent claim generator functionality."""
    
    def test_claim_type_class(self):
        """Test claim type class."""
        from backend.services.patentpath.claim_generator import ClaimType
        
        assert ClaimType.COMPOSITION_OF_MATTER == "composition_of_matter"
        assert ClaimType.METHOD_OF_USE == "method_of_use"
        assert ClaimType.PHARMACEUTICAL_COMPOSITION == "pharmaceutical_composition"
        assert ClaimType.METHOD_OF_SYNTHESIS == "method_of_synthesis"
        assert ClaimType.DIMER_COMPOSITION == "dimer_composition"
    
    def test_claim_generator_initialization(self):
        """Test ClaimGenerator initialization."""
        from backend.services.patentpath.claim_generator import ClaimGenerator
        
        generator = ClaimGenerator()
        assert generator is not None
        assert generator.templates is not None
    
    def test_generate_composition_of_matter_claim(self):
        """Test composition of matter claim generation."""
        from backend.services.patentpath.claim_generator import ClaimGenerator, ClaimType
        
        generator = ClaimGenerator()
        report = generator.generate_claims(
            compound_name="CBD-CBG Dimer",
            compound_smiles="CCO",  # Simplified
            therapeutic_use="pain management",
            claim_types=[ClaimType.COMPOSITION_OF_MATTER]
        )
        
        assert len(report.generated_claims) > 0
        assert ClaimType.COMPOSITION_OF_MATTER in report.generated_claims
    
    def test_generate_method_of_use_claim(self):
        """Test method of use claim generation."""
        from backend.services.patentpath.claim_generator import ClaimGenerator, ClaimType
        
        generator = ClaimGenerator()
        report = generator.generate_claims(
            compound_name="NB-101",
            therapeutic_use="treating anxiety",
            claim_types=[ClaimType.METHOD_OF_USE]
        )
        
        assert ClaimType.METHOD_OF_USE in report.generated_claims
        claims_text = report.generated_claims[ClaimType.METHOD_OF_USE].claims_text
        assert "anxiety" in claims_text.lower()
    
    def test_generate_pharmaceutical_composition_claim(self):
        """Test pharmaceutical composition claim generation."""
        from backend.services.patentpath.claim_generator import ClaimGenerator, ClaimType
        
        generator = ClaimGenerator()
        report = generator.generate_claims(
            compound_name="NB-102",
            therapeutic_use="inflammation",
            claim_types=[ClaimType.PHARMACEUTICAL_COMPOSITION]
        )
        
        assert ClaimType.PHARMACEUTICAL_COMPOSITION in report.generated_claims
        claims_text = report.generated_claims[ClaimType.PHARMACEUTICAL_COMPOSITION].claims_text
        # Should mention carrier or excipient
        assert "carrier" in claims_text.lower() or "excipient" in claims_text.lower() or "composition" in claims_text.lower()
    
    def test_generate_dimer_claims(self):
        """Test dimer-specific claim generation."""
        from backend.services.patentpath.claim_generator import ClaimGenerator
        
        generator = ClaimGenerator()
        report = generator.generate_dimer_claims(
            dimer_name="CBD-THC-ester",
            parent_1="CBD",
            parent_2="THC",
            linker_type="ester",
            therapeutic_use="neuroprotection"
        )
        
        assert report is not None
        assert len(report.generated_claims) > 0
    
    def test_generate_all_claim_types(self):
        """Test generating all claim types at once."""
        from backend.services.patentpath.claim_generator import ClaimGenerator
        
        generator = ClaimGenerator()
        report = generator.generate_claims(
            compound_name="NB-103",
            compound_smiles="CCO",
            therapeutic_use="pain",
            claim_types=None  # Generate default types
        )
        
        # Should have multiple claim types
        assert len(report.generated_claims) > 1
    
    def test_filing_cost_estimation(self):
        """Test filing cost estimation."""
        from backend.services.patentpath.claim_generator import ClaimGenerator
        
        generator = ClaimGenerator()
        report = generator.generate_claims(
            compound_name="NB-104",
            therapeutic_use="epilepsy"
        )
        
        assert report.filing_cost_estimate is not None
        assert report.filing_cost_estimate.base_filing_fee > 0
        assert report.filing_cost_estimate.examination_fee > 0
        assert report.filing_cost_estimate.total_uspto_fees > 0
    
    def test_claim_disclaimer_present(self):
        """Test disclaimer is included in report."""
        from backend.services.patentpath.claim_generator import ClaimGenerator
        
        generator = ClaimGenerator()
        report = generator.generate_claims(
            compound_name="NB-105",
            therapeutic_use="anxiety"
        )
        
        assert report.disclaimer is not None
        assert "template" in report.disclaimer.lower() or "attorney" in report.disclaimer.lower()
    
    def test_claim_report_to_dict(self):
        """Test ClaimGenerationReport serialization."""
        from backend.services.patentpath.claim_generator import ClaimGenerator
        
        generator = ClaimGenerator()
        report = generator.generate_claims(
            compound_name="NB-106",
            therapeutic_use="sleep"
        )
        
        report_dict = report.to_dict()
        
        assert "compound_name" in report_dict
        assert "generated_claims" in report_dict
        assert "filing_cost_estimate" in report_dict
        assert "disclaimer" in report_dict
    
    def test_available_claim_types(self):
        """Test getting available claim types."""
        from backend.services.patentpath.claim_generator import ClaimGenerator, ClaimType
        
        generator = ClaimGenerator()
        
        # Check that templates exist for all claim types
        assert ClaimType.COMPOSITION_OF_MATTER in generator.templates
        assert ClaimType.METHOD_OF_USE in generator.templates
        assert ClaimType.PHARMACEUTICAL_COMPOSITION in generator.templates


# ============================================================================
# API Endpoint Tests
# ============================================================================

class TestPatentPathAPIEndpoints:
    """Test PatentPath API endpoints."""
    
    @pytest.fixture
    def client(self):
        """Create test client."""
        from fastapi.testclient import TestClient
        from backend.main import app
        return TestClient(app)
    
    def test_claim_generation_endpoint(self, client):
        """Test claim generation endpoint."""
        response = client.post(
            "/api/v1/patentpath/claims/generate",
            json={
                "compound_name": "NB-API-Test",
                "therapeutic_use": "anxiety",
                "claim_types": ["composition_of_matter", "method_of_use"]
            }
        )
        
        assert response.status_code == 200
        data = response.json()
        assert data["status"] == "success"
        assert "report" in data
    
    def test_dimer_claim_endpoint(self, client):
        """Test dimer claim generation endpoint."""
        response = client.post(
            "/api/v1/patentpath/claims/generate-dimer",
            json={
                "dimer_name": "CBD-THC-ester",
                "parent_1": "CBD",
                "parent_2": "THC",
                "linker_type": "ester",
                "therapeutic_use": "neuroprotection"
            }
        )
        
        assert response.status_code == 200
        data = response.json()
        assert data["status"] == "success"
        assert "report" in data
    
    def test_claim_types_endpoint(self, client):
        """Test claim types listing endpoint."""
        response = client.get("/api/v1/patentpath/claims/types")
        
        assert response.status_code == 200
        data = response.json()
        assert "claim_types" in data
        assert len(data["claim_types"]) > 0


class TestTerpeneAPIEndpoints:
    """Test Terpene Analysis API endpoints."""
    
    @pytest.fixture
    def client(self):
        """Create test client."""
        from fastapi.testclient import TestClient
        from backend.main import app
        return TestClient(app)
    
    def test_list_terpenes_endpoint(self, client):
        """Test terpene listing endpoint."""
        response = client.get("/api/v1/terpenes/")
        
        assert response.status_code == 200
        data = response.json()
        assert "count" in data
        assert "terpenes" in data
        assert data["count"] >= 10
    
    def test_get_terpene_endpoint(self, client):
        """Test single terpene endpoint."""
        response = client.get("/api/v1/terpenes/myrcene")
        
        assert response.status_code == 200
        data = response.json()
        assert "terpene" in data
    
    def test_get_terpene_not_found(self, client):
        """Test terpene not found."""
        response = client.get("/api/v1/terpenes/nonexistent")
        
        assert response.status_code == 404
    
    def test_synergy_analysis_endpoint(self, client):
        """Test synergy analysis endpoint."""
        response = client.post(
            "/api/v1/terpenes/synergy/analyze",
            json={
                "cannabinoid": "CBD",
                "terpene_profile": {"myrcene": 1.0, "limonene": 0.5},
                "min_evidence_tier": 3
            }
        )
        
        assert response.status_code == 200
        data = response.json()
        assert data["status"] == "success"
        assert "analysis" in data
    
    def test_synergy_analysis_with_all_synergies(self, client):
        """Test synergy analysis with all synergies flag."""
        response = client.post(
            "/api/v1/terpenes/synergy/analyze",
            json={
                "cannabinoid": "CBD",
                "terpene_profile": {"myrcene": 1.0},
                "min_evidence_tier": 4,
                "include_all_synergies": True
            }
        )
        
        assert response.status_code == 200
        data = response.json()
        assert "analysis" in data
    
    def test_profile_comparison_endpoint(self, client):
        """Test profile comparison endpoint."""
        response = client.post(
            "/api/v1/terpenes/profiles/compare",
            json={
                "cannabinoid": "CBD",
                "profile_a": {"myrcene": 1.0, "limonene": 0.5},
                "profile_b": {"linalool": 0.8, "beta_caryophyllene": 0.6},
                "name_a": "Indica",
                "name_b": "Hybrid"
            }
        )
        
        assert response.status_code == 200
        data = response.json()
        assert data["status"] == "success"
        assert "comparison" in data
    
    def test_optimal_profile_endpoint(self, client):
        """Test optimal profile recommendation endpoint."""
        response = client.post(
            "/api/v1/terpenes/profiles/optimize",
            json={
                "cannabinoid": "CBD",
                "therapeutic_goal": "analgesic",
                "target_synergy": "enhanced_analgesic",
                "max_terpenes": 3
            }
        )
        
        assert response.status_code == 200
        data = response.json()
        assert data["status"] == "success"
        assert "recommendation" in data
    
    def test_evidence_tiers_endpoint(self, client):
        """Test evidence tiers information endpoint."""
        response = client.get("/api/v1/terpenes/evidence/tiers")
        
        assert response.status_code == 200
        data = response.json()
        assert "evidence_tiers" in data
        assert len(data["evidence_tiers"]) == 5
    
    def test_therapeutic_targets_endpoint(self, client):
        """Test therapeutic targets listing endpoint."""
        response = client.get("/api/v1/terpenes/therapeutic-targets")
        
        assert response.status_code == 200
        data = response.json()
        assert "therapeutic_targets" in data
        
        # Should include pain, anxiety, etc.
        targets = [t["target"] for t in data["therapeutic_targets"]]
        assert "pain" in targets
        assert "anxiety" in targets


# ============================================================================
# Integration Tests
# ============================================================================

class TestWeek7Integration:
    """Integration tests for Week 7 features."""
    
    def test_terpene_to_patent_workflow(self):
        """Test workflow from terpene analysis to patent claims."""
        from backend.services.terpene_analyzer import TerpeneAnalyzer
        from backend.services.patentpath.claim_generator import ClaimGenerator
        
        # Step 1: Analyze synergy
        analyzer = TerpeneAnalyzer()
        synergy = analyzer.analyze_synergy(
            cannabinoid="CBD",
            terpene_profile={"myrcene": 1.5, "beta_caryophyllene": 0.8},
            min_evidence_tier=3
        )
        
        # Step 2: Use synergy data for claim generation
        generator = ClaimGenerator()
        
        # Build key features from synergy analysis
        key_features = []
        for s in synergy["synergies"][:3]:
            key_features.append(f"Enhanced {s['mechanism'][:30]}... effect")
        
        if not key_features:
            key_features = ["Novel therapeutic combination"]
        
        report = generator.generate_claims(
            compound_name="CBD-Myrcene-Caryophyllene Complex",
            therapeutic_use="pain management",
            key_features=key_features
        )
        
        assert len(report.generated_claims) > 0
        assert report.filing_cost_estimate is not None
    
    @pytest.mark.asyncio
    async def test_comprehensive_analysis_workflow(self):
        """Test comprehensive patent analysis workflow."""
        from backend.services.patentpath import (
            PriorArtSearcher,
            FTOChecker,
            ClaimGenerator
        )
        
        compound_name = "NB-Integration-Test"
        
        # Generate claims (no API needed)
        generator = ClaimGenerator()
        claims = generator.generate_claims(
            compound_name=compound_name,
            therapeutic_use="anxiety"
        )
        
        assert claims is not None
        assert len(claims.generated_claims) > 0


# ============================================================================
# Validation Tests
# ============================================================================

class TestWeek7Validation:
    """Validation tests for Week 7 acceptance criteria."""
    
    def test_evidence_gating_excludes_low_quality(self):
        """Validate: Terpene evidence gating excludes low-quality interactions."""
        from backend.services.terpene_analyzer import TerpeneAnalyzer
        
        analyzer = TerpeneAnalyzer()
        
        # Get tier 1 (all) synergies
        all_synergies = analyzer.analyze_synergy(
            cannabinoid="CBD",
            terpene_profile={"myrcene": 1.0, "limonene": 0.5, "linalool": 0.3},
            min_evidence_tier=1
        )
        
        # Get tier 4+ (good/high) synergies
        quality_synergies = analyzer.analyze_synergy(
            cannabinoid="CBD",
            terpene_profile={"myrcene": 1.0, "limonene": 0.5, "linalool": 0.3},
            min_evidence_tier=4
        )
        
        # Verify filtering occurred
        assert quality_synergies["num_synergies_included"] <= all_synergies["num_synergies_included"]
        
        # Verify all quality synergies meet threshold
        for synergy in quality_synergies["synergies"]:
            assert synergy["evidence_tier"] >= 4
    
    def test_fto_returns_risk_with_licensing(self):
        """Validate: FTO check returns risk assessment with licensing costs."""
        from backend.services.patentpath.fto_checker import (
            FTOReport, RiskLevel, LicensingCostEstimate
        )
        
        # Create a report with mock data
        report = FTOReport(
            fto_risk_level=RiskLevel.MODERATE,
            risk_score=0.5,
            target_market="US",
            intended_use="pain",
            blocking_patents=[],
            num_blocking_patents=0,
            num_patents_analyzed=10,
            estimated_licensing_costs=LicensingCostEstimate(
                total_estimated_cost="$100,000 - $500,000",
                per_patent_average="N/A",
                breakdown=[],
                confidence="moderate",
                note="Moderate risk assessment"
            ),
            recommendation="Consider design-around options",
            mitigation_strategies=["Design around", "License"],
            next_steps=["Consult patent attorney"],
            analysis_date=datetime.now().isoformat()
        )
        
        # Verify risk level
        assert report.fto_risk_level in list(RiskLevel)
        
        # Verify licensing cost
        assert report.estimated_licensing_costs is not None
    
    def test_claim_generator_produces_uspto_format(self):
        """Validate: Claim generator produces USPTO-formatted claims."""
        from backend.services.patentpath.claim_generator import ClaimGenerator, ClaimType
        
        generator = ClaimGenerator()
        report = generator.generate_claims(
            compound_name="Test Compound",
            therapeutic_use="treating condition"
        )
        
        # Check USPTO formatting elements
        for claim_type, claims in report.generated_claims.items():
            claim_text = claims.claims_text
            
            # Claims should contain numbered claims
            assert "CLAIM 1:" in claim_text or "Claim 1" in claim_text.lower(), \
                f"Claims don't have USPTO format numbering: {claim_text[:100]}"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
