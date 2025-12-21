"""
Test Suite for TKPath - Traditional Knowledge Attribution & Compensation

Tests the complete TKPath module including:
- Attribution calculation
- Compensation distribution  
- Certificate generation
- Verification lab
- Ethical pricing integration
"""

import pytest
from datetime import datetime, timedelta
from unittest.mock import patch

import sys
sys.path.insert(0, str(__file__).replace('/tests/test_tkpath.py', '').replace('\\tests\\test_tkpath.py', ''))

from backend.services.tkpath.attribution import (
    TKPathAttribution,
    TKContribution,
    CommunityAttribution,
    KnowledgeElement,
    KnowledgeElementType,
)
from backend.services.tkpath.compensation import (
    TKCompensationEngine,
    CompensationTransaction,
    CompensationStatus,
    PaymentMethod,
)
from backend.services.tkpath.certificate import (
    TKCertificateGenerator,
    AttributionCertificate,
    CertificateStatus,
    CertificateType,
)
from backend.services.tkpath.verification import (
    TKVerificationLab,
    VerificationReport,
    VerificationStatus,
    ChemicalProfile,
    ProvenanceRecord,
    ProvenanceStatus,
)


class TestTKPathAttribution:
    """Tests for TK attribution calculation."""
    
    def setup_method(self):
        """Set up test fixtures."""
        self.attribution = TKPathAttribution()
    
    def test_calculate_basic_attribution(self):
        """Test basic attribution calculation."""
        result = self.attribution.calculate_tk_contribution(
            formulation="Traditional Full-Spectrum CBD Extract",
            compounds=["CBD", "THC", "CBG"],
        )
        
        assert result is not None
        assert isinstance(result, TKContribution)
        assert result.formulation_name == "Traditional Full-Spectrum CBD Extract"
        assert result.verification_hash is not None
    
    def test_detect_knowledge_elements(self):
        """Test detection of TK elements in formulation."""
        result = self.attribution.calculate_tk_contribution(
            formulation="Traditional entourage effect full-spectrum pain relief tincture",
        )
        
        assert result.knowledge_element_count > 0
        # Should detect "traditional", "entourage", "full-spectrum", "pain", "tincture"
    
    def test_community_matching(self):
        """Test matching of TK elements to communities."""
        result = self.attribution.calculate_tk_contribution(
            formulation="Ayurvedic-style CBD preparation for anxiety",
            claimed_communities=["indian_ayurvedic"],
        )
        
        assert len(result.communities) > 0
        assert any("Indian" in c.community_name for c in result.communities)
    
    def test_compensation_suggestions(self):
        """Test compensation suggestion calculation."""
        result = self.attribution.calculate_tk_contribution(
            formulation="Traditional preparation with multiple TK elements: "
                       "entourage effect synergy pain relief microdose sublingual",
        )
        
        assert result.suggested_compensation_percentage > 0
        assert result.suggested_minimum_payment >= 25.0
    
    def test_list_known_communities(self):
        """Test listing of known communities."""
        communities = self.attribution.list_known_communities()
        
        assert len(communities) > 0
        assert all("id" in c for c in communities)
        assert all("name" in c for c in communities)
        assert all("country" in c for c in communities)
    
    def test_verification_hash_generation(self):
        """Test that verification hash is generated."""
        result = self.attribution.calculate_tk_contribution(
            formulation="Test formulation",
        )
        
        assert result.verification_hash is not None
        assert len(result.verification_hash) == 64  # SHA-256 hex


class TestTKCompensationEngine:
    """Tests for TK compensation distribution."""
    
    def setup_method(self):
        """Set up test fixtures."""
        self.compensation = TKCompensationEngine()
        self.attribution = TKPathAttribution()
    
    def test_distribute_compensation(self):
        """Test compensation distribution."""
        # Create attribution first
        attribution_result = self.attribution.calculate_tk_contribution(
            formulation="Traditional CBD preparation",
            claimed_communities=["indian_ayurvedic"],
        )
        
        # Distribute compensation
        transaction = self.compensation.distribute_compensation(
            revenue=1000.0,
            attribution=attribution_result,
        )
        
        assert transaction is not None
        assert isinstance(transaction, CompensationTransaction)
        assert transaction.total_amount > 0
        assert transaction.status == CompensationStatus.PENDING
    
    def test_process_surcharge_contribution(self):
        """Test surcharge contribution to pool."""
        result = self.compensation.process_surcharge_contribution(
            company_id="test-company-001",
            surcharge_amount=500.0,
        )
        
        assert result["success"] is True
        assert result["amount_added"] == 500.0
        assert result["pool_balance"] > 0
    
    def test_pool_accumulation(self):
        """Test that surcharge pool accumulates correctly."""
        initial_status = self.compensation.get_pool_status()
        initial_balance = initial_status["current_balance"]
        
        # Add multiple surcharges
        self.compensation.process_surcharge_contribution("company-a", 100.0)
        self.compensation.process_surcharge_contribution("company-b", 200.0)
        self.compensation.process_surcharge_contribution("company-c", 300.0)
        
        final_status = self.compensation.get_pool_status()
        
        assert final_status["current_balance"] == initial_balance + 600.0
        assert final_status["contributing_companies"] >= 3
    
    def test_company_contribution_report(self):
        """Test company contribution report generation."""
        # Add contribution
        self.compensation.process_surcharge_contribution("audit-test-company", 250.0)
        
        report = self.compensation.get_company_contribution_report("audit-test-company")
        
        assert report["company_id"] == "audit-test-company"
        assert report["total_surcharge_contributions"] == 250.0
        assert "message" in report
    
    def test_transaction_approval_flow(self):
        """Test transaction approval workflow."""
        # Create attribution and transaction
        attribution_result = self.attribution.calculate_tk_contribution(
            formulation="Test formulation",
        )
        transaction = self.compensation.distribute_compensation(
            revenue=500.0,
            attribution=attribution_result,
        )
        
        # Approve transaction
        approval_result = self.compensation.approve_transaction(
            transaction_id=transaction.transaction_id,
            approver_id="admin-001",
        )
        
        assert approval_result["success"] is True
        assert approval_result["status"] == "approved"
        
        # Complete transaction
        complete_result = self.compensation.complete_transaction(
            transaction_id=transaction.transaction_id,
            blockchain_tx_hash="0x1234567890abcdef",
        )
        
        assert complete_result["success"] is True
        assert complete_result["status"] == "completed"


class TestTKCertificateGenerator:
    """Tests for TK certificate generation."""
    
    def setup_method(self):
        """Set up test fixtures."""
        self.attribution = TKPathAttribution()
        self.certificate_gen = TKCertificateGenerator()
    
    def test_generate_basic_certificate(self):
        """Test basic certificate generation."""
        attribution_result = self.attribution.calculate_tk_contribution(
            formulation="Traditional CBD Oil",
            claimed_communities=["indian_ayurvedic"],
        )
        
        certificate = self.certificate_gen.generate_attribution_certificate(
            formulation="Traditional CBD Oil",
            attribution=attribution_result,
        )
        
        assert certificate is not None
        assert isinstance(certificate, AttributionCertificate)
        assert certificate.status == CertificateStatus.ACTIVE
        assert certificate.verification_url is not None
    
    def test_certificate_types(self):
        """Test different certificate types."""
        attribution_result = self.attribution.calculate_tk_contribution(
            formulation="Premium Formulation",
        )
        
        # Generate compensation-verified certificate
        cert = self.certificate_gen.generate_attribution_certificate(
            formulation="Premium Formulation",
            attribution=attribution_result,
            certificate_type=CertificateType.COMPENSATION_VERIFIED,
            compensation_verified=True,
            total_compensated=500.0,
        )
        
        assert cert.certificate_type in [
            CertificateType.COMPENSATION_VERIFIED,
            CertificateType.NAGOYA_COMPLIANT,
        ]
        assert cert.total_compensated == 500.0
    
    def test_badge_text_generation(self):
        """Test display badge text generation."""
        attribution_result = self.attribution.calculate_tk_contribution(
            formulation="Test Product",
            claimed_communities=["indian_ayurvedic", "chinese_tcm"],
        )
        
        certificate = self.certificate_gen.generate_attribution_certificate(
            formulation="Test Product",
            attribution=attribution_result,
        )
        
        assert certificate.display_badge_text != ""
        assert "TK" in certificate.display_badge_text or "Community" in certificate.display_badge_text
    
    def test_certificate_verification(self):
        """Test certificate verification."""
        attribution_result = self.attribution.calculate_tk_contribution(
            formulation="Verifiable Product",
        )
        
        certificate = self.certificate_gen.generate_attribution_certificate(
            formulation="Verifiable Product",
            attribution=attribution_result,
        )
        
        # Verify certificate
        verification = self.certificate_gen.verify_certificate(
            certificate_id=certificate.certificate_id,
        )
        
        assert verification["valid"] is True
        assert verification["certificate_id"] == certificate.certificate_id
    
    def test_certificate_expiration(self):
        """Test certificate validity checking."""
        attribution_result = self.attribution.calculate_tk_contribution(
            formulation="Expiring Product",
        )
        
        certificate = self.certificate_gen.generate_attribution_certificate(
            formulation="Expiring Product",
            attribution=attribution_result,
        )
        
        assert certificate.is_valid() is True
        assert certificate.expires_at > datetime.utcnow()
    
    def test_certificate_revocation(self):
        """Test certificate revocation."""
        attribution_result = self.attribution.calculate_tk_contribution(
            formulation="Revocable Product",
        )
        
        certificate = self.certificate_gen.generate_attribution_certificate(
            formulation="Revocable Product",
            attribution=attribution_result,
        )
        
        # Revoke certificate
        result = self.certificate_gen.revoke_certificate(
            certificate_id=certificate.certificate_id,
            reason="Compliance violation",
            revoked_by="admin-001",
        )
        
        assert result["success"] is True
        assert result["status"] == "revoked"
        
        # Verify now fails
        verification = self.certificate_gen.verify_certificate(
            certificate_id=certificate.certificate_id,
        )
        
        assert verification["valid"] is False


class TestTKVerificationLab:
    """Tests for TK verification lab."""
    
    def setup_method(self):
        """Set up test fixtures."""
        self.lab = TKVerificationLab()
    
    def test_list_traditional_preparations(self):
        """Test listing available traditional preparations."""
        preparations = self.lab.list_traditional_preparations()
        
        assert len(preparations) > 0
        assert all("key" in p for p in preparations)
        assert all("name" in p for p in preparations)
    
    def test_verify_preparation_with_profile(self):
        """Test verification with chemical profile."""
        # Create chemical profile matching Himalayan Charras
        profile = ChemicalProfile(
            cannabinoids={"THC": 20.0, "CBD": 1.5, "CBN": 0.5},
            terpenes={"myrcene": 0.5, "limonene": 0.3, "pinene": 0.2},
        )
        
        report = self.lab.verify_traditional_preparation(
            claimed_formulation="Premium Himalayan Hash",
            traditional_recipe="charras_himalayan",
            chemical_profile=profile,
        )
        
        assert report is not None
        assert isinstance(report, VerificationReport)
        assert report.chemical_match_score > 0
    
    def test_verify_with_provenance(self):
        """Test verification with provenance record."""
        provenance = ProvenanceRecord(
            source_community="Hindu Kush Traditional Cultivators",
            source_region="Hindu Kush Mountains",
            source_country="Afghanistan",
            cultivation_method="traditional_outdoor",
            status=ProvenanceStatus.BLOCKCHAIN_VERIFIED,
        )
        
        report = self.lab.verify_traditional_preparation(
            claimed_formulation="Authentic Charras",
            traditional_recipe="charras_himalayan",
            provenance=provenance,
        )
        
        assert report.provenance_score > 0
    
    def test_verification_status_levels(self):
        """Test different verification status outcomes."""
        # Good match
        good_profile = ChemicalProfile(
            cannabinoids={"THC": 18.0, "CBD": 1.8, "CBN": 0.3},
            terpenes={"myrcene": 0.5, "limonene": 0.35, "pinene": 0.25},
        )
        
        good_report = self.lab.verify_traditional_preparation(
            claimed_formulation="Well-matched Product",
            traditional_recipe="charras_himalayan",
            chemical_profile=good_profile,
        )
        
        # Bad match
        bad_profile = ChemicalProfile(
            cannabinoids={"THC": 5.0, "CBD": 20.0},  # Way off
            terpenes={"linalool": 2.0},  # Wrong terpenes
        )
        
        bad_report = self.lab.verify_traditional_preparation(
            claimed_formulation="Poorly-matched Product",
            traditional_recipe="charras_himalayan",
            chemical_profile=bad_profile,
        )
        
        assert good_report.overall_match_score > bad_report.overall_match_score
    
    def test_record_provenance(self):
        """Test provenance recording."""
        provenance = ProvenanceRecord(
            source_community="Test Community",
            source_region="Test Region",
            source_country="Test Country",
            cultivation_method="organic",
            harvest_date=datetime.utcnow(),
            gps_coordinates=(35.0, 70.0),
        )
        
        result = self.lab.record_provenance(provenance)
        
        assert result["success"] is True
        assert result["record_id"] is not None
        assert result["blockchain_hash"] is not None


class TestEthicalPricingIntegration:
    """Tests for integration between TKPath and ethical pricing."""
    
    def test_surcharge_funds_compensation(self):
        """Test that surcharges fund TK compensation."""
        compensation = TKCompensationEngine()
        attribution_service = TKPathAttribution()
        
        # Add surcharge (simulating company opting out of TK)
        compensation.process_surcharge_contribution(
            company_id="big-pharma-corp",
            surcharge_amount=1000.0,
        )
        
        # Verify pool has funds
        pool_status = compensation.get_pool_status()
        assert pool_status["current_balance"] >= 1000.0
        
        # These funds would be distributed to communities
        # via monthly distribution
    
    def test_compensation_tracks_surcharge_funding(self):
        """Test that transactions track whether funded by surcharges."""
        compensation = TKCompensationEngine()
        attribution_service = TKPathAttribution()
        
        attribution = attribution_service.calculate_tk_contribution(
            formulation="Test Product",
        )
        
        # Transaction funded by regular subscription
        regular_tx = compensation.distribute_compensation(
            revenue=500.0,
            attribution=attribution,
            revenue_source="subscription",
            surcharge_funded=False,
        )
        
        # Transaction funded by surcharge pool
        surcharge_tx = compensation.distribute_compensation(
            revenue=500.0,
            attribution=attribution,
            revenue_source="surcharge",
            surcharge_funded=True,
        )
        
        assert regular_tx.surcharge_funded_percentage == 0.0
        assert surcharge_tx.surcharge_funded_percentage == 100.0


class TestTKPathEndToEnd:
    """End-to-end tests for complete TKPath workflow."""
    
    def test_full_attribution_to_certificate_flow(self):
        """Test complete flow from attribution to certificate."""
        # 1. Calculate attribution
        attribution = TKPathAttribution()
        attribution_result = attribution.calculate_tk_contribution(
            formulation="Premium Ayurvedic CBD Formula with traditional entourage synergy",
            compounds=["CBD", "THC", "CBG", "myrcene", "limonene"],
            claimed_communities=["indian_ayurvedic"],
        )
        
        assert attribution_result.total_tk_contribution_percentage > 0
        
        # 2. Create compensation transaction
        compensation = TKCompensationEngine()
        transaction = compensation.distribute_compensation(
            revenue=1000.0,
            attribution=attribution_result,
        )
        
        assert transaction.total_amount > 0
        
        # 3. Approve and complete transaction
        compensation.approve_transaction(transaction.transaction_id, "admin")
        compensation.complete_transaction(
            transaction.transaction_id,
            "0xblockchain_hash_here",
        )
        
        # 4. Generate certificate
        cert_gen = TKCertificateGenerator()
        certificate = cert_gen.generate_attribution_certificate(
            formulation="Premium Ayurvedic CBD Formula",
            attribution=attribution_result,
            compensation_verified=True,
            total_compensated=transaction.total_amount,
            compensation_tx_ids=[transaction.transaction_id],
        )
        
        assert certificate.compensation_verified is True
        assert certificate.is_valid() is True
        
        # 5. Verify certificate
        verification = cert_gen.verify_certificate(certificate.certificate_id)
        assert verification["valid"] is True
        assert verification["compensation_verified"] is True


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
