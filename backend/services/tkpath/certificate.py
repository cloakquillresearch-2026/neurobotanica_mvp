"""
TKPath Certificate Generator

Generates verifiable certificates for marketing/regulatory use.
Companies can display: "TK-Attributed Formulation - 3 Communities Compensated"

Trade Secret: Certificate generation algorithms, verification protocols.
"""

from enum import Enum
from typing import Dict, List, Optional
from pydantic import BaseModel, Field
from datetime import datetime, timedelta
import uuid
import hashlib
import base64

from .attribution import TKContribution, CommunityAttribution


class CertificateStatus(str, Enum):
    """Status of attribution certificate."""
    DRAFT = "draft"
    ACTIVE = "active"
    EXPIRED = "expired"
    REVOKED = "revoked"
    SUSPENDED = "suspended"


class CertificateType(str, Enum):
    """Types of attribution certificates."""
    BASIC = "basic"              # Standard attribution only
    COMPENSATION_VERIFIED = "compensation_verified"  # With payment verification
    NAGOYA_COMPLIANT = "nagoya_compliant"  # Full Nagoya Protocol compliance
    PREMIUM = "premium"          # All verifications + extended validity


class AttributionCertificate(BaseModel):
    """
    Verifiable certificate of TK attribution and compensation.
    Can be displayed in marketing materials and regulatory submissions.
    """
    certificate_id: str = Field(default_factory=lambda: str(uuid.uuid4()))
    certificate_type: CertificateType = CertificateType.BASIC
    
    # Formulation details
    formulation_id: str
    formulation_name: str
    attribution_id: str
    
    # Attribution summary
    communities_attributed: List[str] = Field(default_factory=list)
    community_count: int = 0
    total_tk_contribution_percentage: float = 0.0
    primary_knowledge_types: List[str] = Field(default_factory=list)
    
    # Compensation verification
    compensation_verified: bool = False
    total_compensated: float = 0.0
    compensation_transaction_ids: List[str] = Field(default_factory=list)
    
    # Compliance
    nagoya_compliant: bool = False
    undrip_compliant: bool = False
    
    # Certificate metadata
    status: CertificateStatus = CertificateStatus.DRAFT
    issued_at: Optional[datetime] = None
    expires_at: Optional[datetime] = None
    issued_by: str = "NeuroBotanica TKPath"
    
    # Verification
    verification_hash: Optional[str] = None
    verification_url: Optional[str] = None
    qr_code_data: Optional[str] = None
    
    # Display text
    display_badge_text: str = ""
    display_full_text: str = ""
    
    def generate_verification_hash(self) -> str:
        """Generate tamper-evident hash for certificate verification."""
        data = (
            f"{self.certificate_id}:{self.formulation_id}:{self.attribution_id}:"
            f"{','.join(self.communities_attributed)}:{self.total_compensated}"
        )
        self.verification_hash = hashlib.sha256(data.encode()).hexdigest()
        return self.verification_hash
    
    def is_valid(self) -> bool:
        """Check if certificate is currently valid."""
        if self.status != CertificateStatus.ACTIVE:
            return False
        if self.expires_at and datetime.utcnow() > self.expires_at:
            return False
        return True


class TKCertificateGenerator:
    """
    Generate verifiable certificates for marketing/regulatory use.
    
    Trade Secret: Certificate generation algorithms, verification URLs,
    badge text generation, QR code encoding.
    """
    
    # Certificate validity periods by type
    _VALIDITY_PERIODS = {
        CertificateType.BASIC: timedelta(days=365),
        CertificateType.COMPENSATION_VERIFIED: timedelta(days=365),
        CertificateType.NAGOYA_COMPLIANT: timedelta(days=730),  # 2 years
        CertificateType.PREMIUM: timedelta(days=1095),  # 3 years
    }
    
    # Badge text templates
    _BADGE_TEMPLATES = {
        CertificateType.BASIC: "TK-Attributed | {count} {communities}",
        CertificateType.COMPENSATION_VERIFIED: "TK-Compensated | {count} {communities} | ${amount}",
        CertificateType.NAGOYA_COMPLIANT: "Nagoya Certified | {count} {communities} Compensated",
        CertificateType.PREMIUM: "Premium TK Certified | Full Compliance | {count} {communities}",
    }
    
    def __init__(self, base_verification_url: str = "https://verify.neurobotanica.ai/tk"):
        self._base_url = base_verification_url
        self._certificates: Dict[str, AttributionCertificate] = {}
    
    def generate_attribution_certificate(
        self,
        formulation: str,
        attribution: TKContribution,
        certificate_type: CertificateType = CertificateType.BASIC,
        compensation_verified: bool = False,
        total_compensated: float = 0.0,
        compensation_tx_ids: Optional[List[str]] = None,
    ) -> AttributionCertificate:
        """
        Generate verifiable certificate for marketing/regulatory use.
        
        Args:
            formulation: Formulation name
            attribution: TKContribution with community details
            certificate_type: Level of certification
            compensation_verified: Whether payments have been verified
            total_compensated: Total amount compensated to communities
            compensation_tx_ids: Transaction IDs for compensation payments
            
        Returns:
            AttributionCertificate ready for use
        """
        certificate = AttributionCertificate(
            formulation_id=attribution.formulation_id,
            formulation_name=formulation,
            attribution_id=attribution.contribution_id,
            certificate_type=certificate_type,
        )
        
        # Extract community information
        certificate.communities_attributed = [
            c.community_name for c in attribution.communities
        ]
        certificate.community_count = len(attribution.communities)
        certificate.total_tk_contribution_percentage = attribution.total_tk_contribution_percentage
        
        # Extract knowledge types
        knowledge_types = set()
        for community in attribution.communities:
            for elem in community.knowledge_elements:
                knowledge_types.add(elem.element_type)
        certificate.primary_knowledge_types = list(knowledge_types)
        
        # Compensation details
        certificate.compensation_verified = compensation_verified
        certificate.total_compensated = total_compensated
        certificate.compensation_transaction_ids = compensation_tx_ids or []
        
        # Compliance flags
        certificate.nagoya_compliant = attribution.nagoya_compliant
        certificate.undrip_compliant = attribution.undrip_compliant
        
        # Upgrade certificate type if compliance allows
        if compensation_verified and certificate_type == CertificateType.BASIC:
            certificate.certificate_type = CertificateType.COMPENSATION_VERIFIED
        if attribution.nagoya_compliant and compensation_verified:
            certificate.certificate_type = CertificateType.NAGOYA_COMPLIANT
        
        # Set validity period
        certificate.issued_at = datetime.utcnow()
        validity = self._VALIDITY_PERIODS.get(
            certificate.certificate_type, timedelta(days=365)
        )
        certificate.expires_at = certificate.issued_at + validity
        certificate.status = CertificateStatus.ACTIVE
        
        # Generate verification data
        certificate.generate_verification_hash()
        certificate.verification_url = (
            f"{self._base_url}/{certificate.certificate_id}"
        )
        certificate.qr_code_data = self._generate_qr_data(certificate)
        
        # Generate display text
        certificate.display_badge_text = self._generate_badge_text(certificate)
        certificate.display_full_text = self._generate_full_text(certificate)
        
        # Store certificate
        self._certificates[certificate.certificate_id] = certificate
        
        return certificate
    
    def _generate_badge_text(self, certificate: AttributionCertificate) -> str:
        """Generate short badge text for marketing materials."""
        template = self._BADGE_TEMPLATES.get(
            certificate.certificate_type,
            self._BADGE_TEMPLATES[CertificateType.BASIC]
        )
        
        communities_word = "Community" if certificate.community_count == 1 else "Communities"
        
        return template.format(
            count=certificate.community_count,
            communities=communities_word,
            amount=f"{certificate.total_compensated:.0f}",
        )
    
    def _generate_full_text(self, certificate: AttributionCertificate) -> str:
        """Generate full certification text for regulatory documents."""
        lines = [
            f"TRADITIONAL KNOWLEDGE ATTRIBUTION CERTIFICATE",
            f"Certificate ID: {certificate.certificate_id}",
            f"",
            f"This certifies that the formulation '{certificate.formulation_name}' "
            f"incorporates traditional knowledge from the following indigenous communities:",
            "",
        ]
        
        for i, community in enumerate(certificate.communities_attributed, 1):
            lines.append(f"  {i}. {community}")
        
        lines.extend([
            "",
            f"Total TK Contribution: {certificate.total_tk_contribution_percentage:.1f}%",
            f"Knowledge Types: {', '.join(certificate.primary_knowledge_types)}",
            "",
        ])
        
        if certificate.compensation_verified:
            lines.extend([
                f"COMPENSATION VERIFIED",
                f"Total Compensated: ${certificate.total_compensated:.2f}",
                "",
            ])
        
        if certificate.nagoya_compliant:
            lines.append("✓ Nagoya Protocol Compliant")
        if certificate.undrip_compliant:
            lines.append("✓ UNDRIP Compliant (Articles 31 & 32)")
        
        lines.extend([
            "",
            f"Issued: {certificate.issued_at.strftime('%Y-%m-%d')}",
            f"Expires: {certificate.expires_at.strftime('%Y-%m-%d')}",
            f"Verify: {certificate.verification_url}",
            "",
            f"Issued by: {certificate.issued_by}",
        ])
        
        return "\n".join(lines)
    
    def _generate_qr_data(self, certificate: AttributionCertificate) -> str:
        """Generate QR code data for certificate verification."""
        # Compact verification payload
        data = {
            "id": certificate.certificate_id[:8],  # Short ID
            "v": certificate.verification_hash[:16],  # Short hash
            "c": certificate.community_count,
            "t": certificate.certificate_type.value[0],  # First letter
        }
        
        # Encode as base64 for QR efficiency
        import json
        return base64.urlsafe_b64encode(
            json.dumps(data).encode()
        ).decode()
    
    def verify_certificate(
        self,
        certificate_id: str,
        verification_hash: Optional[str] = None,
    ) -> Dict:
        """
        Verify a certificate is valid and unaltered.
        
        Returns verification status and certificate details.
        """
        certificate = self._certificates.get(certificate_id)
        
        if not certificate:
            return {
                "valid": False,
                "error": "Certificate not found",
                "certificate_id": certificate_id,
            }
        
        # Check hash if provided
        if verification_hash and verification_hash != certificate.verification_hash:
            return {
                "valid": False,
                "error": "Verification hash mismatch - certificate may be altered",
                "certificate_id": certificate_id,
            }
        
        # Check validity
        if not certificate.is_valid():
            return {
                "valid": False,
                "error": f"Certificate status: {certificate.status.value}",
                "certificate_id": certificate_id,
                "expired": certificate.expires_at < datetime.utcnow() if certificate.expires_at else False,
            }
        
        return {
            "valid": True,
            "certificate_id": certificate_id,
            "certificate_type": certificate.certificate_type.value,
            "formulation": certificate.formulation_name,
            "communities": certificate.communities_attributed,
            "community_count": certificate.community_count,
            "compensation_verified": certificate.compensation_verified,
            "total_compensated": certificate.total_compensated,
            "nagoya_compliant": certificate.nagoya_compliant,
            "undrip_compliant": certificate.undrip_compliant,
            "issued_at": certificate.issued_at.isoformat(),
            "expires_at": certificate.expires_at.isoformat(),
            "verification_hash": certificate.verification_hash,
        }
    
    def revoke_certificate(
        self,
        certificate_id: str,
        reason: str,
        revoked_by: str,
    ) -> Dict:
        """Revoke a certificate (e.g., for compliance violations)."""
        certificate = self._certificates.get(certificate_id)
        
        if not certificate:
            return {"success": False, "error": "Certificate not found"}
        
        certificate.status = CertificateStatus.REVOKED
        
        return {
            "success": True,
            "certificate_id": certificate_id,
            "status": certificate.status.value,
            "revoked_by": revoked_by,
            "reason": reason,
        }
    
    def get_certificate(self, certificate_id: str) -> Optional[AttributionCertificate]:
        """Retrieve certificate by ID."""
        return self._certificates.get(certificate_id)
    
    def list_certificates(
        self,
        formulation_id: Optional[str] = None,
        status: Optional[CertificateStatus] = None,
    ) -> List[AttributionCertificate]:
        """List certificates with optional filtering."""
        certificates = list(self._certificates.values())
        
        if formulation_id:
            certificates = [c for c in certificates if c.formulation_id == formulation_id]
        if status:
            certificates = [c for c in certificates if c.status == status]
        
        return sorted(certificates, key=lambda c: c.issued_at or datetime.min, reverse=True)
