"""
TKPath API Router

REST API endpoints for Traditional Knowledge attribution, compensation,
and verification services.

Endpoints:
- POST /api/v1/tkpath/attribute - Calculate TK attribution for formulation
- POST /api/v1/tkpath/compensate - Execute compensation transaction
- GET /api/v1/tkpath/certificate/{formulation_id} - Get attribution certificate
- GET /api/v1/tkpath/audit/{company_id} - Audit trail of TK payments
- POST /api/v1/tkpath/verify - Verify TK attribution claim
- GET /api/v1/tkpath/communities - List known indigenous communities
"""

from typing import Dict, List, Optional
from datetime import datetime
from fastapi import APIRouter, HTTPException, Query, Depends
from pydantic import BaseModel, Field

from backend.services.tkpath import (
    TKPathAttribution,
    TKCompensationEngine,
    TKCertificateGenerator,
    TKVerificationLab,
    CertificateType,
    ChemicalProfile,
    ProvenanceRecord,
    ProvenanceStatus,
)


router = APIRouter(prefix="/api/v1/tkpath", tags=["TKPath - Traditional Knowledge"])


# =============================================================================
# Request/Response Models
# =============================================================================

class AttributionRequest(BaseModel):
    """Request for TK attribution analysis."""
    formulation_name: str = Field(..., description="Name of the formulation")
    formulation_id: Optional[str] = Field(None, description="Unique formulation ID")
    compounds: Optional[List[str]] = Field(None, description="List of compounds in formulation")
    preparation_method: Optional[str] = Field(None, description="Description of preparation method")
    claimed_communities: Optional[List[str]] = Field(None, description="Claimed TK source communities")


class CompensationRequest(BaseModel):
    """Request for TK compensation distribution."""
    attribution_id: str = Field(..., description="Attribution ID to compensate")
    revenue_amount: float = Field(..., ge=0, description="Revenue to calculate compensation from")
    revenue_source: str = Field("subscription", description="Source of funds")
    surcharge_funded: bool = Field(False, description="Whether funded by TK opt-out surcharges")


class SurchargeContributionRequest(BaseModel):
    """Request to add surcharge to compensation pool."""
    company_id: str = Field(..., description="Company paying the surcharge")
    surcharge_amount: float = Field(..., gt=0, description="Surcharge amount")


class CertificateRequest(BaseModel):
    """Request for TK attribution certificate."""
    formulation_name: str
    attribution_id: str
    certificate_type: str = Field("basic", description="Certificate type: basic, compensation_verified, nagoya_compliant, premium")
    compensation_verified: bool = False
    total_compensated: float = 0.0
    compensation_tx_ids: Optional[List[str]] = None


class VerificationRequest(BaseModel):
    """Request for TK verification analysis."""
    claimed_formulation: str = Field(..., description="Formulation claiming TK attribution")
    traditional_recipe: str = Field(..., description="Traditional preparation to verify against")
    chemical_profile: Optional[Dict[str, Dict[str, float]]] = Field(None, description="Chemical analysis results")
    provenance: Optional[Dict] = Field(None, description="Supply chain provenance data")


class ProvenanceRecordRequest(BaseModel):
    """Request to record provenance information."""
    source_community: str
    source_region: str
    source_country: str
    gps_coordinates: Optional[List[float]] = None
    cultivation_method: str
    harvest_date: Optional[str] = None
    harvest_timing_notes: Optional[str] = None
    processing_location: Optional[str] = None
    processing_method: Optional[str] = None


# =============================================================================
# Service Instances
# =============================================================================

_attribution_service = TKPathAttribution()
_compensation_service = TKCompensationEngine()
_certificate_service = TKCertificateGenerator()
_verification_service = TKVerificationLab()


# =============================================================================
# Attribution Endpoints
# =============================================================================

@router.post("/attribute", response_model=Dict)
async def calculate_tk_attribution(request: AttributionRequest) -> Dict:
    """
    Calculate TK attribution for a formulation.
    
    Returns:
    - Community attribution percentages
    - Specific knowledge elements used
    - Suggested compensation amounts
    - Nagoya/UNDRIP compliance status
    """
    contribution = _attribution_service.calculate_tk_contribution(
        formulation=request.formulation_name,
        formulation_id=request.formulation_id,
        compounds=request.compounds,
        preparation_method=request.preparation_method,
        claimed_communities=request.claimed_communities,
    )
    
    return {
        "contribution_id": contribution.contribution_id,
        "formulation_name": contribution.formulation_name,
        "communities": [
            {
                "community_id": c.community_id,
                "community_name": c.community_name,
                "region": c.region,
                "country": c.country,
                "contribution_percentage": c.total_contribution_percentage,
                "knowledge_elements": [
                    {
                        "type": e.element_type,
                        "description": e.description,
                        "weight": e.contribution_weight,
                    }
                    for e in c.knowledge_elements
                ],
            }
            for c in contribution.communities
        ],
        "total_tk_contribution_percentage": contribution.total_tk_contribution_percentage,
        "primary_community": contribution.primary_community,
        "suggested_compensation": {
            "percentage_of_revenue": contribution.suggested_compensation_percentage,
            "minimum_monthly_payment": contribution.suggested_minimum_payment,
        },
        "compliance": {
            "nagoya_compliant": contribution.nagoya_compliant,
            "undrip_compliant": contribution.undrip_compliant,
        },
        "verification_hash": contribution.verification_hash,
    }


@router.get("/communities", response_model=List[Dict])
async def list_known_communities() -> List[Dict]:
    """
    List all known indigenous communities with TK contributions.
    
    Returns community IDs, names, regions, and primary knowledge types.
    """
    return _attribution_service.list_known_communities()


# =============================================================================
# Compensation Endpoints
# =============================================================================

@router.post("/compensate", response_model=Dict)
async def distribute_compensation(request: CompensationRequest) -> Dict:
    """
    Execute compensation distribution to TK holders.
    
    Calculates fair distribution based on attribution percentages
    and initiates blockchain-verified payments.
    """
    # Get attribution
    attribution = _attribution_service.get_attribution_audit(request.attribution_id)
    if not attribution:
        raise HTTPException(status_code=404, detail="Attribution not found")
    
    transaction = _compensation_service.distribute_compensation(
        revenue=request.revenue_amount,
        attribution=attribution,
        revenue_source=request.revenue_source,
        surcharge_funded=request.surcharge_funded,
    )
    
    return {
        "transaction_id": transaction.transaction_id,
        "attribution_id": transaction.attribution_id,
        "total_amount": transaction.total_amount,
        "status": transaction.status.value,
        "recipients": [
            {
                "community_id": r.community_id,
                "community_name": r.community_name,
                "amount": transaction.recipient_amounts.get(r.recipient_id, 0),
            }
            for r in transaction.recipients
        ],
        "revenue_source": transaction.revenue_source,
        "surcharge_funded_percentage": transaction.surcharge_funded_percentage,
        "verification_hash": transaction.verification_hash,
    }


@router.post("/surcharge/contribute", response_model=Dict)
async def add_surcharge_contribution(request: SurchargeContributionRequest) -> Dict:
    """
    Add TK opt-out surcharge to compensation pool.
    
    When companies opt OUT of TK features, their surcharge payments
    fund compensation to indigenous communities.
    """
    return _compensation_service.process_surcharge_contribution(
        company_id=request.company_id,
        surcharge_amount=request.surcharge_amount,
    )


@router.get("/pool/status", response_model=Dict)
async def get_compensation_pool_status() -> Dict:
    """
    Get current status of TK compensation pool.
    
    Shows balance, total collected, total distributed, and metrics.
    """
    return _compensation_service.get_pool_status()


@router.post("/pool/distribute", response_model=Dict)
async def execute_monthly_distribution() -> Dict:
    """
    Execute monthly distribution of surcharge pool to communities.
    
    Distributes accumulated surcharge payments to TK holders.
    """
    transactions = _compensation_service.execute_monthly_distribution()
    
    return {
        "success": True,
        "transactions_created": len(transactions),
        "total_distributed": sum(tx.total_amount for tx in transactions),
        "transaction_ids": [tx.transaction_id for tx in transactions],
    }


@router.get("/audit/{company_id}", response_model=Dict)
async def get_company_contribution_report(company_id: str) -> Dict:
    """
    Get audit report of TK contributions from a company.
    
    Shows how opt-out surcharges fund TK compensation.
    """
    return _compensation_service.get_company_contribution_report(company_id)


# =============================================================================
# Certificate Endpoints
# =============================================================================

@router.post("/certificate/generate", response_model=Dict)
async def generate_attribution_certificate(request: CertificateRequest) -> Dict:
    """
    Generate verifiable TK attribution certificate.
    
    Certificates can be displayed in marketing materials and
    regulatory submissions.
    """
    # Get attribution
    attribution = _attribution_service.get_attribution_audit(request.attribution_id)
    if not attribution:
        raise HTTPException(status_code=404, detail="Attribution not found")
    
    # Map certificate type
    cert_type_map = {
        "basic": CertificateType.BASIC,
        "compensation_verified": CertificateType.COMPENSATION_VERIFIED,
        "nagoya_compliant": CertificateType.NAGOYA_COMPLIANT,
        "premium": CertificateType.PREMIUM,
    }
    cert_type = cert_type_map.get(request.certificate_type, CertificateType.BASIC)
    
    certificate = _certificate_service.generate_attribution_certificate(
        formulation=request.formulation_name,
        attribution=attribution,
        certificate_type=cert_type,
        compensation_verified=request.compensation_verified,
        total_compensated=request.total_compensated,
        compensation_tx_ids=request.compensation_tx_ids,
    )
    
    return {
        "certificate_id": certificate.certificate_id,
        "certificate_type": certificate.certificate_type.value,
        "status": certificate.status.value,
        "formulation_name": certificate.formulation_name,
        "communities_attributed": certificate.communities_attributed,
        "community_count": certificate.community_count,
        "compensation_verified": certificate.compensation_verified,
        "total_compensated": certificate.total_compensated,
        "nagoya_compliant": certificate.nagoya_compliant,
        "undrip_compliant": certificate.undrip_compliant,
        "issued_at": certificate.issued_at.isoformat() if certificate.issued_at else None,
        "expires_at": certificate.expires_at.isoformat() if certificate.expires_at else None,
        "verification_url": certificate.verification_url,
        "qr_code_data": certificate.qr_code_data,
        "display_badge_text": certificate.display_badge_text,
        "display_full_text": certificate.display_full_text,
    }


@router.get("/certificate/{certificate_id}", response_model=Dict)
async def get_certificate(certificate_id: str) -> Dict:
    """Get certificate by ID."""
    certificate = _certificate_service.get_certificate(certificate_id)
    if not certificate:
        raise HTTPException(status_code=404, detail="Certificate not found")
    
    return {
        "certificate_id": certificate.certificate_id,
        "certificate_type": certificate.certificate_type.value,
        "status": certificate.status.value,
        "formulation_name": certificate.formulation_name,
        "communities_attributed": certificate.communities_attributed,
        "display_badge_text": certificate.display_badge_text,
        "verification_url": certificate.verification_url,
        "is_valid": certificate.is_valid(),
    }


@router.get("/certificate/verify/{certificate_id}", response_model=Dict)
async def verify_certificate(
    certificate_id: str,
    verification_hash: Optional[str] = Query(None),
) -> Dict:
    """
    Verify certificate authenticity and validity.
    
    Checks that certificate is unaltered and currently valid.
    """
    return _certificate_service.verify_certificate(
        certificate_id=certificate_id,
        verification_hash=verification_hash,
    )


# =============================================================================
# Verification Endpoints
# =============================================================================

@router.post("/verify", response_model=Dict)
async def verify_tk_attribution(request: VerificationRequest) -> Dict:
    """
    Verify TK attribution claim against traditional preparation.
    
    Compares chemical profiles and provenance to validate claims.
    """
    # Build chemical profile if provided
    chemical_profile = None
    if request.chemical_profile:
        chemical_profile = ChemicalProfile(
            cannabinoids=request.chemical_profile.get("cannabinoids", {}),
            terpenes=request.chemical_profile.get("terpenes", {}),
            flavonoids=request.chemical_profile.get("flavonoids", {}),
        )
    
    # Build provenance record if provided
    provenance = None
    if request.provenance:
        provenance = ProvenanceRecord(
            source_community=request.provenance.get("source_community", ""),
            source_region=request.provenance.get("source_region", ""),
            source_country=request.provenance.get("source_country", ""),
            cultivation_method=request.provenance.get("cultivation_method", ""),
            status=ProvenanceStatus.SELF_ATTESTED,
        )
    
    report = _verification_service.verify_traditional_preparation(
        claimed_formulation=request.claimed_formulation,
        traditional_recipe=request.traditional_recipe,
        chemical_profile=chemical_profile,
        provenance=provenance,
    )
    
    return {
        "report_id": report.report_id,
        "status": report.status.value,
        "overall_match_score": report.overall_match_score,
        "scores": {
            "chemical_match": report.chemical_match_score,
            "method_match": report.method_match_score,
            "provenance": report.provenance_score,
        },
        "claimed_formulation": report.claimed_formulation,
        "claimed_traditional_method": report.claimed_traditional_method,
        "claimed_community": report.claimed_community,
        "discrepancies": report.discrepancies,
        "recommendations": report.recommendations,
        "verification_hash": report.verification_hash,
    }


@router.get("/preparations", response_model=List[Dict])
async def list_traditional_preparations() -> List[Dict]:
    """
    List available traditional preparations for verification.
    
    Returns preparation keys, names, communities, and methods.
    """
    return _verification_service.list_traditional_preparations()


@router.post("/provenance/record", response_model=Dict)
async def record_provenance(request: ProvenanceRecordRequest) -> Dict:
    """
    Record provenance information for future verification.
    
    Creates blockchain-verified supply chain record.
    """
    coords = None
    if request.gps_coordinates and len(request.gps_coordinates) == 2:
        coords = (request.gps_coordinates[0], request.gps_coordinates[1])
    
    harvest_date = None
    if request.harvest_date:
        harvest_date = datetime.fromisoformat(request.harvest_date)
    
    provenance = ProvenanceRecord(
        source_community=request.source_community,
        source_region=request.source_region,
        source_country=request.source_country,
        gps_coordinates=coords,
        cultivation_method=request.cultivation_method,
        harvest_date=harvest_date,
        harvest_timing_notes=request.harvest_timing_notes,
        processing_location=request.processing_location,
        processing_method=request.processing_method,
        status=ProvenanceStatus.BLOCKCHAIN_VERIFIED,
    )
    
    return _verification_service.record_provenance(provenance)


# =============================================================================
# Health Check
# =============================================================================

@router.get("/health", response_model=Dict)
async def tkpath_health() -> Dict:
    """Health check for TKPath services."""
    return {
        "status": "healthy",
        "service": "TKPath",
        "components": {
            "attribution": "operational",
            "compensation": "operational",
            "certificate": "operational",
            "verification": "operational",
        },
        "pool_balance": _compensation_service.get_pool_status()["current_balance"],
        "known_communities": len(_attribution_service.list_known_communities()),
        "traditional_preparations": len(_verification_service.list_traditional_preparations()),
    }
