"""
TKPath - Traditional Knowledge Attribution & Compensation Engine

The first economically self-sustaining model for indigenous knowledge compensation.
TK features are FREE by default; companies opting OUT pay surcharges that fund
EquiPath compensation to indigenous communities.

Trade Secret: Attribution algorithms, community weighting, contribution scoring.
"""

from .attribution import (
    TKPathAttribution,
    TKContribution,
    CommunityAttribution,
    KnowledgeElement,
)
from .compensation import (
    TKCompensationEngine,
    CompensationTransaction,
    CompensationStatus,
)
from .certificate import (
    TKCertificateGenerator,
    AttributionCertificate,
    CertificateStatus,
)
from .verification import (
    TKVerificationLab,
    VerificationReport,
    ProvenanceRecord,
)

__all__ = [
    # Attribution
    "TKPathAttribution",
    "TKContribution",
    "CommunityAttribution",
    "KnowledgeElement",
    # Compensation
    "TKCompensationEngine",
    "CompensationTransaction",
    "CompensationStatus",
    # Certificates
    "TKCertificateGenerator",
    "AttributionCertificate",
    "CertificateStatus",
    # Verification
    "TKVerificationLab",
    "VerificationReport",
    "ProvenanceRecord",
]
