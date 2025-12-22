"""
EthnoPath - Traditional Knowledge Digitization and Preservation Engine

Trade Secret ID: TS-ETH-001
Classification: TOP SECRET - Proprietary Trade Secret
Estimated Value: $4.6 Billion
Competitive Advantage Duration: 12-14 years

Core Technology Components:
- Federated Learning Infrastructure: Train models without centralizing sensitive TK
- DAO Governance Framework: Elder veto authority, quadratic voting
- Community Validation Networks: Multi-stakeholder consensus for TK entries
- EquiPath Integration: Perpetual rights enforcement, automated compensation

Ecosystem Position: Traditional Knowledge Foundation Layer

Downstream Pathway Dependencies:
- ChemPath: Traditional preparation methods, compound combinations
- BioPath: Therapeutic efficacy signals from traditional use
- ToxPath: Historical safety data from community knowledge
- ClinPath: Culturally informed trial design inputs
- RegPath: Traditional use documentation for regulatory submissions
- EquiPath: Attribution data for compensation calculations
- RetroPath: Living knowledge validation for historical recovery

Cross-Sector Applications:
- Therapeutic development
- Climate adaptation
- Sustainable agriculture
- Cultural heritage preservation
"""

from backend.services.ethnopath.digitizer import (
    EthnoPathDigitizer,
    KnowledgeDomain,
    KnowledgeType,
    SensitivityLevel,
    ConsentStatus,
    VerificationStatus,
    ConsentRecord,
    FederatedNode,
    CommunityProfile,
    KnowledgeHolder,
    TraditionalKnowledgeEntry,
    DigitizationSession,
    FederatedLearningConfig,
)

from backend.services.ethnopath.validator import (
    CommunityValidator,
    ValidationLevel,
    ValidatorRole,
    ConsensusType,
    ValidationStatus,
    ValidatorProfile,
    ValidationVote,
    ValidationRound,
    CommunityConsensus,
    ValidationMetrics,
)

from backend.services.ethnopath.governance import (
    EthnoPathGovernance,
    GovernanceRole,
    ProposalType,
    ProposalStatus,
    VetoAuthority,
    VoteType,
    TKComplianceLevel,
    PricingTier,
    GovernanceParticipant,
    GovernanceVote,
    GovernanceProposal,
    AccessGrant,
    GovernanceMetrics,
)

__all__ = [
    # Digitizer
    "EthnoPathDigitizer",
    "KnowledgeDomain",
    "KnowledgeType",
    "SensitivityLevel",
    "ConsentStatus",
    "VerificationStatus",
    "ConsentRecord",
    "FederatedNode",
    "CommunityProfile",
    "KnowledgeHolder",
    "TraditionalKnowledgeEntry",
    "DigitizationSession",
    "FederatedLearningConfig",
    # Validator
    "CommunityValidator",
    "ValidationLevel",
    "ValidatorRole",
    "ConsensusType",
    "ValidationStatus",
    "ValidatorProfile",
    "ValidationVote",
    "ValidationRound",
    "CommunityConsensus",
    "ValidationMetrics",
    # Governance
    "EthnoPathGovernance",
    "GovernanceRole",
    "ProposalType",
    "ProposalStatus",
    "VetoAuthority",
    "VoteType",
    "TKComplianceLevel",
    "PricingTier",
    "GovernanceParticipant",
    "GovernanceVote",
    "GovernanceProposal",
    "AccessGrant",
    "GovernanceMetrics",
]
