"""
EthnoPath API Router - Traditional Knowledge Digitization Endpoints

Exposes federated TK digitization, community validation, and DAO governance
capabilities through RESTful API.
"""

from fastapi import APIRouter, HTTPException, Depends, Body
from typing import List, Optional, Dict, Any
from pydantic import BaseModel, Field
from datetime import datetime

from backend.services.ethnopath import (
    EthnoPathDigitizer,
    CommunityValidator,
    EthnoPathGovernance,
)
from backend.services.ethnopath.digitizer import (
    KnowledgeDomain,
    KnowledgeType,
    SensitivityLevel,
    ConsentStatus,
    CommunityProfile,
    KnowledgeHolder,
    TraditionalKnowledgeEntry,
)
from backend.services.ethnopath.validator import (
    ValidationLevel,
    ValidatorRole,
    ValidationRound,
    ValidatorProfile,
)
from backend.services.ethnopath.governance import (
    GovernanceRole,
    ProposalType,
    VoteType,
    GovernanceProposal,
    GovernanceParticipant,
    AccessGrant,
)


router = APIRouter(prefix="/ethnopath", tags=["EthnoPath"])

# Global service instances
digitizer = EthnoPathDigitizer()
validator = CommunityValidator()
governance = EthnoPathGovernance()


# ==========================================================================
# Request/Response Models
# ==========================================================================

class CommunityRegistrationRequest(BaseModel):
    """Request to register community."""
    name: str
    location: str
    indigenous_group: str
    contact_email: str
    language: Optional[str] = None
    population_size: Optional[int] = None


class KnowledgeHolderRegistrationRequest(BaseModel):
    """Request to register knowledge holder."""
    community_id: str
    name: str
    role: str
    areas_of_expertise: List[str]
    years_of_experience: Optional[int] = None


class KnowledgeEntryRequest(BaseModel):
    """Request to create knowledge entry."""
    community_id: str
    holder_id: str
    domain: KnowledgeDomain
    knowledge_type: KnowledgeType
    sensitivity_level: SensitivityLevel
    title: str
    description: str
    traditional_name: Optional[str] = None
    tags: Optional[List[str]] = None


class ConsentGrantRequest(BaseModel):
    """Request to grant consent."""
    entry_id: str
    holder_id: str
    consent_type: str
    restrictions: Optional[List[str]] = None


class ValidationRoundRequest(BaseModel):
    """Request to initiate validation round."""
    entry_id: str
    community_id: str
    validation_level: ValidationLevel
    invited_validators: Optional[List[str]] = None


class ValidationVoteRequest(BaseModel):
    """Request to submit validation vote."""
    round_id: str
    validator_id: str
    approved: bool
    confidence_score: float = Field(ge=0.0, le=1.0)
    comments: Optional[str] = ""
    suggested_revisions: Optional[List[str]] = None
    cultural_concerns: Optional[List[str]] = None


class ValidatorRegistrationRequest(BaseModel):
    """Request to register validator."""
    community_id: str
    display_name: str
    role: ValidatorRole
    expertise_domains: Optional[List[str]] = None
    years_of_knowledge: Optional[int] = None


class ProposalCreationRequest(BaseModel):
    """Request to create governance proposal."""
    community_id: str
    created_by: str
    proposal_type: ProposalType
    title: str
    description: str
    requested_entry_ids: Optional[List[str]] = None
    requested_access_level: Optional[str] = None
    intended_use: Optional[str] = ""
    benefit_sharing_terms: Optional[str] = ""


class GovernanceVoteRequest(BaseModel):
    """Request to cast governance vote."""
    proposal_id: str
    participant_id: str
    vote_type: VoteType
    voting_power_used: float = 1.0
    rationale: Optional[str] = ""
    cultural_considerations: Optional[List[str]] = None


class ParticipantRegistrationRequest(BaseModel):
    """Request to register governance participant."""
    community_id: str
    display_name: str
    role: GovernanceRole


# ==========================================================================
# Community Management Endpoints
# ==========================================================================

@router.post("/communities/register", response_model=CommunityProfile)
async def register_community(request: CommunityRegistrationRequest):
    """
    Register a new community for traditional knowledge digitization.
    
    Returns community profile with federated node configuration.
    """
    try:
        community = digitizer.register_community(
            name=request.name,
            location=request.location,
            indigenous_group=request.indigenous_group,
            contact_email=request.contact_email,
            language=request.language,
            population_size=request.population_size,
        )
        return community
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))


@router.get("/communities/{community_id}", response_model=CommunityProfile)
async def get_community(community_id: str):
    """Get community profile by ID."""
    community = digitizer._communities.get(community_id)
    if not community:
        raise HTTPException(status_code=404, detail="Community not found")
    return community


@router.get("/communities", response_model=List[CommunityProfile])
async def list_communities():
    """List all registered communities."""
    return list(digitizer._communities.values())


@router.get("/communities/{community_id}/statistics")
async def get_community_statistics(community_id: str):
    """Get statistics for a community."""
    try:
        stats = digitizer.get_community_statistics(community_id)
        return stats
    except ValueError as e:
        raise HTTPException(status_code=404, detail=str(e))


# ==========================================================================
# Knowledge Holder Endpoints
# ==========================================================================

@router.post("/knowledge-holders/register", response_model=KnowledgeHolder)
async def register_knowledge_holder(request: KnowledgeHolderRegistrationRequest):
    """Register a knowledge holder within a community."""
    try:
        holder = digitizer.register_knowledge_holder(
            community_id=request.community_id,
            name=request.name,
            role=request.role,
            areas_of_expertise=request.areas_of_expertise,
            years_of_experience=request.years_of_experience,
        )
        return holder
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))


@router.get("/knowledge-holders/{holder_id}", response_model=KnowledgeHolder)
async def get_knowledge_holder(holder_id: str):
    """Get knowledge holder by ID."""
    holder = digitizer._knowledge_holders.get(holder_id)
    if not holder:
        raise HTTPException(status_code=404, detail="Knowledge holder not found")
    return holder


@router.get("/communities/{community_id}/knowledge-holders", response_model=List[KnowledgeHolder])
async def list_community_knowledge_holders(community_id: str):
    """List all knowledge holders in a community."""
    holders = [
        h for h in digitizer._knowledge_holders.values()
        if h.community_id == community_id
    ]
    return holders


# ==========================================================================
# Knowledge Entry Endpoints
# ==========================================================================

@router.post("/knowledge-entries/create", response_model=TraditionalKnowledgeEntry)
async def create_knowledge_entry(request: KnowledgeEntryRequest):
    """Create a new traditional knowledge entry."""
    try:
        entry = digitizer.create_knowledge_entry(
            community_id=request.community_id,
            holder_id=request.holder_id,
            domain=request.domain,
            knowledge_type=request.knowledge_type,
            sensitivity_level=request.sensitivity_level,
            title=request.title,
            description=request.description,
            traditional_name=request.traditional_name,
            tags=request.tags,
        )
        return entry
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))


@router.get("/knowledge-entries/{entry_id}", response_model=TraditionalKnowledgeEntry)
async def get_knowledge_entry(entry_id: str):
    """Get knowledge entry by ID."""
    entry = digitizer._knowledge_entries.get(entry_id)
    if not entry:
        raise HTTPException(status_code=404, detail="Knowledge entry not found")
    return entry


@router.get("/communities/{community_id}/knowledge-entries", response_model=List[TraditionalKnowledgeEntry])
async def list_community_knowledge_entries(community_id: str):
    """List all knowledge entries for a community."""
    entries = [
        e for e in digitizer._knowledge_entries.values()
        if e.community_id == community_id
    ]
    return entries


@router.post("/knowledge-entries/{entry_id}/consent")
async def grant_consent(entry_id: str, request: ConsentGrantRequest):
    """Grant consent for knowledge entry usage."""
    try:
        digitizer.grant_consent(
            entry_id=request.entry_id,
            holder_id=request.holder_id,
            consent_type=request.consent_type,
            restrictions=request.restrictions,
        )
        return {"status": "success", "message": "Consent granted"}
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))


@router.get("/knowledge-entries/{entry_id}/access-check")
async def check_entry_access(
    entry_id: str,
    community_id: str,
    sensitivity_level: SensitivityLevel,
):
    """Check if community can access knowledge entry."""
    can_access = digitizer.can_access_entry(entry_id, community_id, sensitivity_level)
    return {"entry_id": entry_id, "can_access": can_access}


# ==========================================================================
# Federated Aggregation Endpoints
# ==========================================================================

@router.post("/federated/aggregate-patterns")
async def aggregate_knowledge_patterns(
    domain: KnowledgeDomain,
    minimum_communities: int = 3,
):
    """Aggregate knowledge patterns across communities using federated learning."""
    try:
        patterns = digitizer.aggregate_knowledge_patterns(
            domain=domain,
            minimum_communities=minimum_communities,
        )
        return {"patterns": patterns}
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))


@router.get("/digitization/metrics")
async def get_digitization_metrics():
    """Get digitization system metrics."""
    metrics = digitizer.get_digitization_metrics()
    return metrics


# ==========================================================================
# Validator Management Endpoints
# ==========================================================================

@router.post("/validators/register", response_model=ValidatorProfile)
async def register_validator(request: ValidatorRegistrationRequest):
    """Register a community validator."""
    try:
        validator_profile = validator.register_validator(
            community_id=request.community_id,
            display_name=request.display_name,
            role=request.role,
            expertise_domains=request.expertise_domains,
            years_of_knowledge=request.years_of_knowledge,
        )
        return validator_profile
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))


@router.get("/validators/{validator_id}", response_model=ValidatorProfile)
async def get_validator(validator_id: str):
    """Get validator profile by ID."""
    validator_profile = validator._validators.get(validator_id)
    if not validator_profile:
        raise HTTPException(status_code=404, detail="Validator not found")
    return validator_profile


# ==========================================================================
# Validation Round Endpoints
# ==========================================================================

@router.post("/validation/rounds/initiate", response_model=ValidationRound)
async def initiate_validation_round(request: ValidationRoundRequest):
    """Initiate a validation round for knowledge entry."""
    try:
        round_obj = validator.initiate_validation_round(
            entry_id=request.entry_id,
            community_id=request.community_id,
            validation_level=request.validation_level,
            invited_validators=request.invited_validators,
        )
        return round_obj
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))


@router.post("/validation/rounds/{round_id}/vote")
async def submit_validation_vote(round_id: str, request: ValidationVoteRequest):
    """Submit validation vote for a round."""
    try:
        vote = validator.submit_validation_vote(
            round_id=request.round_id,
            validator_id=request.validator_id,
            approved=request.approved,
            confidence_score=request.confidence_score,
            comments=request.comments,
            suggested_revisions=request.suggested_revisions,
            cultural_concerns=request.cultural_concerns,
        )
        return vote
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))


@router.get("/validation/rounds/{round_id}", response_model=ValidationRound)
async def get_validation_round(round_id: str):
    """Get validation round details."""
    round_obj = validator._rounds.get(round_id)
    if not round_obj:
        raise HTTPException(status_code=404, detail="Validation round not found")
    return round_obj


@router.get("/validation/metrics")
async def get_validation_metrics():
    """Get validation system metrics."""
    metrics = validator.get_validation_metrics()
    return metrics


# ==========================================================================
# Multi-Community Consensus Endpoints
# ==========================================================================

@router.post("/validation/multi-community/initiate")
async def initiate_multi_community_consensus(
    entry_id: str,
    communities_invited: List[str] = Body(...),
):
    """Initiate multi-community consensus validation."""
    try:
        consensus = validator.initiate_multi_community_consensus(
            entry_id=entry_id,
            communities_invited=communities_invited,
        )
        return consensus
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))


@router.get("/validation/multi-community/{consensus_id}")
async def get_multi_community_consensus(consensus_id: str):
    """Get multi-community consensus status."""
    consensus = validator._consensus_records.get(consensus_id)
    if not consensus:
        raise HTTPException(status_code=404, detail="Consensus record not found")
    
    # Update consensus before returning
    updated_consensus = validator.update_multi_community_consensus(consensus_id)
    return updated_consensus


# ==========================================================================
# Governance Participant Endpoints
# ==========================================================================

@router.post("/governance/participants/register", response_model=GovernanceParticipant)
async def register_governance_participant(request: ParticipantRegistrationRequest):
    """Register governance participant."""
    try:
        participant = governance.register_participant(
            community_id=request.community_id,
            display_name=request.display_name,
            role=request.role,
        )
        return participant
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))


@router.get("/governance/participants/{participant_id}", response_model=GovernanceParticipant)
async def get_governance_participant(participant_id: str):
    """Get governance participant by ID."""
    participant = governance._participants.get(participant_id)
    if not participant:
        raise HTTPException(status_code=404, detail="Participant not found")
    return participant


# ==========================================================================
# Governance Proposal Endpoints
# ==========================================================================

@router.post("/governance/proposals/create", response_model=GovernanceProposal)
async def create_governance_proposal(request: ProposalCreationRequest):
    """Create governance proposal."""
    try:
        proposal = governance.create_proposal(
            community_id=request.community_id,
            created_by=request.created_by,
            proposal_type=request.proposal_type,
            title=request.title,
            description=request.description,
            requested_entry_ids=request.requested_entry_ids,
            requested_access_level=request.requested_access_level,
            intended_use=request.intended_use,
            benefit_sharing_terms=request.benefit_sharing_terms,
        )
        return proposal
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))


@router.post("/governance/proposals/{proposal_id}/submit")
async def submit_governance_proposal(proposal_id: str):
    """Submit proposal for voting."""
    try:
        proposal = governance.submit_proposal(proposal_id)
        return proposal
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))


@router.post("/governance/proposals/{proposal_id}/vote")
async def cast_governance_vote(proposal_id: str, request: GovernanceVoteRequest):
    """Cast vote on governance proposal."""
    try:
        vote = governance.cast_vote(
            proposal_id=request.proposal_id,
            participant_id=request.participant_id,
            vote_type=request.vote_type,
            voting_power_used=request.voting_power_used,
            rationale=request.rationale,
            cultural_considerations=request.cultural_considerations,
        )
        return vote
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))


@router.get("/governance/proposals/{proposal_id}", response_model=GovernanceProposal)
async def get_governance_proposal(proposal_id: str):
    """Get governance proposal by ID."""
    proposal = governance._proposals.get(proposal_id)
    if not proposal:
        raise HTTPException(status_code=404, detail="Proposal not found")
    return proposal


@router.get("/governance/proposals", response_model=List[GovernanceProposal])
async def list_governance_proposals(
    community_id: Optional[str] = None,
    status: Optional[str] = None,
):
    """List governance proposals with optional filters."""
    proposals = list(governance._proposals.values())
    
    if community_id:
        proposals = [p for p in proposals if p.community_id == community_id]
    
    if status:
        proposals = [p for p in proposals if p.status == status]
    
    return proposals


# ==========================================================================
# Access Grant Endpoints
# ==========================================================================

@router.get("/governance/grants/{grant_id}", response_model=AccessGrant)
async def get_access_grant(grant_id: str):
    """Get access grant by ID."""
    grant = governance._grants.get(grant_id)
    if not grant:
        raise HTTPException(status_code=404, detail="Access grant not found")
    return grant


@router.post("/governance/grants/{grant_id}/revoke")
async def revoke_access_grant(
    grant_id: str,
    revocation_reason: str = Body(...),
):
    """Revoke access grant."""
    try:
        grant = governance.revoke_grant(grant_id, revocation_reason)
        return grant
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))


@router.get("/governance/grants/check-access")
async def check_access_permission(
    participant_id: str,
    entry_id: str,
):
    """Check if participant has access permission for entry."""
    has_access = governance.check_access_permission(participant_id, entry_id)
    return {"participant_id": participant_id, "entry_id": entry_id, "has_access": has_access}


@router.get("/governance/metrics")
async def get_governance_metrics():
    """Get governance system metrics."""
    metrics = governance.get_governance_metrics()
    return metrics


# ==========================================================================
# Health Check
# ==========================================================================

@router.get("/health")
async def health_check():
    """EthnoPath system health check."""
    return {
        "status": "healthy",
        "components": {
            "digitizer": "operational",
            "validator": "operational",
            "governance": "operational",
        },
        "timestamp": datetime.utcnow().isoformat(),
    }
