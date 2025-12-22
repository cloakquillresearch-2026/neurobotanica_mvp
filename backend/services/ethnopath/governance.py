"""
EthnoPath Governance - DAO Governance with Elder Veto Authority

Trade Secret: Elder veto mechanisms, quadratic voting algorithms, access request
approval workflows, cultural protocol enforcement, community governance models.

Core Capabilities:
- DAO governance with elder council oversight
- Access request and approval workflows
- Cultural protocol enforcement
- Community governance and decision-making
"""

from enum import Enum
from typing import Dict, List, Optional, Set, Any
from pydantic import BaseModel, Field
from datetime import datetime, timedelta
import uuid
import math


class GovernanceRole(str, Enum):
    """Governance participant role."""
    COMMUNITY_MEMBER = "community_member"
    ELDER = "elder"
    KNOWLEDGE_HOLDER = "knowledge_holder"
    ADMINISTRATOR = "administrator"
    RESEARCHER = "researcher"
    EXTERNAL_PARTY = "external_party"


class ProposalType(str, Enum):
    """Governance proposal types."""
    ACCESS_REQUEST = "access_request"           # Request to access TK
    KNOWLEDGE_UPDATE = "knowledge_update"       # Update existing TK
    GOVERNANCE_CHANGE = "governance_change"     # Change governance rules
    BENEFIT_SHARING = "benefit_sharing"         # Benefit sharing agreement
    CULTURAL_PROTOCOL = "cultural_protocol"     # Cultural protocol decision
    EMERGENCY_RESPONSE = "emergency_response"   # Emergency decision


class VoteType(str, Enum):
    """Vote type."""
    APPROVE = "approve"
    REJECT = "reject"
    ABSTAIN = "abstain"
    VETO = "veto"  # Elder veto


class ProposalStatus(str, Enum):
    """Proposal status."""
    DRAFT = "draft"
    SUBMITTED = "submitted"
    IN_VOTING = "in_voting"
    APPROVED = "approved"
    REJECTED = "rejected"
    VETOED = "vetoed"
    IMPLEMENTED = "implemented"
    WITHDRAWN = "withdrawn"


class VetoAuthority(str, Enum):
    """Veto authority level."""
    NONE = "none"
    ELDER_COUNCIL = "elder_council"
    CEREMONIAL_LEADER = "ceremonial_leader"
    UNANIMOUS_ELDERS = "unanimous_elders"


class TKComplianceLevel(str, Enum):
    """Traditional Knowledge compliance level."""
    FULL_COMPLIANCE = "full_compliance"           # Default: Full TK attribution + compensation
    NO_ATTRIBUTION = "no_attribution"             # Opt-out: No source attribution
    NO_COMPENSATION = "no_compensation"           # Opt-out: No TK holder compensation  
    NO_TK_FULL = "no_tk_full"                     # Opt-out: Complete TK opt-out


class PricingTier(str, Enum):
    """Pricing tier for access grants."""
    BASIC = "basic"           # $99/mo - Community members only
    STANDARD = "standard"     # $299/mo - Commercial use
    PREMIUM = "premium"       # $799/mo - Full IP rights
    ENTERPRISE = "enterprise" # Custom - Large organizations


class GovernanceParticipant(BaseModel):
    """Governance participant profile."""
    participant_id: str = Field(default_factory=lambda: str(uuid.uuid4()))
    community_id: str
    
    # Identity
    display_name: str
    role: GovernanceRole
    
    # Voting power (for quadratic voting)
    base_voting_power: float = 1.0
    reputation_score: float = 1.0  # 0-1 scale
    
    # Participation
    proposals_created: int = 0
    votes_cast: int = 0
    participation_rate: float = 0.0
    
    # Status
    is_active: bool = True
    joined_governance: datetime = Field(default_factory=datetime.utcnow)


class GovernanceVote(BaseModel):
    """Individual governance vote."""
    vote_id: str = Field(default_factory=lambda: str(uuid.uuid4()))
    proposal_id: str
    participant_id: str
    
    # Vote
    vote_type: VoteType
    voting_power_used: float = 1.0  # Can use partial power
    
    # Rationale
    rationale: str = ""
    cultural_considerations: List[str] = Field(default_factory=list)
    
    # Context
    voted_at: datetime = Field(default_factory=datetime.utcnow)
    voter_role: GovernanceRole = GovernanceRole.COMMUNITY_MEMBER


class GovernanceProposal(BaseModel):
    """Governance proposal."""
    proposal_id: str = Field(default_factory=lambda: str(uuid.uuid4()))
    community_id: str
    
    # Proposal details
    proposal_type: ProposalType
    title: str
    description: str
    
    # Requester
    created_by: str  # participant_id
    
    # Access request specifics (if applicable)
    requested_entry_ids: List[str] = Field(default_factory=list)
    requested_access_level: Optional[str] = None
    intended_use: str = ""
    benefit_sharing_terms: str = ""
    
    # Voting
    votes: List[GovernanceVote] = Field(default_factory=list)
    voting_starts: datetime = Field(default_factory=datetime.utcnow)
    voting_ends: datetime = Field(default_factory=lambda: datetime.utcnow() + timedelta(days=7))
    
    # Requirements
    minimum_participation: int = 3
    approval_threshold: float = 0.66  # Supermajority
    veto_authority: VetoAuthority = VetoAuthority.ELDER_COUNCIL
    
    # Status
    status: ProposalStatus = ProposalStatus.DRAFT
    final_decision: Optional[bool] = None
    decision_rationale: str = ""
    
    # Timestamps
    submitted_at: Optional[datetime] = None
    decided_at: Optional[datetime] = None


class AccessGrant(BaseModel):
    """Access grant record."""
    grant_id: str = Field(default_factory=lambda: str(uuid.uuid4()))
    proposal_id: str
    
    # Grant details
    grantee_id: str  # participant_id
    granted_entry_ids: List[str]
    access_level: str
    
    # Terms
    intended_use: str
    benefit_sharing_terms: str
    cultural_protocols: List[str] = Field(default_factory=list)
    
    # TRADE SECRET: Ethical Pricing & TK Compliance
    pricing_tier: PricingTier = PricingTier.STANDARD
    tk_compliance_level: TKComplianceLevel = TKComplianceLevel.FULL_COMPLIANCE
    base_price_usd: float = 0.0
    surcharge_amount_usd: float = 0.0
    total_price_usd: float = 0.0
    compensation_pool_contribution_usd: float = 0.0
    
    # Validity
    granted_at: datetime = Field(default_factory=datetime.utcnow)
    expires_at: Optional[datetime] = None
    is_revoked: bool = False
    revoked_at: Optional[datetime] = None
    revocation_reason: str = ""


class GovernanceMetrics(BaseModel):
    """Governance system metrics."""
    timestamp: datetime = Field(default_factory=datetime.utcnow)
    
    # Participation
    total_participants: int = 0
    active_participants: int = 0
    average_participation_rate: float = 0.0
    
    # Proposals
    total_proposals: int = 0
    active_proposals: int = 0
    approved_proposals: int = 0
    rejected_proposals: int = 0
    vetoed_proposals: int = 0
    
    # Access grants
    total_grants: int = 0
    active_grants: int = 0
    revoked_grants: int = 0
    
    # TRADE SECRET: Revenue & Compensation Tracking
    total_revenue_usd: float = 0.0
    total_tk_surcharges_usd: float = 0.0
    compensation_pool_usd: float = 0.0
    full_compliance_grants: int = 0
    opt_out_grants: int = 0
    
    # Efficiency
    average_decision_time_days: float = 0.0
    approval_rate: float = 0.0


class EthnoPathGovernance:
    """
    DAO Governance with Elder Veto Authority.
    
    Trade Secret: Elder veto mechanisms, quadratic voting, access approval
    workflows, cultural protocol enforcement, community governance models.
    
    Key Features:
    - Decentralized autonomous governance with elder oversight
    - Quadratic voting for equitable decision-making
    - Access request and approval workflows
    - Cultural protocol enforcement
    - Benefit-sharing agreement management
    """
    
    # ==========================================================================
    # TRADE SECRET: Proposal Type Requirements Matrix
    # ==========================================================================
    _PROPOSAL_REQUIREMENTS = {
        ProposalType.ACCESS_REQUEST: {
            "minimum_participation": 3,
            "approval_threshold": 0.66,
            "veto_authority": VetoAuthority.ELDER_COUNCIL,
            "voting_period_days": 7,
        },
        ProposalType.KNOWLEDGE_UPDATE: {
            "minimum_participation": 2,
            "approval_threshold": 0.60,
            "veto_authority": VetoAuthority.ELDER_COUNCIL,
            "voting_period_days": 5,
        },
        ProposalType.GOVERNANCE_CHANGE: {
            "minimum_participation": 5,
            "approval_threshold": 0.75,
            "veto_authority": VetoAuthority.UNANIMOUS_ELDERS,
            "voting_period_days": 14,
        },
        ProposalType.BENEFIT_SHARING: {
            "minimum_participation": 4,
            "approval_threshold": 0.70,
            "veto_authority": VetoAuthority.ELDER_COUNCIL,
            "voting_period_days": 10,
        },
        ProposalType.CULTURAL_PROTOCOL: {
            "minimum_participation": 3,
            "approval_threshold": 0.80,
            "veto_authority": VetoAuthority.CEREMONIAL_LEADER,
            "voting_period_days": 7,
        },
        ProposalType.EMERGENCY_RESPONSE: {
            "minimum_participation": 2,
            "approval_threshold": 0.66,
            "veto_authority": VetoAuthority.ELDER_COUNCIL,
            "voting_period_days": 1,
        },
    }
    
    # ==========================================================================
    # TRADE SECRET: Role Voting Power Multipliers
    # ==========================================================================
    _ROLE_VOTING_POWER = {
        GovernanceRole.COMMUNITY_MEMBER: 1.0,
        GovernanceRole.ELDER: 3.0,           # Elders get 3x power
        GovernanceRole.KNOWLEDGE_HOLDER: 2.0,
        GovernanceRole.ADMINISTRATOR: 1.5,
        GovernanceRole.RESEARCHER: 0.5,      # External researchers limited
        GovernanceRole.EXTERNAL_PARTY: 0.0,  # No voting power
    }
    
    # ==========================================================================
    # TRADE SECRET: Ethical Pricing - TK Opt-Out Surcharges
    # ==========================================================================
    _BASE_PRICING = {
        PricingTier.BASIC: 99.0,
        PricingTier.STANDARD: 299.0,
        PricingTier.PREMIUM: 799.0,
        PricingTier.ENTERPRISE: 2499.0,
    }
    
    _TK_SURCHARGES = {
        TKComplianceLevel.FULL_COMPLIANCE: {"percentage": 0.0, "flat_fee": 0.0},
        TKComplianceLevel.NO_ATTRIBUTION: {"percentage": 0.15, "flat_fee": 150.0},
        TKComplianceLevel.NO_COMPENSATION: {"percentage": 0.25, "flat_fee": 250.0},
        TKComplianceLevel.NO_TK_FULL: {"percentage": 0.40, "flat_fee": 500.0},
    }
    
    def __init__(self):
        self._participants: Dict[str, GovernanceParticipant] = {}
        self._proposals: Dict[str, GovernanceProposal] = {}
        self._grants: Dict[str, AccessGrant] = {}
        self._metrics = GovernanceMetrics()
    
    # ==========================================================================
    # Participant Management
    # ==========================================================================
    
    def register_participant(
        self,
        community_id: str,
        display_name: str,
        role: GovernanceRole,
    ) -> GovernanceParticipant:
        """Register governance participant."""
        base_power = self._ROLE_VOTING_POWER.get(role, 1.0)
        
        participant = GovernanceParticipant(
            community_id=community_id,
            display_name=display_name,
            role=role,
            base_voting_power=base_power,
        )
        
        self._participants[participant.participant_id] = participant
        self._metrics.total_participants += 1
        self._metrics.active_participants += 1
        
        return participant
    
    # ==========================================================================
    # TRADE SECRET: Proposal Creation and Submission
    # ==========================================================================
    
    def create_proposal(
        self,
        community_id: str,
        created_by: str,
        proposal_type: ProposalType,
        title: str,
        description: str,
        requested_entry_ids: Optional[List[str]] = None,
        requested_access_level: Optional[str] = None,
        intended_use: str = "",
        benefit_sharing_terms: str = "",
    ) -> GovernanceProposal:
        """
        Create governance proposal.
        Trade Secret: Automatic requirement configuration based on type.
        """
        # Get requirements for this proposal type
        requirements = self._PROPOSAL_REQUIREMENTS.get(
            proposal_type,
            {
                "minimum_participation": 3,
                "approval_threshold": 0.66,
                "veto_authority": VetoAuthority.ELDER_COUNCIL,
                "voting_period_days": 7,
            }
        )
        
        voting_ends = datetime.utcnow() + timedelta(days=requirements["voting_period_days"])
        
        proposal = GovernanceProposal(
            community_id=community_id,
            created_by=created_by,
            proposal_type=proposal_type,
            title=title,
            description=description,
            requested_entry_ids=requested_entry_ids or [],
            requested_access_level=requested_access_level,
            intended_use=intended_use,
            benefit_sharing_terms=benefit_sharing_terms,
            voting_ends=voting_ends,
            minimum_participation=requirements["minimum_participation"],
            approval_threshold=requirements["approval_threshold"],
            veto_authority=requirements["veto_authority"],
        )
        
        self._proposals[proposal.proposal_id] = proposal
        
        # Update creator stats
        creator = self._participants.get(created_by)
        if creator:
            creator.proposals_created += 1
        
        return proposal
    
    def submit_proposal(self, proposal_id: str) -> GovernanceProposal:
        """Submit proposal for voting."""
        proposal = self._proposals.get(proposal_id)
        if not proposal:
            raise ValueError(f"Proposal not found: {proposal_id}")
        
        if proposal.status != ProposalStatus.DRAFT:
            raise ValueError(f"Proposal not in draft state: {proposal.status}")
        
        proposal.status = ProposalStatus.IN_VOTING
        proposal.submitted_at = datetime.utcnow()
        
        self._metrics.total_proposals += 1
        self._metrics.active_proposals += 1
        
        return proposal
    
    # ==========================================================================
    # TRADE SECRET: Quadratic Voting Algorithm
    # ==========================================================================
    
    def cast_vote(
        self,
        proposal_id: str,
        participant_id: str,
        vote_type: VoteType,
        voting_power_used: float = 1.0,
        rationale: str = "",
        cultural_considerations: Optional[List[str]] = None,
    ) -> GovernanceVote:
        """
        Cast vote on proposal.
        Trade Secret: Quadratic voting power calculation.
        """
        proposal = self._proposals.get(proposal_id)
        if not proposal:
            raise ValueError(f"Proposal not found: {proposal_id}")
        
        if proposal.status != ProposalStatus.IN_VOTING:
            raise ValueError(f"Proposal not accepting votes: {proposal.status}")
        
        if datetime.utcnow() > proposal.voting_ends:
            raise ValueError("Voting period has ended")
        
        participant = self._participants.get(participant_id)
        if not participant:
            raise ValueError(f"Participant not found: {participant_id}")
        
        # Calculate quadratic voting power
        # Power = sqrt(tokens_used) * base_power * reputation
        quadratic_power = (
            math.sqrt(voting_power_used)
            * participant.base_voting_power
            * participant.reputation_score
        )
        
        vote = GovernanceVote(
            proposal_id=proposal_id,
            participant_id=participant_id,
            vote_type=vote_type,
            voting_power_used=quadratic_power,
            rationale=rationale,
            cultural_considerations=cultural_considerations or [],
            voter_role=participant.role,
        )
        
        proposal.votes.append(vote)
        
        # Update participant stats
        participant.votes_cast += 1
        
        # Check if decision reached
        self._check_proposal_decision(proposal)
        
        return vote
    
    # ==========================================================================
    # TRADE SECRET: Elder Veto Mechanism
    # ==========================================================================
    
    def _check_proposal_decision(self, proposal: GovernanceProposal) -> None:
        """
        Check if proposal has reached decision.
        Trade Secret: Elder veto logic, quadratic vote aggregation.
        """
        # Check for elder veto
        if self._check_elder_veto(proposal):
            proposal.status = ProposalStatus.VETOED
            proposal.final_decision = False
            proposal.decision_rationale = "Elder council exercised veto authority"
            proposal.decided_at = datetime.utcnow()
            self._metrics.vetoed_proposals += 1
            self._metrics.active_proposals -= 1
            return
        
        # TRADE SECRET: If proposal has veto authority, wait for at least one elder vote
        if proposal.veto_authority != VetoAuthority.NONE:
            elder_votes = [v for v in proposal.votes if v.voter_role == GovernanceRole.ELDER]
            if not elder_votes:
                return  # Wait for elder participation before deciding
        
        # Check if voting period ended or minimum participation met
        voting_ended = datetime.utcnow() > proposal.voting_ends
        min_participation_met = len(proposal.votes) >= proposal.minimum_participation
        
        if not (voting_ended or min_participation_met):
            return  # Still collecting votes
        
        # Calculate approval percentage
        total_power = sum(v.voting_power_used for v in proposal.votes if v.vote_type != VoteType.ABSTAIN)
        approved_power = sum(
            v.voting_power_used for v in proposal.votes
            if v.vote_type == VoteType.APPROVE
        )
        
        approval_percentage = approved_power / total_power if total_power > 0 else 0.0
        
        # Check threshold
        if approval_percentage >= proposal.approval_threshold:
            proposal.status = ProposalStatus.APPROVED
            proposal.final_decision = True
            proposal.decision_rationale = f"Approved with {approval_percentage:.1%} support"
            self._metrics.approved_proposals += 1
            
            # Create access grant if applicable
            if proposal.proposal_type == ProposalType.ACCESS_REQUEST:
                self._create_access_grant(proposal)
        else:
            proposal.status = ProposalStatus.REJECTED
            proposal.final_decision = False
            proposal.decision_rationale = f"Rejected: {approval_percentage:.1%} support (needed {proposal.approval_threshold:.1%})"
            self._metrics.rejected_proposals += 1
        
        proposal.decided_at = datetime.utcnow()
        self._metrics.active_proposals -= 1
    
    def _check_elder_veto(self, proposal: GovernanceProposal) -> bool:
        """
        Check if elders have exercised veto.
        Trade Secret: Elder veto authority logic.
        """
        if proposal.veto_authority == VetoAuthority.NONE:
            return False
        
        # Find elder votes
        elder_votes = [
            v for v in proposal.votes
            if v.voter_role == GovernanceRole.ELDER
        ]
        
        if not elder_votes:
            return False  # No elder votes yet
        
        # Check veto authority type
        if proposal.veto_authority == VetoAuthority.ELDER_COUNCIL:
            # Any elder can veto
            return any(v.vote_type == VoteType.VETO for v in elder_votes)
        
        elif proposal.veto_authority == VetoAuthority.UNANIMOUS_ELDERS:
            # All elders must veto
            veto_count = len([v for v in elder_votes if v.vote_type == VoteType.VETO])
            return veto_count == len(elder_votes) and veto_count > 0
        
        elif proposal.veto_authority == VetoAuthority.CEREMONIAL_LEADER:
            # Check for ceremonial leader veto (would need role tracking)
            # Simplified: Any elder veto for now
            return any(v.vote_type == VoteType.VETO for v in elder_votes)
        
        return False
    
    # ==========================================================================
    # Access Grant Management
    # ==========================================================================
    
    def calculate_grant_pricing(
        self,
        pricing_tier: PricingTier,
        tk_compliance_level: TKComplianceLevel,
    ) -> Dict[str, float]:
        """
        Calculate pricing for access grant with TK opt-out surcharges.
        Trade Secret: Ethical pricing algorithm.
        
        Returns:
            Dict with base_price, surcharge, total, and compensation_pool
        """
        base_price = self._BASE_PRICING[pricing_tier]
        surcharge_config = self._TK_SURCHARGES[tk_compliance_level]
        
        # Calculate surcharge
        percentage_surcharge = base_price * surcharge_config["percentage"]
        flat_fee = surcharge_config["flat_fee"]
        total_surcharge = percentage_surcharge + flat_fee
        
        # Total price
        total_price = base_price + total_surcharge
        
        # Compensation pool: All surcharges go to TK holders
        compensation_contribution = total_surcharge
        
        return {
            "base_price_usd": base_price,
            "surcharge_amount_usd": total_surcharge,
            "total_price_usd": total_price,
            "compensation_pool_contribution_usd": compensation_contribution,
        }
    
    def _create_access_grant(
        self,
        proposal: GovernanceProposal,
        pricing_tier: PricingTier = PricingTier.STANDARD,
        tk_compliance_level: TKComplianceLevel = TKComplianceLevel.FULL_COMPLIANCE,
    ) -> AccessGrant:
        """Create access grant from approved proposal with ethical pricing."""
        # Calculate pricing
        pricing = self.calculate_grant_pricing(pricing_tier, tk_compliance_level)
        
        grant = AccessGrant(
            proposal_id=proposal.proposal_id,
            grantee_id=proposal.created_by,
            granted_entry_ids=proposal.requested_entry_ids,
            access_level=proposal.requested_access_level or "read",
            intended_use=proposal.intended_use,
            benefit_sharing_terms=proposal.benefit_sharing_terms,
            expires_at=datetime.utcnow() + timedelta(days=365),  # 1 year default
            # Pricing fields
            pricing_tier=pricing_tier,
            tk_compliance_level=tk_compliance_level,
            base_price_usd=pricing["base_price_usd"],
            surcharge_amount_usd=pricing["surcharge_amount_usd"],
            total_price_usd=pricing["total_price_usd"],
            compensation_pool_contribution_usd=pricing["compensation_pool_contribution_usd"],
        )
        
        self._grants[grant.grant_id] = grant
        
        # Update metrics
        self._metrics.total_grants += 1
        self._metrics.active_grants += 1
        self._metrics.total_revenue_usd += grant.total_price_usd
        self._metrics.total_tk_surcharges_usd += grant.surcharge_amount_usd
        self._metrics.compensation_pool_usd += grant.compensation_pool_contribution_usd
        
        if tk_compliance_level == TKComplianceLevel.FULL_COMPLIANCE:
            self._metrics.full_compliance_grants += 1
        else:
            self._metrics.opt_out_grants += 1
        
        return grant
    
    def revoke_grant(
        self,
        grant_id: str,
        revocation_reason: str,
    ) -> AccessGrant:
        """Revoke access grant."""
        grant = self._grants.get(grant_id)
        if not grant:
            raise ValueError(f"Grant not found: {grant_id}")
        
        if grant.is_revoked:
            raise ValueError("Grant already revoked")
        
        grant.is_revoked = True
        grant.revoked_at = datetime.utcnow()
        grant.revocation_reason = revocation_reason
        
        self._metrics.active_grants -= 1
        self._metrics.revoked_grants += 1
        
        return grant
    
    def check_access_permission(
        self,
        participant_id: str,
        entry_id: str,
    ) -> bool:
        """Check if participant has access to knowledge entry."""
        # Find active grants for this participant
        active_grants = [
            g for g in self._grants.values()
            if g.grantee_id == participant_id
            and not g.is_revoked
            and (g.expires_at is None or g.expires_at > datetime.utcnow())
            and entry_id in g.granted_entry_ids
        ]
        
        return len(active_grants) > 0
    
    def get_governance_metrics(self) -> GovernanceMetrics:
        """Get governance system metrics."""
        # Active participants (voted in last 30 days)
        thirty_days_ago = datetime.utcnow() - timedelta(days=30)
        active_count = len([
            p for p in self._participants.values()
            if p.is_active and p.votes_cast > 0
        ])
        
        # Average decision time
        decided_proposals = [
            p for p in self._proposals.values()
            if p.decided_at is not None
        ]
        
        avg_decision_time = 0.0
        if decided_proposals:
            total_time = sum(
                (p.decided_at - p.submitted_at).total_seconds() / 86400
                for p in decided_proposals
                if p.submitted_at
            )
            avg_decision_time = total_time / len(decided_proposals)
        
        # Approval rate
        total_decided = self._metrics.approved_proposals + self._metrics.rejected_proposals
        approval_rate = (
            self._metrics.approved_proposals / total_decided
            if total_decided > 0 else 0.0
        )
        
        self._metrics.active_participants = active_count
        self._metrics.average_decision_time_days = avg_decision_time
        self._metrics.approval_rate = approval_rate
        self._metrics.timestamp = datetime.utcnow()
        
        return self._metrics
