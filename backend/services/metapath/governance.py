"""
MetaPath DAO Governance Module

Trade Secret: Quadratic voting algorithm with TK-weighted multipliers,
proposal eligibility system, impact assessment algorithms.

Key Features:
- Quadratic voting with diminishing returns
- TK holder 2.8x vote multiplier
- Indigenous community 1.8x multiplier
- Validator 1.5x multiplier
- Elder endorsement requirements for TK-sensitive proposals
"""

from enum import Enum
from typing import Dict, List, Optional, Set, Tuple, Any
from pydantic import BaseModel, Field
from datetime import datetime, timedelta
import uuid
import math


class MemberRole(str, Enum):
    """DAO member roles with associated voting multipliers."""
    TK_HOLDER = "tk_holder"          # 2.8x multiplier
    INDIGENOUS_COMMUNITY = "indigenous_community"  # 1.8x multiplier
    VALIDATOR = "validator"          # 1.5x multiplier
    RESEARCHER = "researcher"        # 1.2x multiplier
    DEVELOPER = "developer"          # 1.1x multiplier
    STANDARD = "standard"            # 1.0x multiplier


class ProposalType(str, Enum):
    """Types of governance proposals."""
    PATHWAY_ACTIVATION = "pathway_activation"
    RESOURCE_ALLOCATION = "resource_allocation"
    TK_POLICY = "tk_policy"
    COMPENSATION_DISTRIBUTION = "compensation_distribution"
    PROTOCOL_UPGRADE = "protocol_upgrade"
    EMERGENCY_RESPONSE = "emergency_response"
    PARAMETER_CHANGE = "parameter_change"
    COMMUNITY_GRANT = "community_grant"


class ProposalStatus(str, Enum):
    """Proposal lifecycle status."""
    DRAFT = "draft"
    PENDING_ENDORSEMENT = "pending_endorsement"
    OPEN_VOTING = "open_voting"
    PASSED = "passed"
    REJECTED = "rejected"
    EXECUTED = "executed"
    CANCELLED = "cancelled"
    EXPIRED = "expired"


class ProposalImpactLevel(str, Enum):
    """Impact assessment levels."""
    MINIMAL = "minimal"        # < 1% system impact
    LOW = "low"                # 1-5% system impact
    MODERATE = "moderate"      # 5-15% system impact
    HIGH = "high"              # 15-30% system impact
    CRITICAL = "critical"      # > 30% system impact


class DAOMember(BaseModel):
    """DAO member with roles and voting power."""
    member_id: str = Field(default_factory=lambda: str(uuid.uuid4()))
    wallet_address: str
    display_name: str
    
    # Roles and multipliers
    roles: List[MemberRole] = Field(default_factory=lambda: [MemberRole.STANDARD])
    primary_role: MemberRole = MemberRole.STANDARD
    
    # Token holdings
    governance_tokens: int = 0
    staked_tokens: int = 0
    
    # TK associations
    tk_communities: List[str] = Field(default_factory=list)
    is_tk_holder: bool = False
    tk_contribution_score: float = 0.0  # 0-100
    
    # Validator status
    is_validator: bool = False
    validation_count: int = 0
    validation_accuracy: float = 0.0  # 0-1
    
    # Activity
    joined_at: datetime = Field(default_factory=datetime.utcnow)
    last_active_at: datetime = Field(default_factory=datetime.utcnow)
    proposals_created: int = 0
    votes_cast: int = 0


class ElderEndorsement(BaseModel):
    """Elder endorsement for TK-sensitive proposals."""
    endorsement_id: str = Field(default_factory=lambda: str(uuid.uuid4()))
    proposal_id: str
    elder_member_id: str
    community_name: str
    
    # Endorsement details
    endorsed: bool
    rationale: str = ""
    conditions: List[str] = Field(default_factory=list)
    
    # Metadata
    created_at: datetime = Field(default_factory=datetime.utcnow)


class Vote(BaseModel):
    """Individual vote on a proposal."""
    vote_id: str = Field(default_factory=lambda: str(uuid.uuid4()))
    proposal_id: str
    member_id: str
    
    # Vote details
    support: bool  # True = for, False = against
    tokens_used: int
    
    # Calculated values
    base_vote_power: float = 0.0
    multiplier: float = 1.0
    final_vote_power: float = 0.0
    
    # Metadata
    cast_at: datetime = Field(default_factory=datetime.utcnow)
    rationale: Optional[str] = None


class Proposal(BaseModel):
    """Governance proposal."""
    proposal_id: str = Field(default_factory=lambda: str(uuid.uuid4()))
    title: str
    description: str
    proposal_type: ProposalType
    
    # Creator
    creator_id: str
    creator_name: str = ""
    
    # Status
    status: ProposalStatus = ProposalStatus.DRAFT
    created_at: datetime = Field(default_factory=datetime.utcnow)
    
    # Timing
    voting_starts_at: Optional[datetime] = None
    voting_ends_at: Optional[datetime] = None
    execution_deadline: Optional[datetime] = None
    
    # TK sensitivity
    involves_tk: bool = False
    affected_communities: List[str] = Field(default_factory=list)
    requires_elder_endorsement: bool = False
    elder_endorsements: List[ElderEndorsement] = Field(default_factory=list)
    minimum_endorsements_required: int = 3
    
    # Impact assessment
    impact_level: ProposalImpactLevel = ProposalImpactLevel.MINIMAL
    affected_pathways: List[str] = Field(default_factory=list)
    estimated_resource_impact: Dict[str, float] = Field(default_factory=dict)
    
    # Voting results
    votes: List[Vote] = Field(default_factory=list)
    total_for_power: float = 0.0
    total_against_power: float = 0.0
    quorum_reached: bool = False
    
    # Execution
    execution_data: Dict[str, Any] = Field(default_factory=dict)
    executed_at: Optional[datetime] = None
    execution_tx_hash: Optional[str] = None


class GovernanceConfig(BaseModel):
    """DAO governance configuration."""
    # Voting parameters
    quorum_percentage: float = 0.20  # 20% of voting power
    approval_threshold: float = 0.51  # 51% for passage
    voting_period_days: int = 7
    execution_delay_days: int = 2
    
    # Role multipliers (Trade Secret)
    role_multipliers: Dict[MemberRole, float] = Field(default_factory=lambda: {
        MemberRole.TK_HOLDER: 2.8,
        MemberRole.INDIGENOUS_COMMUNITY: 1.8,
        MemberRole.VALIDATOR: 1.5,
        MemberRole.RESEARCHER: 1.2,
        MemberRole.DEVELOPER: 1.1,
        MemberRole.STANDARD: 1.0,
    })
    
    # Elder endorsement requirements by proposal type
    elder_endorsement_required: Dict[ProposalType, bool] = Field(default_factory=lambda: {
        ProposalType.TK_POLICY: True,
        ProposalType.COMPENSATION_DISTRIBUTION: True,
        ProposalType.PATHWAY_ACTIVATION: False,
        ProposalType.RESOURCE_ALLOCATION: False,
        ProposalType.PROTOCOL_UPGRADE: False,
        ProposalType.EMERGENCY_RESPONSE: False,
        ProposalType.PARAMETER_CHANGE: False,
        ProposalType.COMMUNITY_GRANT: False,
    })


class DAOGovernance:
    """
    DAO Governance Engine for MetaPath.
    
    Trade Secret: Quadratic voting algorithm with culturally-weighted multipliers,
    elder endorsement system, impact assessment algorithms.
    
    Key Multipliers:
    - TK Holder: 2.8x vote weight
    - Indigenous Community: 1.8x vote weight
    - Validator: 1.5x vote weight
    - Researcher: 1.2x vote weight
    - Developer: 1.1x vote weight
    - Standard: 1.0x vote weight
    
    Decentralization Target: 0.91 coefficient
    """
    
    def __init__(self, config: Optional[GovernanceConfig] = None):
        self._config = config or GovernanceConfig()
        self._members: Dict[str, DAOMember] = {}
        self._proposals: Dict[str, Proposal] = {}
        self._total_voting_power: float = 0.0
    
    # ==========================================================================
    # TRADE SECRET: Quadratic Voting Algorithm
    # ==========================================================================
    def _calculate_quadratic_vote_power(
        self,
        tokens_used: int,
    ) -> float:
        """
        Calculate base vote power using quadratic formula.
        Trade Secret: Square root diminishing returns.
        
        Formula: vote_power = sqrt(tokens_used)
        """
        return math.sqrt(tokens_used)
    
    def _calculate_member_multiplier(
        self,
        member: DAOMember,
    ) -> float:
        """
        Calculate vote multiplier based on member roles.
        Trade Secret: Role-weighted multiplier stacking algorithm.
        
        Logic: Uses highest applicable multiplier (no stacking)
        """
        multipliers = self._config.role_multipliers
        
        # Get highest multiplier from member roles
        max_multiplier = 1.0
        for role in member.roles:
            role_mult = multipliers.get(role, 1.0)
            max_multiplier = max(max_multiplier, role_mult)
        
        # Additional TK contribution bonus (up to 0.2x)
        if member.is_tk_holder and member.tk_contribution_score > 0:
            tk_bonus = (member.tk_contribution_score / 100) * 0.2
            max_multiplier += tk_bonus
        
        return max_multiplier
    
    def _calculate_final_vote_power(
        self,
        member: DAOMember,
        tokens_used: int,
    ) -> Tuple[float, float, float]:
        """
        Calculate final vote power with all modifiers.
        Trade Secret: Complete vote power calculation pipeline.
        
        Returns: (base_power, multiplier, final_power)
        """
        base_power = self._calculate_quadratic_vote_power(tokens_used)
        multiplier = self._calculate_member_multiplier(member)
        final_power = base_power * multiplier
        
        return base_power, multiplier, final_power
    
    # ==========================================================================
    # TRADE SECRET: Impact Assessment Algorithm
    # ==========================================================================
    def _assess_proposal_impact(
        self,
        proposal: Proposal,
    ) -> ProposalImpactLevel:
        """
        Assess impact level of proposal.
        Trade Secret: Multi-factor impact scoring algorithm.
        """
        impact_score = 0.0
        
        # Factor 1: Proposal type weight
        type_weights = {
            ProposalType.PATHWAY_ACTIVATION: 25,
            ProposalType.RESOURCE_ALLOCATION: 20,
            ProposalType.TK_POLICY: 30,
            ProposalType.COMPENSATION_DISTRIBUTION: 20,
            ProposalType.PROTOCOL_UPGRADE: 35,
            ProposalType.EMERGENCY_RESPONSE: 40,
            ProposalType.PARAMETER_CHANGE: 15,
            ProposalType.COMMUNITY_GRANT: 10,
        }
        impact_score += type_weights.get(proposal.proposal_type, 10)
        
        # Factor 2: Affected pathways
        pathway_impact = len(proposal.affected_pathways) * 3
        impact_score += min(pathway_impact, 30)  # Cap at 30
        
        # Factor 3: TK involvement
        if proposal.involves_tk:
            impact_score += 15
            impact_score += len(proposal.affected_communities) * 5
        
        # Factor 4: Resource impact
        resource_impact = sum(proposal.estimated_resource_impact.values())
        impact_score += min(resource_impact / 10, 20)  # Cap at 20
        
        # Determine impact level
        if impact_score < 15:
            return ProposalImpactLevel.MINIMAL
        elif impact_score < 35:
            return ProposalImpactLevel.LOW
        elif impact_score < 55:
            return ProposalImpactLevel.MODERATE
        elif impact_score < 75:
            return ProposalImpactLevel.HIGH
        else:
            return ProposalImpactLevel.CRITICAL
    
    # ==========================================================================
    # Public API
    # ==========================================================================
    
    def register_member(
        self,
        wallet_address: str,
        display_name: str,
        governance_tokens: int = 0,
        roles: Optional[List[MemberRole]] = None,
        tk_communities: Optional[List[str]] = None,
        is_validator: bool = False,
    ) -> DAOMember:
        """
        Register a new DAO member.
        """
        member = DAOMember(
            wallet_address=wallet_address,
            display_name=display_name,
            governance_tokens=governance_tokens,
            roles=roles or [MemberRole.STANDARD],
            primary_role=(roles or [MemberRole.STANDARD])[0],
            tk_communities=tk_communities or [],
            is_tk_holder=bool(tk_communities),
            is_validator=is_validator,
        )
        
        # Update total voting power
        base_power = self._calculate_quadratic_vote_power(governance_tokens)
        multiplier = self._calculate_member_multiplier(member)
        self._total_voting_power += base_power * multiplier
        
        self._members[member.member_id] = member
        return member
    
    def update_member_tokens(
        self,
        member_id: str,
        new_token_balance: int,
    ) -> DAOMember:
        """Update member token balance."""
        member = self._members.get(member_id)
        if not member:
            raise ValueError(f"Member not found: {member_id}")
        
        # Recalculate voting power
        old_power = self._calculate_quadratic_vote_power(member.governance_tokens)
        new_power = self._calculate_quadratic_vote_power(new_token_balance)
        multiplier = self._calculate_member_multiplier(member)
        
        self._total_voting_power -= old_power * multiplier
        self._total_voting_power += new_power * multiplier
        
        member.governance_tokens = new_token_balance
        return member
    
    def create_proposal(
        self,
        creator_id: str,
        title: str,
        description: str,
        proposal_type: ProposalType,
        involves_tk: bool = False,
        affected_communities: Optional[List[str]] = None,
        affected_pathways: Optional[List[str]] = None,
        execution_data: Optional[Dict[str, Any]] = None,
    ) -> Proposal:
        """
        Create a new governance proposal.
        """
        creator = self._members.get(creator_id)
        if not creator:
            raise ValueError(f"Creator not found: {creator_id}")
        
        # Check eligibility
        if creator.governance_tokens < 100:
            raise ValueError("Minimum 100 tokens required to create proposal")
        
        proposal = Proposal(
            title=title,
            description=description,
            proposal_type=proposal_type,
            creator_id=creator_id,
            creator_name=creator.display_name,
            involves_tk=involves_tk,
            affected_communities=affected_communities or [],
            affected_pathways=affected_pathways or [],
            execution_data=execution_data or {},
        )
        
        # Determine if elder endorsement required
        if involves_tk or self._config.elder_endorsement_required.get(proposal_type, False):
            proposal.requires_elder_endorsement = True
            proposal.status = ProposalStatus.PENDING_ENDORSEMENT
        
        # Assess impact
        proposal.impact_level = self._assess_proposal_impact(proposal)
        
        # Update creator stats
        creator.proposals_created += 1
        
        self._proposals[proposal.proposal_id] = proposal
        return proposal
    
    def add_elder_endorsement(
        self,
        proposal_id: str,
        elder_member_id: str,
        community_name: str,
        endorsed: bool,
        rationale: str = "",
        conditions: Optional[List[str]] = None,
    ) -> ElderEndorsement:
        """
        Add elder endorsement to proposal.
        """
        proposal = self._proposals.get(proposal_id)
        if not proposal:
            raise ValueError(f"Proposal not found: {proposal_id}")
        
        if not proposal.requires_elder_endorsement:
            raise ValueError("Proposal does not require elder endorsement")
        
        elder = self._members.get(elder_member_id)
        if not elder or MemberRole.TK_HOLDER not in elder.roles:
            raise ValueError("Only TK holders can provide elder endorsements")
        
        endorsement = ElderEndorsement(
            proposal_id=proposal_id,
            elder_member_id=elder_member_id,
            community_name=community_name,
            endorsed=endorsed,
            rationale=rationale,
            conditions=conditions or [],
        )
        
        proposal.elder_endorsements.append(endorsement)
        
        # Check if minimum endorsements reached
        positive_endorsements = sum(
            1 for e in proposal.elder_endorsements if e.endorsed
        )
        
        if positive_endorsements >= proposal.minimum_endorsements_required:
            # Move to voting
            proposal.status = ProposalStatus.OPEN_VOTING
            proposal.voting_starts_at = datetime.utcnow()
            proposal.voting_ends_at = datetime.utcnow() + timedelta(
                days=self._config.voting_period_days
            )
        
        return endorsement
    
    def open_voting(
        self,
        proposal_id: str,
    ) -> Proposal:
        """
        Open voting on a proposal (for non-TK proposals).
        """
        proposal = self._proposals.get(proposal_id)
        if not proposal:
            raise ValueError(f"Proposal not found: {proposal_id}")
        
        if proposal.requires_elder_endorsement:
            raise ValueError("Proposal requires elder endorsement before voting")
        
        proposal.status = ProposalStatus.OPEN_VOTING
        proposal.voting_starts_at = datetime.utcnow()
        proposal.voting_ends_at = datetime.utcnow() + timedelta(
            days=self._config.voting_period_days
        )
        
        return proposal
    
    def cast_vote(
        self,
        proposal_id: str,
        member_id: str,
        support: bool,
        tokens_to_use: int,
        rationale: Optional[str] = None,
    ) -> Vote:
        """
        Cast vote on a proposal using quadratic voting.
        """
        proposal = self._proposals.get(proposal_id)
        if not proposal:
            raise ValueError(f"Proposal not found: {proposal_id}")
        
        if proposal.status != ProposalStatus.OPEN_VOTING:
            raise ValueError(f"Proposal not open for voting: {proposal.status}")
        
        member = self._members.get(member_id)
        if not member:
            raise ValueError(f"Member not found: {member_id}")
        
        if tokens_to_use > member.governance_tokens:
            raise ValueError("Insufficient tokens")
        
        # Check for duplicate vote
        existing_vote = next(
            (v for v in proposal.votes if v.member_id == member_id),
            None
        )
        if existing_vote:
            raise ValueError("Member has already voted on this proposal")
        
        # Calculate vote power
        base_power, multiplier, final_power = self._calculate_final_vote_power(
            member, tokens_to_use
        )
        
        vote = Vote(
            proposal_id=proposal_id,
            member_id=member_id,
            support=support,
            tokens_used=tokens_to_use,
            base_vote_power=base_power,
            multiplier=multiplier,
            final_vote_power=final_power,
            rationale=rationale,
        )
        
        proposal.votes.append(vote)
        
        # Update totals
        if support:
            proposal.total_for_power += final_power
        else:
            proposal.total_against_power += final_power
        
        # Check quorum
        total_power_voted = proposal.total_for_power + proposal.total_against_power
        if total_power_voted >= self._total_voting_power * self._config.quorum_percentage:
            proposal.quorum_reached = True
        
        # Update member stats
        member.votes_cast += 1
        member.last_active_at = datetime.utcnow()
        
        return vote
    
    def finalize_proposal(
        self,
        proposal_id: str,
    ) -> Proposal:
        """
        Finalize voting and determine outcome.
        """
        proposal = self._proposals.get(proposal_id)
        if not proposal:
            raise ValueError(f"Proposal not found: {proposal_id}")
        
        if proposal.status != ProposalStatus.OPEN_VOTING:
            raise ValueError(f"Proposal not in voting state: {proposal.status}")
        
        # Check if voting period ended
        if proposal.voting_ends_at and datetime.utcnow() < proposal.voting_ends_at:
            raise ValueError("Voting period has not ended")
        
        # Determine outcome
        if not proposal.quorum_reached:
            proposal.status = ProposalStatus.EXPIRED
        else:
            total_power = proposal.total_for_power + proposal.total_against_power
            approval_rate = proposal.total_for_power / total_power if total_power > 0 else 0
            
            if approval_rate >= self._config.approval_threshold:
                proposal.status = ProposalStatus.PASSED
                proposal.execution_deadline = datetime.utcnow() + timedelta(
                    days=self._config.execution_delay_days
                )
            else:
                proposal.status = ProposalStatus.REJECTED
        
        return proposal
    
    def execute_proposal(
        self,
        proposal_id: str,
        tx_hash: Optional[str] = None,
    ) -> Proposal:
        """
        Mark proposal as executed.
        """
        proposal = self._proposals.get(proposal_id)
        if not proposal:
            raise ValueError(f"Proposal not found: {proposal_id}")
        
        if proposal.status != ProposalStatus.PASSED:
            raise ValueError(f"Proposal not passed: {proposal.status}")
        
        proposal.status = ProposalStatus.EXECUTED
        proposal.executed_at = datetime.utcnow()
        proposal.execution_tx_hash = tx_hash
        
        return proposal
    
    def get_proposal(self, proposal_id: str) -> Optional[Proposal]:
        """Get proposal by ID."""
        return self._proposals.get(proposal_id)
    
    def get_member(self, member_id: str) -> Optional[DAOMember]:
        """Get member by ID."""
        return self._members.get(member_id)
    
    def get_active_proposals(self) -> List[Proposal]:
        """Get all active proposals."""
        active_statuses = {
            ProposalStatus.DRAFT,
            ProposalStatus.PENDING_ENDORSEMENT,
            ProposalStatus.OPEN_VOTING,
            ProposalStatus.PASSED,
        }
        return [
            p for p in self._proposals.values()
            if p.status in active_statuses
        ]
    
    def get_member_voting_power(
        self,
        member_id: str,
    ) -> Dict[str, Any]:
        """Get member's voting power breakdown."""
        member = self._members.get(member_id)
        if not member:
            return {}
        
        base_power = self._calculate_quadratic_vote_power(member.governance_tokens)
        multiplier = self._calculate_member_multiplier(member)
        
        return {
            "member_id": member_id,
            "governance_tokens": member.governance_tokens,
            "base_vote_power": base_power,
            "multiplier": multiplier,
            "max_vote_power": base_power * multiplier,
            "roles": [r.value for r in member.roles],
            "is_tk_holder": member.is_tk_holder,
            "tk_contribution_bonus": (member.tk_contribution_score / 100) * 0.2 if member.is_tk_holder else 0,
        }
    
    def calculate_decentralization_coefficient(self) -> float:
        """
        Calculate decentralization coefficient.
        Trade Secret: Gini-based decentralization metric.
        
        Target: 0.91 (higher = more decentralized)
        """
        if not self._members:
            return 0.0
        
        # Get voting powers
        powers = []
        for member in self._members.values():
            base_power = self._calculate_quadratic_vote_power(member.governance_tokens)
            multiplier = self._calculate_member_multiplier(member)
            powers.append(base_power * multiplier)
        
        if not powers or sum(powers) == 0:
            return 0.0
        
        # Calculate Gini coefficient
        n = len(powers)
        sorted_powers = sorted(powers)
        
        gini_sum = sum((2 * (i + 1) - n - 1) * x for i, x in enumerate(sorted_powers))
        gini = gini_sum / (n * sum(powers))
        
        # Invert: 1 - Gini gives decentralization (1 = perfect equality)
        decentralization = 1 - gini
        
        return round(decentralization, 3)
    
    def get_governance_stats(self) -> Dict[str, Any]:
        """Get overall governance statistics."""
        return {
            "total_members": len(self._members),
            "total_voting_power": self._total_voting_power,
            "total_proposals": len(self._proposals),
            "active_proposals": len(self.get_active_proposals()),
            "decentralization_coefficient": self.calculate_decentralization_coefficient(),
            "config": {
                "quorum_percentage": self._config.quorum_percentage,
                "approval_threshold": self._config.approval_threshold,
                "voting_period_days": self._config.voting_period_days,
                "role_multipliers": {
                    k.value: v for k, v in self._config.role_multipliers.items()
                },
            },
        }
