"""
EthnoPath Validator - Community Validation Networks

Trade Secret: Multi-stakeholder consensus algorithms, validation scoring,
accuracy thresholds, healer and elder network validation protocols.

Core Capabilities:
- Distributed verification of traditional knowledge accuracy
- Multi-stakeholder consensus requirements
- Cultural appropriateness validation
- Continuous validation and update cycles
"""

from enum import Enum
from typing import Dict, List, Optional, Set, Any
from pydantic import BaseModel, Field
from datetime import datetime, timedelta
import uuid


class ValidationLevel(str, Enum):
    """Validation level requirements."""
    BASIC = "basic"                    # Single community member
    COMMUNITY = "community"            # Multiple community members
    ELDER_COUNCIL = "elder_council"    # Elder council approval
    MULTI_COMMUNITY = "multi_community"  # Multiple communities agree
    SACRED_PROTOCOL = "sacred_protocol"  # Special ceremonial validation


class ValidatorRole(str, Enum):
    """Validator role in community."""
    COMMUNITY_MEMBER = "community_member"
    KNOWLEDGE_HOLDER = "knowledge_holder"
    ELDER = "elder"
    HEALER = "healer"
    CEREMONIAL_LEADER = "ceremonial_leader"
    CULTURAL_EXPERT = "cultural_expert"
    YOUTH_REPRESENTATIVE = "youth_representative"


class ConsensusType(str, Enum):
    """Consensus mechanism type."""
    SIMPLE_MAJORITY = "simple_majority"      # >50%
    SUPERMAJORITY = "supermajority"          # >66%
    UNANIMOUS = "unanimous"                  # 100%
    ELDER_WEIGHTED = "elder_weighted"        # Elders get 2x weight
    HEALER_VETO = "healer_veto"              # Healers can veto


class ValidationStatus(str, Enum):
    """Validation status."""
    PENDING = "pending"
    IN_PROGRESS = "in_progress"
    APPROVED = "approved"
    REJECTED = "rejected"
    NEEDS_REVISION = "needs_revision"
    DISPUTED = "disputed"
    WITHDRAWN = "withdrawn"


class ValidatorProfile(BaseModel):
    """Community validator profile."""
    validator_id: str = Field(default_factory=lambda: str(uuid.uuid4()))
    community_id: str
    
    # Identity
    display_name: str
    role: ValidatorRole
    
    # Expertise
    expertise_domains: List[str] = Field(default_factory=list)
    years_of_knowledge: Optional[int] = None
    
    # Validation history
    validations_performed: int = 0
    validations_approved: int = 0
    validations_rejected: int = 0
    accuracy_score: float = 1.0  # 0-1 scale
    
    # Status
    is_active: bool = True
    last_validation: Optional[datetime] = None
    joined_date: datetime = Field(default_factory=datetime.utcnow)


class ValidationVote(BaseModel):
    """Individual validation vote."""
    vote_id: str = Field(default_factory=lambda: str(uuid.uuid4()))
    round_id: str
    validator_id: str
    
    # Vote
    approved: bool
    confidence_score: float = Field(ge=0.0, le=1.0)  # How confident in validation
    
    # Feedback
    comments: str = ""
    suggested_revisions: List[str] = Field(default_factory=list)
    cultural_concerns: List[str] = Field(default_factory=list)
    
    # Context
    voted_at: datetime = Field(default_factory=datetime.utcnow)
    validator_role: ValidatorRole = ValidatorRole.COMMUNITY_MEMBER


class ValidationRound(BaseModel):
    """Validation round for knowledge entry."""
    round_id: str = Field(default_factory=lambda: str(uuid.uuid4()))
    entry_id: str
    community_id: str
    
    # Round configuration
    validation_level: ValidationLevel
    consensus_type: ConsensusType
    minimum_validators: int = 3
    
    # Votes
    votes: List[ValidationVote] = Field(default_factory=list)
    invited_validators: List[str] = Field(default_factory=list)
    
    # Status
    status: ValidationStatus = ValidationStatus.PENDING
    started_at: datetime = Field(default_factory=datetime.utcnow)
    ended_at: Optional[datetime] = None
    
    # Results
    approval_percentage: Optional[float] = None
    consensus_reached: bool = False
    final_decision: Optional[bool] = None
    decision_rationale: str = ""


class CommunityConsensus(BaseModel):
    """Multi-community consensus result."""
    consensus_id: str = Field(default_factory=lambda: str(uuid.uuid4()))
    entry_id: str
    
    # Participating communities
    communities_invited: List[str] = Field(default_factory=list)
    communities_responded: List[str] = Field(default_factory=list)
    
    # Consensus results
    communities_approved: List[str] = Field(default_factory=list)
    communities_rejected: List[str] = Field(default_factory=list)
    
    # Final outcome
    consensus_reached: bool = False
    approval_percentage: float = 0.0
    final_decision: Optional[bool] = None
    
    # Timing
    initiated_at: datetime = Field(default_factory=datetime.utcnow)
    completed_at: Optional[datetime] = None


class ValidationMetrics(BaseModel):
    """Validation system metrics."""
    timestamp: datetime = Field(default_factory=datetime.utcnow)
    
    # Validation volume
    total_validations: int = 0
    active_validations: int = 0
    completed_validations: int = 0
    
    # Outcomes
    approved_count: int = 0
    rejected_count: int = 0
    disputed_count: int = 0
    
    # Quality
    average_approval_time_hours: float = 0.0
    average_validator_confidence: float = 0.0
    consensus_achievement_rate: float = 0.0
    
    # Participation
    active_validators: int = 0
    average_validators_per_round: float = 0.0


class CommunityValidator:
    """
    Community Validation Networks for Traditional Knowledge.
    
    Trade Secret: Multi-stakeholder consensus algorithms, validation scoring,
    cultural appropriateness checks, healer/elder network protocols.
    
    Key Features:
    - Distributed verification of TK accuracy
    - Multi-stakeholder consensus requirements
    - Cultural context validation
    - Continuous validation cycles with community input
    """
    
    # ==========================================================================
    # TRADE SECRET: Validation Level Requirements Matrix
    # ==========================================================================
    _VALIDATION_REQUIREMENTS = {
        ValidationLevel.BASIC: {
            "min_validators": 1,
            "required_roles": [ValidatorRole.COMMUNITY_MEMBER],
            "consensus_type": ConsensusType.SIMPLE_MAJORITY,
        },
        ValidationLevel.COMMUNITY: {
            "min_validators": 3,
            "required_roles": [ValidatorRole.COMMUNITY_MEMBER, ValidatorRole.KNOWLEDGE_HOLDER],
            "consensus_type": ConsensusType.SIMPLE_MAJORITY,
        },
        ValidationLevel.ELDER_COUNCIL: {
            "min_validators": 2,
            "required_roles": [ValidatorRole.ELDER],
            "consensus_type": ConsensusType.UNANIMOUS,
        },
        ValidationLevel.MULTI_COMMUNITY: {
            "min_validators": 5,
            "min_communities": 3,
            "required_roles": [ValidatorRole.ELDER, ValidatorRole.KNOWLEDGE_HOLDER],
            "consensus_type": ConsensusType.SUPERMAJORITY,
        },
        ValidationLevel.SACRED_PROTOCOL: {
            "min_validators": 3,
            "required_roles": [ValidatorRole.ELDER, ValidatorRole.CEREMONIAL_LEADER],
            "consensus_type": ConsensusType.UNANIMOUS,
            "requires_ceremony": True,
        },
    }
    
    # ==========================================================================
    # TRADE SECRET: Consensus Threshold Configuration
    # ==========================================================================
    _CONSENSUS_THRESHOLDS = {
        ConsensusType.SIMPLE_MAJORITY: 0.50,
        ConsensusType.SUPERMAJORITY: 0.66,
        ConsensusType.UNANIMOUS: 1.00,
        ConsensusType.ELDER_WEIGHTED: 0.60,  # With elder weighting
        ConsensusType.HEALER_VETO: 0.50,     # Base, but healers can override
    }
    
    # ==========================================================================
    # TRADE SECRET: Role Voting Weights
    # ==========================================================================
    _ROLE_WEIGHTS = {
        ValidatorRole.COMMUNITY_MEMBER: 1.0,
        ValidatorRole.KNOWLEDGE_HOLDER: 1.2,
        ValidatorRole.ELDER: 2.0,
        ValidatorRole.HEALER: 1.8,
        ValidatorRole.CEREMONIAL_LEADER: 2.5,
        ValidatorRole.CULTURAL_EXPERT: 1.5,
        ValidatorRole.YOUTH_REPRESENTATIVE: 1.0,
    }
    
    def __init__(self):
        self._validators: Dict[str, ValidatorProfile] = {}
        self._rounds: Dict[str, ValidationRound] = {}
        self._consensus_records: Dict[str, CommunityConsensus] = {}
        self._metrics = ValidationMetrics()
    
    # ==========================================================================
    # Validator Management
    # ==========================================================================
    
    def register_validator(
        self,
        community_id: str,
        display_name: str,
        role: ValidatorRole,
        expertise_domains: Optional[List[str]] = None,
        years_of_knowledge: Optional[int] = None,
    ) -> ValidatorProfile:
        """Register a community validator."""
        validator = ValidatorProfile(
            community_id=community_id,
            display_name=display_name,
            role=role,
            expertise_domains=expertise_domains or [],
            years_of_knowledge=years_of_knowledge,
        )
        
        self._validators[validator.validator_id] = validator
        return validator
    
    # ==========================================================================
    # TRADE SECRET: Validation Round Orchestration
    # ==========================================================================
    
    def initiate_validation_round(
        self,
        entry_id: str,
        community_id: str,
        validation_level: ValidationLevel,
        invited_validators: Optional[List[str]] = None,
    ) -> ValidationRound:
        """
        Initiate validation round for knowledge entry.
        Trade Secret: Automatic validator selection, role requirements.
        """
        requirements = self._VALIDATION_REQUIREMENTS.get(validation_level, {})
        
        # Auto-select validators if not provided
        if not invited_validators:
            invited_validators = self._select_validators(
                community_id,
                validation_level,
                requirements,
            )
        
        round_obj = ValidationRound(
            entry_id=entry_id,
            community_id=community_id,
            validation_level=validation_level,
            consensus_type=requirements.get("consensus_type", ConsensusType.SIMPLE_MAJORITY),
            minimum_validators=requirements.get("min_validators", 3),
            invited_validators=invited_validators,
            status=ValidationStatus.IN_PROGRESS,
        )
        
        self._rounds[round_obj.round_id] = round_obj
        self._metrics.total_validations += 1
        self._metrics.active_validations += 1
        
        return round_obj
    
    def _select_validators(
        self,
        community_id: str,
        validation_level: ValidationLevel,
        requirements: Dict[str, Any],
    ) -> List[str]:
        """
        Auto-select validators based on requirements.
        Trade Secret: Validator selection algorithm.
        """
        required_roles = requirements.get("required_roles", [])
        min_validators = requirements.get("min_validators", 3)
        
        # Get eligible validators from community
        eligible = [
            v for v in self._validators.values()
            if v.community_id == community_id
            and v.is_active
            and v.role in required_roles
        ]
        
        # Sort by accuracy score and experience
        eligible.sort(key=lambda v: (v.accuracy_score, v.validations_performed), reverse=True)
        
        # Select top validators
        selected = [v.validator_id for v in eligible[:min_validators]]
        
        return selected
    
    def submit_validation_vote(
        self,
        round_id: str,
        validator_id: str,
        approved: bool,
        confidence_score: float,
        comments: str = "",
        suggested_revisions: Optional[List[str]] = None,
        cultural_concerns: Optional[List[str]] = None,
    ) -> ValidationVote:
        """Submit validation vote."""
        round_obj = self._rounds.get(round_id)
        if not round_obj:
            raise ValueError(f"Validation round not found: {round_id}")
        
        if round_obj.status != ValidationStatus.IN_PROGRESS:
            raise ValueError(f"Round not accepting votes: {round_obj.status}")
        
        # Verify validator is invited
        if validator_id not in round_obj.invited_validators:
            raise ValueError(f"Validator not invited to this round: {validator_id}")
        
        # Get validator profile
        validator = self._validators.get(validator_id)
        if not validator:
            raise ValueError(f"Validator not found: {validator_id}")
        
        vote = ValidationVote(
            round_id=round_id,
            validator_id=validator_id,
            approved=approved,
            confidence_score=confidence_score,
            comments=comments,
            suggested_revisions=suggested_revisions or [],
            cultural_concerns=cultural_concerns or [],
            validator_role=validator.role,
        )
        
        round_obj.votes.append(vote)
        
        # Update validator stats
        validator.validations_performed += 1
        if approved:
            validator.validations_approved += 1
        else:
            validator.validations_rejected += 1
        validator.last_validation = datetime.utcnow()
        
        # Check if round complete
        self._check_round_completion(round_obj)
        
        return vote
    
    # ==========================================================================
    # TRADE SECRET: Consensus Calculation Algorithm
    # ==========================================================================
    
    def _check_round_completion(self, round_obj: ValidationRound) -> None:
        """
        Check if validation round has reached consensus.
        Trade Secret: Weighted consensus algorithm, role-based voting.
        """
        # Check if minimum validators have voted
        if len(round_obj.votes) < round_obj.minimum_validators:
            return  # Not enough votes yet
        
        # Calculate consensus based on type
        consensus_type = round_obj.consensus_type
        threshold = self._CONSENSUS_THRESHOLDS.get(consensus_type, 0.50)
        
        if consensus_type == ConsensusType.ELDER_WEIGHTED:
            approval_rate = self._calculate_weighted_approval(round_obj.votes)
        elif consensus_type == ConsensusType.HEALER_VETO:
            approval_rate = self._calculate_healer_veto_approval(round_obj.votes)
        else:
            # Simple vote count
            approved_votes = len([v for v in round_obj.votes if v.approved])
            approval_rate = approved_votes / len(round_obj.votes)
        
        round_obj.approval_percentage = approval_rate
        
        # Check consensus threshold
        if approval_rate >= threshold:
            round_obj.consensus_reached = True
            round_obj.final_decision = True
            round_obj.status = ValidationStatus.APPROVED
            round_obj.decision_rationale = f"Consensus reached: {approval_rate:.1%} approval"
            self._metrics.approved_count += 1
        elif len(round_obj.votes) == len(round_obj.invited_validators):
            # All invited validators have voted, no consensus
            round_obj.consensus_reached = False
            round_obj.final_decision = False
            round_obj.status = ValidationStatus.REJECTED
            round_obj.decision_rationale = f"Consensus not reached: {approval_rate:.1%} approval (needed {threshold:.1%})"
            self._metrics.rejected_count += 1
        else:
            # Still waiting for more votes
            return
        
        # Round complete
        round_obj.ended_at = datetime.utcnow()
        self._metrics.active_validations -= 1
        self._metrics.completed_validations += 1
    
    def _calculate_weighted_approval(self, votes: List[ValidationVote]) -> float:
        """
        Calculate approval with elder weighting.
        Trade Secret: Role-based vote weighting algorithm.
        """
        total_weight = 0.0
        approved_weight = 0.0
        
        for vote in votes:
            weight = self._ROLE_WEIGHTS.get(vote.validator_role, 1.0)
            total_weight += weight
            if vote.approved:
                approved_weight += weight
        
        return approved_weight / total_weight if total_weight > 0 else 0.0
    
    def _calculate_healer_veto_approval(self, votes: List[ValidationVote]) -> float:
        """
        Calculate approval with healer veto.
        Trade Secret: Healer veto logic.
        """
        # Check for healer veto
        healer_votes = [v for v in votes if v.validator_role == ValidatorRole.HEALER]
        if any(not v.approved for v in healer_votes):
            return 0.0  # Healer vetoed
        
        # Standard majority
        approved_votes = len([v for v in votes if v.approved])
        return approved_votes / len(votes)
    
    # ==========================================================================
    # Multi-Community Consensus
    # ==========================================================================
    
    def initiate_multi_community_consensus(
        self,
        entry_id: str,
        communities_invited: List[str],
    ) -> CommunityConsensus:
        """Initiate multi-community consensus validation."""
        consensus = CommunityConsensus(
            entry_id=entry_id,
            communities_invited=communities_invited,
        )
        
        self._consensus_records[consensus.consensus_id] = consensus
        
        # Initiate validation rounds in each community
        for community_id in communities_invited:
            self.initiate_validation_round(
                entry_id=entry_id,
                community_id=community_id,
                validation_level=ValidationLevel.COMMUNITY,
            )
        
        return consensus
    
    def update_multi_community_consensus(
        self,
        consensus_id: str,
    ) -> CommunityConsensus:
        """Update multi-community consensus based on round results."""
        consensus = self._consensus_records.get(consensus_id)
        if not consensus:
            raise ValueError(f"Consensus not found: {consensus_id}")
        
        # Find related rounds
        related_rounds = [
            r for r in self._rounds.values()
            if r.entry_id == consensus.entry_id
            and r.community_id in consensus.communities_invited
        ]
        
        # Update consensus status
        consensus.communities_responded = []
        consensus.communities_approved = []
        consensus.communities_rejected = []
        
        for round_obj in related_rounds:
            if round_obj.status in [ValidationStatus.APPROVED, ValidationStatus.REJECTED]:
                consensus.communities_responded.append(round_obj.community_id)
                
                if round_obj.status == ValidationStatus.APPROVED:
                    consensus.communities_approved.append(round_obj.community_id)
                else:
                    consensus.communities_rejected.append(round_obj.community_id)
        
        # Calculate approval percentage
        if len(consensus.communities_responded) > 0:
            consensus.approval_percentage = (
                len(consensus.communities_approved) / len(consensus.communities_responded)
            )
        
        # Check if consensus reached (supermajority of communities)
        if len(consensus.communities_responded) >= len(consensus.communities_invited):
            if consensus.approval_percentage >= 0.66:  # Supermajority
                consensus.consensus_reached = True
                consensus.final_decision = True
            else:
                consensus.consensus_reached = False
                consensus.final_decision = False
            
            consensus.completed_at = datetime.utcnow()
        
        return consensus
    
    def get_validation_metrics(self) -> ValidationMetrics:
        """Get validation system metrics."""
        # Update metrics
        active_validators = len([v for v in self._validators.values() if v.is_active])
        
        # Average confidence
        all_votes = []
        for round_obj in self._rounds.values():
            all_votes.extend(round_obj.votes)
        
        avg_confidence = (
            sum(v.confidence_score for v in all_votes) / len(all_votes)
            if all_votes else 0.0
        )
        
        # Consensus achievement rate
        completed_rounds = [r for r in self._rounds.values() if r.ended_at is not None]
        consensus_achieved = len([r for r in completed_rounds if r.consensus_reached])
        consensus_rate = (
            consensus_achieved / len(completed_rounds)
            if completed_rounds else 0.0
        )
        
        self._metrics.active_validators = active_validators
        self._metrics.average_validator_confidence = avg_confidence
        self._metrics.consensus_achievement_rate = consensus_rate
        self._metrics.timestamp = datetime.utcnow()
        
        return self._metrics
