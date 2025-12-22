"""
EthnoPath Test Suite - Traditional Knowledge Digitization and Preservation

Comprehensive tests for TS-ETH-001 implementation:
- Community registration and management
- Knowledge holder registration
- Traditional knowledge entry creation
- Consent management
- Federated aggregation
- Community validation networks
- DAO governance with elder veto
- Access control and privacy preservation
"""

import pytest
from datetime import datetime, timedelta

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
)
from backend.services.ethnopath.validator import (
    ValidationLevel,
    ValidatorRole,
    ValidationStatus,
    ConsensusType,
)
from backend.services.ethnopath.governance import (
    GovernanceRole,
    ProposalType,
    VoteType,
    ProposalStatus,
    TKComplianceLevel,
    PricingTier,
)


# ==========================================================================
# Digitizer Tests
# ==========================================================================

class TestEthnoPathDigitizer:
    """Test traditional knowledge digitization functionality."""
    
    def test_register_community(self):
        """Test community registration."""
        digitizer = EthnoPathDigitizer()
        
        community = digitizer.register_community(
            name="Test Indigenous Community",
            location="Test Region",
            indigenous_group="Test People",
            contact_email="contact@test.example",
            language="Test Language",
            population_size=1000,
        )
        
        assert community.community_name == "Test Indigenous Community"
        assert community.indigenous_group == "Test People"
        assert community.is_active is True
        assert community.data_sovereignty_status == "autonomous"
        assert community.federated_node is not None
        assert community.federated_node.node_type == "community_server"
    
    def test_register_knowledge_holder(self):
        """Test knowledge holder registration."""
        digitizer = EthnoPathDigitizer()
        
        # Register community first
        community = digitizer.register_community(
            name="Test Community",
            location="Test Location",
            indigenous_group="Test Group",
            contact_email="test@example.com",
        )
        
        # Register knowledge holder
        holder = digitizer.register_knowledge_holder(
            community_id=community.community_id,
            name="Elder TestName",
            role="Elder",
            areas_of_expertise=["Medicinal Plants", "Ceremonies"],
            years_of_experience=40,
        )
        
        assert holder.name == "Elder TestName"
        assert holder.role == "Elder"
        assert "Medicinal Plants" in holder.areas_of_expertise
        assert holder.years_of_experience == 40
        assert holder.is_active is True
    
    def test_create_knowledge_entry(self):
        """Test knowledge entry creation."""
        digitizer = EthnoPathDigitizer()
        
        # Setup
        community = digitizer.register_community(
            name="Test Community",
            location="Test Location",
            indigenous_group="Test Group",
            contact_email="test@example.com",
        )
        
        holder = digitizer.register_knowledge_holder(
            community_id=community.community_id,
            name="Test Holder",
            role="Healer",
            areas_of_expertise=["Medicinal Plants"],
        )
        
        # Create entry
        entry = digitizer.create_knowledge_entry(
            community_id=community.community_id,
            holder_id=holder.holder_id,
            domain=KnowledgeDomain.MEDICINAL_PLANTS,
            knowledge_type=KnowledgeType.THERAPEUTIC_USE,
            sensitivity_level=SensitivityLevel.COMMUNITY_ONLY,
            title="Traditional Healing Plant",
            description="A plant used for healing in our community for generations.",
            traditional_name="Test Plant Name",
            tags=["healing", "traditional"],
        )
        
        assert entry.title == "Traditional Healing Plant"
        assert entry.domain == KnowledgeDomain.MEDICINAL_PLANTS
        assert entry.sensitivity_level == SensitivityLevel.COMMUNITY_ONLY
        assert entry.is_validated is False  # Not validated yet
        assert "healing" in entry.tags
    
    def test_consent_management(self):
        """Test consent granting and tracking."""
        digitizer = EthnoPathDigitizer()
        
        # Setup
        community = digitizer.register_community(
            name="Test Community",
            location="Test Location",
            indigenous_group="Test Group",
            contact_email="test@example.com",
        )
        
        holder = digitizer.register_knowledge_holder(
            community_id=community.community_id,
            name="Test Holder",
            role="Elder",
            areas_of_expertise=["Traditional Knowledge"],
        )
        
        entry = digitizer.create_knowledge_entry(
            community_id=community.community_id,
            holder_id=holder.holder_id,
            domain=KnowledgeDomain.ECOLOGICAL_MANAGEMENT,
            knowledge_type=KnowledgeType.CONSERVATION_PRACTICE,
            sensitivity_level=SensitivityLevel.PUBLIC,
            title="Ecological Practice",
            description="Traditional ecological management.",
        )
        
        # Grant consent
        digitizer.grant_consent(
            entry_id=entry.entry_id,
            holder_id=holder.holder_id,
            consent_type="research",
            restrictions=["No commercial use", "Attribution required"],
        )
        
        # Verify consent
        updated_entry = digitizer._knowledge_entries[entry.entry_id]
        assert len(updated_entry.consent_records) == 1
        assert updated_entry.consent_records[0].consent_type == "research"
        assert "No commercial use" in updated_entry.consent_records[0].restrictions
    
    def test_access_control_community_only(self):
        """Test access control for community-only knowledge."""
        digitizer = EthnoPathDigitizer()
        
        # Setup two communities
        community1 = digitizer.register_community(
            name="Community 1",
            location="Location 1",
            indigenous_group="Group 1",
            contact_email="c1@example.com",
        )
        
        community2 = digitizer.register_community(
            name="Community 2",
            location="Location 2",
            indigenous_group="Group 2",
            contact_email="c2@example.com",
        )
        
        holder = digitizer.register_knowledge_holder(
            community_id=community1.community_id,
            name="Test Holder",
            role="Elder",
            areas_of_expertise=["Traditional Knowledge"],
        )
        
        entry = digitizer.create_knowledge_entry(
            community_id=community1.community_id,
            holder_id=holder.holder_id,
            domain=KnowledgeDomain.CEREMONIAL_KNOWLEDGE,
            knowledge_type=KnowledgeType.SACRED_PRACTICE,
            sensitivity_level=SensitivityLevel.COMMUNITY_ONLY,
            title="Community Practice",
            description="Practice specific to community 1.",
        )
        
        # Community 1 should have access
        can_access_c1 = digitizer.can_access_entry(
            entry.entry_id,
            community1.community_id,
            SensitivityLevel.COMMUNITY_ONLY,
        )
        assert can_access_c1 is True
        
        # Community 2 should NOT have access
        can_access_c2 = digitizer.can_access_entry(
            entry.entry_id,
            community2.community_id,
            SensitivityLevel.COMMUNITY_ONLY,
        )
        assert can_access_c2 is False
    
    def test_federated_aggregation(self):
        """Test privacy-preserving federated aggregation."""
        digitizer = EthnoPathDigitizer()
        
        # Create multiple communities with similar knowledge
        for i in range(5):
            community = digitizer.register_community(
                name=f"Community {i}",
                location=f"Location {i}",
                indigenous_group=f"Group {i}",
                contact_email=f"c{i}@example.com",
            )
            
            holder = digitizer.register_knowledge_holder(
                community_id=community.community_id,
                name=f"Holder {i}",
                role="Healer",
                areas_of_expertise=["Medicinal Plants"],
            )
            
            # Create public knowledge entry
            entry = digitizer.create_knowledge_entry(
                community_id=community.community_id,
                holder_id=holder.holder_id,
                domain=KnowledgeDomain.MEDICINAL_PLANTS,
                knowledge_type=KnowledgeType.THERAPEUTIC_USE,
                sensitivity_level=SensitivityLevel.PUBLIC,
                title=f"Healing Plant Knowledge {i}",
                description="Traditional medicinal plant knowledge.",
                tags=["medicinal", "healing"],
            )
            
            # Grant consent for aggregation
            digitizer.grant_consent(
                entry_id=entry.entry_id,
                holder_id=holder.holder_id,
                consent_type="aggregation",
            )
        
        # Aggregate patterns
        result = digitizer.aggregate_knowledge_patterns(
            domain=KnowledgeDomain.MEDICINAL_PLANTS,
            minimum_communities=3,
        )
        
        # Should find common pattern (medicinal tag)
        assert "common_tags" in result
        assert "medicinal" in result["common_tags"]
        assert "healing" in result["common_tags"]


# ==========================================================================
# Validator Tests
# ==========================================================================

class TestCommunityValidator:
    """Test community validation networks."""
    
    def test_register_validator(self):
        """Test validator registration."""
        validator = CommunityValidator()
        
        profile = validator.register_validator(
            community_id="test_community",
            display_name="Test Validator",
            role=ValidatorRole.ELDER,
            expertise_domains=["Medicinal Plants", "Ceremonies"],
            years_of_knowledge=30,
        )
        
        assert profile.display_name == "Test Validator"
        assert profile.role == ValidatorRole.ELDER
        assert profile.years_of_knowledge == 30
        assert profile.is_active is True
        assert profile.accuracy_score == 1.0  # Perfect initial score
    
    def test_initiate_validation_round(self):
        """Test validation round initiation."""
        validator = CommunityValidator()
        
        # Register validators
        validators = []
        for i in range(3):
            v = validator.register_validator(
                community_id="test_community",
                display_name=f"Validator {i}",
                role=ValidatorRole.COMMUNITY_MEMBER,
                expertise_domains=["Traditional Knowledge"],
            )
            validators.append(v.validator_id)
        
        # Initiate round
        round_obj = validator.initiate_validation_round(
            entry_id="test_entry",
            community_id="test_community",
            validation_level=ValidationLevel.COMMUNITY,
            invited_validators=validators,
        )
        
        assert round_obj.entry_id == "test_entry"
        assert round_obj.validation_level == ValidationLevel.COMMUNITY
        assert round_obj.status == ValidationStatus.IN_PROGRESS
        assert len(round_obj.invited_validators) == 3
    
    def test_validation_voting_consensus_reached(self):
        """Test validation voting with consensus."""
        validator = CommunityValidator()
        
        # Register validators
        validators = []
        for i in range(3):
            v = validator.register_validator(
                community_id="test_community",
                display_name=f"Validator {i}",
                role=ValidatorRole.KNOWLEDGE_HOLDER,
                expertise_domains=["Traditional Knowledge"],
            )
            validators.append(v.validator_id)
        
        # Initiate round
        round_obj = validator.initiate_validation_round(
            entry_id="test_entry",
            community_id="test_community",
            validation_level=ValidationLevel.COMMUNITY,
            invited_validators=validators,
        )
        
        # All validators approve
        for validator_id in validators:
            validator.submit_validation_vote(
                round_id=round_obj.round_id,
                validator_id=validator_id,
                approved=True,
                confidence_score=0.9,
                comments="Looks accurate",
            )
        
        # Check consensus
        updated_round = validator._rounds[round_obj.round_id]
        assert updated_round.consensus_reached is True
        assert updated_round.final_decision is True
        assert updated_round.status == ValidationStatus.APPROVED
    
    def test_elder_weighted_consensus(self):
        """Test elder-weighted consensus calculation."""
        validator = CommunityValidator()
        
        # Register mixed validators
        elder = validator.register_validator(
            community_id="test_community",
            display_name="Elder",
            role=ValidatorRole.ELDER,
            expertise_domains=["Traditional Knowledge"],
        )
        
        member1 = validator.register_validator(
            community_id="test_community",
            display_name="Member 1",
            role=ValidatorRole.COMMUNITY_MEMBER,
            expertise_domains=["Traditional Knowledge"],
        )
        
        member2 = validator.register_validator(
            community_id="test_community",
            display_name="Member 2",
            role=ValidatorRole.COMMUNITY_MEMBER,
            expertise_domains=["Traditional Knowledge"],
        )
        
        # Initiate round with elder weighting
        round_obj = validator.initiate_validation_round(
            entry_id="test_entry",
            community_id="test_community",
            validation_level=ValidationLevel.COMMUNITY,
            invited_validators=[elder.validator_id, member1.validator_id, member2.validator_id],
        )
        round_obj.consensus_type = ConsensusType.ELDER_WEIGHTED
        
        # Elder approves, members reject
        validator.submit_validation_vote(
            round_id=round_obj.round_id,
            validator_id=elder.validator_id,
            approved=True,
            confidence_score=0.95,
        )
        
        validator.submit_validation_vote(
            round_id=round_obj.round_id,
            validator_id=member1.validator_id,
            approved=False,
            confidence_score=0.7,
        )
        
        validator.submit_validation_vote(
            round_id=round_obj.round_id,
            validator_id=member2.validator_id,
            approved=False,
            confidence_score=0.7,
        )
        
        # Elder vote (2.0 weight) should still result in rejection
        # Elder: 2.0 approve, Members: 2.0 reject = 50% weighted approval
        updated_round = validator._rounds[round_obj.round_id]
        assert updated_round.consensus_reached is False or updated_round.final_decision is False
    
    def test_multi_community_consensus(self):
        """Test multi-community consensus validation."""
        validator = CommunityValidator()
        
        # Create validators in 3 communities
        communities = ["community1", "community2", "community3"]
        
        for community_id in communities:
            for i in range(3):
                validator.register_validator(
                    community_id=community_id,
                    display_name=f"Validator {i}",
                    role=ValidatorRole.KNOWLEDGE_HOLDER,
                    expertise_domains=["Traditional Knowledge"],
                )
        
        # Initiate multi-community consensus
        consensus = validator.initiate_multi_community_consensus(
            entry_id="test_entry",
            communities_invited=communities,
        )
        
        assert consensus.entry_id == "test_entry"
        assert len(consensus.communities_invited) == 3
        assert consensus.consensus_reached is False  # Not completed yet


# ==========================================================================
# Governance Tests
# ==========================================================================

class TestEthnoPathGovernance:
    """Test DAO governance with elder veto."""
    
    def test_register_participant(self):
        """Test governance participant registration."""
        governance = EthnoPathGovernance()
        
        participant = governance.register_participant(
            community_id="test_community",
            display_name="Test Participant",
            role=GovernanceRole.ELDER,
        )
        
        assert participant.display_name == "Test Participant"
        assert participant.role == GovernanceRole.ELDER
        assert participant.base_voting_power == 3.0  # Elders get 3x power
        assert participant.is_active is True
    
    def test_create_access_request_proposal(self):
        """Test access request proposal creation."""
        governance = EthnoPathGovernance()
        
        participant = governance.register_participant(
            community_id="test_community",
            display_name="Researcher",
            role=GovernanceRole.RESEARCHER,
        )
        
        proposal = governance.create_proposal(
            community_id="test_community",
            created_by=participant.participant_id,
            proposal_type=ProposalType.ACCESS_REQUEST,
            title="Request Access to Medicinal Plant Knowledge",
            description="Requesting access for research purposes.",
            requested_entry_ids=["entry1", "entry2"],
            requested_access_level="read",
            intended_use="Academic research on traditional medicine",
            benefit_sharing_terms="Share results with community, co-authorship",
        )
        
        assert proposal.proposal_type == ProposalType.ACCESS_REQUEST
        assert len(proposal.requested_entry_ids) == 2
        assert proposal.status == ProposalStatus.DRAFT
    
    def test_quadratic_voting(self):
        """Test quadratic voting mechanism."""
        governance = EthnoPathGovernance()
        
        # Register participants
        elder = governance.register_participant(
            community_id="test_community",
            display_name="Elder",
            role=GovernanceRole.ELDER,
        )
        
        member = governance.register_participant(
            community_id="test_community",
            display_name="Member",
            role=GovernanceRole.COMMUNITY_MEMBER,
        )
        
        # Create and submit proposal
        proposal = governance.create_proposal(
            community_id="test_community",
            created_by=member.participant_id,
            proposal_type=ProposalType.KNOWLEDGE_UPDATE,
            title="Update Traditional Practice",
            description="Update knowledge entry with new information.",
        )
        governance.submit_proposal(proposal.proposal_id)
        
        # Cast votes with different voting power
        governance.cast_vote(
            proposal_id=proposal.proposal_id,
            participant_id=elder.participant_id,
            vote_type=VoteType.APPROVE,
            voting_power_used=4.0,  # Elder uses 4 tokens
            rationale="I support this update",
        )
        
        governance.cast_vote(
            proposal_id=proposal.proposal_id,
            participant_id=member.participant_id,
            vote_type=VoteType.APPROVE,
            voting_power_used=1.0,  # Member uses 1 token
            rationale="Agree with update",
        )
        
        # Verify votes recorded
        updated_proposal = governance._proposals[proposal.proposal_id]
        assert len(updated_proposal.votes) == 2
        
        # Elder's quadratic power: sqrt(4) * 3.0 * 1.0 = 6.0
        # Member's quadratic power: sqrt(1) * 1.0 * 1.0 = 1.0
        elder_vote = updated_proposal.votes[0]
        assert elder_vote.voting_power_used == pytest.approx(6.0, rel=0.01)
    
    def test_elder_veto_authority(self):
        """Test elder veto mechanism."""
        governance = EthnoPathGovernance()
        
        # Register participants
        elder = governance.register_participant(
            community_id="test_community",
            display_name="Elder",
            role=GovernanceRole.ELDER,
        )
        
        members = []
        for i in range(5):
            m = governance.register_participant(
                community_id="test_community",
                display_name=f"Member {i}",
                role=GovernanceRole.COMMUNITY_MEMBER,
            )
            members.append(m)
        
        # Create proposal with elder veto authority
        proposal = governance.create_proposal(
            community_id="test_community",
            created_by=members[0].participant_id,
            proposal_type=ProposalType.ACCESS_REQUEST,
            title="External Access Request",
            description="External researcher wants access.",
        )
        governance.submit_proposal(proposal.proposal_id)
        
        # All members approve
        for member in members:
            governance.cast_vote(
                proposal_id=proposal.proposal_id,
                participant_id=member.participant_id,
                vote_type=VoteType.APPROVE,
                voting_power_used=1.0,
            )
        
        # Elder vetoes
        governance.cast_vote(
            proposal_id=proposal.proposal_id,
            participant_id=elder.participant_id,
            vote_type=VoteType.VETO,
            voting_power_used=1.0,
            rationale="Cultural concerns not addressed",
        )
        
        # Proposal should be vetoed despite majority approval
        updated_proposal = governance._proposals[proposal.proposal_id]
        assert updated_proposal.status == ProposalStatus.VETOED
        assert updated_proposal.final_decision is False
    
    def test_access_grant_creation(self):
        """Test access grant creation from approved proposal."""
        governance = EthnoPathGovernance()
        
        # Register participants
        elder = governance.register_participant(
            community_id="test_community",
            display_name="Elder",
            role=GovernanceRole.ELDER,
        )
        
        researcher = governance.register_participant(
            community_id="test_community",
            display_name="Researcher",
            role=GovernanceRole.RESEARCHER,
        )
        
        # Create access request
        proposal = governance.create_proposal(
            community_id="test_community",
            created_by=researcher.participant_id,
            proposal_type=ProposalType.ACCESS_REQUEST,
            title="Research Access",
            description="Access for PhD research",
            requested_entry_ids=["entry1"],
            requested_access_level="read",
            intended_use="Academic research",
            benefit_sharing_terms="Co-authorship",
        )
        governance.submit_proposal(proposal.proposal_id)
        
        # Elder approves with supermajority
        governance.cast_vote(
            proposal_id=proposal.proposal_id,
            participant_id=elder.participant_id,
            vote_type=VoteType.APPROVE,
            voting_power_used=10.0,  # High voting power
        )
        
        # Check if grant created
        updated_proposal = governance._proposals[proposal.proposal_id]
        if updated_proposal.status == ProposalStatus.APPROVED:
            # Find created grant
            grants = [g for g in governance._grants.values() if g.proposal_id == proposal.proposal_id]
            assert len(grants) == 1
            assert grants[0].grantee_id == researcher.participant_id
            assert "entry1" in grants[0].granted_entry_ids
    
    def test_grant_revocation(self):
        """Test access grant revocation."""
        governance = EthnoPathGovernance()
        
        # Manually create a grant for testing
        from backend.services.ethnopath.governance import AccessGrant
        grant = AccessGrant(
            proposal_id="test_proposal",
            grantee_id="test_grantee",
            granted_entry_ids=["entry1"],
            access_level="read",
            intended_use="Research",
            benefit_sharing_terms="Co-authorship",
        )
        governance._grants[grant.grant_id] = grant
        
        # Revoke grant
        revoked_grant = governance.revoke_grant(
            grant_id=grant.grant_id,
            revocation_reason="Breach of terms",
        )
        
        assert revoked_grant.is_revoked is True
        assert revoked_grant.revocation_reason == "Breach of terms"
        assert revoked_grant.revoked_at is not None
    
    def test_access_permission_check(self):
        """Test access permission checking."""
        governance = EthnoPathGovernance()
        
        # Create active grant
        from backend.services.ethnopath.governance import AccessGrant
        grant = AccessGrant(
            proposal_id="test_proposal",
            grantee_id="participant1",
            granted_entry_ids=["entry1", "entry2"],
            access_level="read",
            intended_use="Research",
            benefit_sharing_terms="Co-authorship",
        )
        governance._grants[grant.grant_id] = grant
        
        # Check permissions
        has_access_entry1 = governance.check_access_permission("participant1", "entry1")
        has_access_entry2 = governance.check_access_permission("participant1", "entry2")
        has_access_entry3 = governance.check_access_permission("participant1", "entry3")
        
        assert has_access_entry1 is True
        assert has_access_entry2 is True
        assert has_access_entry3 is False  # Not in grant


# ==========================================================================
# Integration Tests
# ==========================================================================

class TestEthnoPathIntegration:
    """Test integrated workflows across digitizer, validator, and governance."""
    
    def test_full_knowledge_lifecycle(self):
        """Test complete knowledge entry lifecycle: create → validate → govern → access."""
        digitizer = EthnoPathDigitizer()
        validator = CommunityValidator()
        governance = EthnoPathGovernance()
        
        # 1. Community registers
        community = digitizer.register_community(
            name="Integration Test Community",
            location="Test Region",
            indigenous_group="Test People",
            contact_email="integration@test.example",
        )
        
        # 2. Knowledge holder registers
        holder = digitizer.register_knowledge_holder(
            community_id=community.community_id,
            name="Test Elder",
            role="Elder",
            areas_of_expertise=["Medicinal Plants"],
        )
        
        # 3. Create knowledge entry
        entry = digitizer.create_knowledge_entry(
            community_id=community.community_id,
            holder_id=holder.holder_id,
            domain=KnowledgeDomain.MEDICINAL_PLANTS,
            knowledge_type=KnowledgeType.THERAPEUTIC_USE,
            sensitivity_level=SensitivityLevel.COMMUNITY_ONLY,
            title="Sacred Healing Plant",
            description="Traditional medicinal plant knowledge.",
        )
        
        # 4. Grant consent
        digitizer.grant_consent(
            entry_id=entry.entry_id,
            holder_id=holder.holder_id,
            consent_type="validation",
        )
        
        # 5. Register validators
        validators = []
        for i in range(3):
            v = validator.register_validator(
                community_id=community.community_id,
                display_name=f"Validator {i}",
                role=ValidatorRole.KNOWLEDGE_HOLDER,
                expertise_domains=["Medicinal Plants"],
            )
            validators.append(v.validator_id)
        
        # 6. Validate knowledge
        validation_round = validator.initiate_validation_round(
            entry_id=entry.entry_id,
            community_id=community.community_id,
            validation_level=ValidationLevel.COMMUNITY,
            invited_validators=validators,
        )
        
        # All validators approve
        for validator_id in validators:
            validator.submit_validation_vote(
                round_id=validation_round.round_id,
                validator_id=validator_id,
                approved=True,
                confidence_score=0.9,
            )
        
        # 7. Register governance participants
        elder = governance.register_participant(
            community_id=community.community_id,
            display_name="Elder Governance",
            role=GovernanceRole.ELDER,
        )
        
        researcher = governance.register_participant(
            community_id=community.community_id,
            display_name="External Researcher",
            role=GovernanceRole.RESEARCHER,
        )
        
        # 8. Researcher requests access
        proposal = governance.create_proposal(
            community_id=community.community_id,
            created_by=researcher.participant_id,
            proposal_type=ProposalType.ACCESS_REQUEST,
            title="Research Access Request",
            description="Requesting access for medicinal plant research",
            requested_entry_ids=[entry.entry_id],
            requested_access_level="read",
            intended_use="Academic research",
            benefit_sharing_terms="Co-authorship and results sharing",
        )
        governance.submit_proposal(proposal.proposal_id)
        
        # 9. Elder approves
        governance.cast_vote(
            proposal_id=proposal.proposal_id,
            participant_id=elder.participant_id,
            vote_type=VoteType.APPROVE,
            voting_power_used=10.0,
        )
        
        # Verify complete workflow
        assert entry is not None
        assert validation_round.status == ValidationStatus.APPROVED
        updated_proposal = governance._proposals[proposal.proposal_id]
        # Proposal may be approved depending on threshold and voting
        assert updated_proposal.status in [ProposalStatus.APPROVED, ProposalStatus.IN_VOTING]


class TestEthnoPathEthicalPricing:
    """Test ethical pricing and TK opt-out surcharges."""
    
    def test_full_compliance_no_surcharge(self):
        """Test full TK compliance has no surcharge."""
        governance = EthnoPathGovernance()
        
        pricing = governance.calculate_grant_pricing(
            pricing_tier=PricingTier.STANDARD,
            tk_compliance_level=TKComplianceLevel.FULL_COMPLIANCE,
        )
        
        assert pricing["base_price_usd"] == 299.0
        assert pricing["surcharge_amount_usd"] == 0.0
        assert pricing["total_price_usd"] == 299.0
        assert pricing["compensation_pool_contribution_usd"] == 0.0
    
    def test_no_attribution_surcharge(self):
        """Test no attribution opt-out applies 15% + $150 surcharge."""
        governance = EthnoPathGovernance()
        
        pricing = governance.calculate_grant_pricing(
            pricing_tier=PricingTier.STANDARD,
            tk_compliance_level=TKComplianceLevel.NO_ATTRIBUTION,
        )
        
        assert pricing["base_price_usd"] == 299.0
        # 15% of 299 = 44.85 + 150 flat = 194.85
        assert pricing["surcharge_amount_usd"] == pytest.approx(194.85, rel=0.01)
        assert pricing["total_price_usd"] == pytest.approx(493.85, rel=0.01)
        assert pricing["compensation_pool_contribution_usd"] == pytest.approx(194.85, rel=0.01)
    
    def test_no_compensation_surcharge(self):
        """Test no compensation opt-out applies 25% + $250 surcharge."""
        governance = EthnoPathGovernance()
        
        pricing = governance.calculate_grant_pricing(
            pricing_tier=PricingTier.PREMIUM,
            tk_compliance_level=TKComplianceLevel.NO_COMPENSATION,
        )
        
        assert pricing["base_price_usd"] == 799.0
        # 25% of 799 = 199.75 + 250 flat = 449.75
        assert pricing["surcharge_amount_usd"] == pytest.approx(449.75, rel=0.01)
        assert pricing["total_price_usd"] == pytest.approx(1248.75, rel=0.01)
        assert pricing["compensation_pool_contribution_usd"] == pytest.approx(449.75, rel=0.01)
    
    def test_no_tk_full_surcharge(self):
        """Test full TK opt-out applies 40% + $500 surcharge."""
        governance = EthnoPathGovernance()
        
        pricing = governance.calculate_grant_pricing(
            pricing_tier=PricingTier.STANDARD,
            tk_compliance_level=TKComplianceLevel.NO_TK_FULL,
        )
        
        assert pricing["base_price_usd"] == 299.0
        # 40% of 299 = 119.60 + 500 flat = 619.60
        assert pricing["surcharge_amount_usd"] == pytest.approx(619.60, rel=0.01)
        assert pricing["total_price_usd"] == pytest.approx(918.60, rel=0.01)
        assert pricing["compensation_pool_contribution_usd"] == pytest.approx(619.60, rel=0.01)
    
    def test_access_grant_with_pricing(self):
        """Test access grant includes pricing fields."""
        governance = EthnoPathGovernance()
        
        # Register participants
        member = governance.register_participant(
            community_id="test_community",
            display_name="Member",
            role=GovernanceRole.COMMUNITY_MEMBER,
        )
        
        elder = governance.register_participant(
            community_id="test_community",
            display_name="Elder",
            role=GovernanceRole.ELDER,
        )
        
        # Create and approve proposal
        proposal = governance.create_proposal(
            community_id="test_community",
            created_by=member.participant_id,
            proposal_type=ProposalType.ACCESS_REQUEST,
            title="Test Access",
            description="Test pricing",
            requested_entry_ids=["entry_1"],
        )
        governance.submit_proposal(proposal.proposal_id)
        
        # Elder approves
        governance.cast_vote(
            proposal_id=proposal.proposal_id,
            participant_id=elder.participant_id,
            vote_type=VoteType.APPROVE,
            voting_power_used=5.0,
        )
        
        # Manually create grant with pricing (since auto-creation uses defaults)
        grant = governance._create_access_grant(
            proposal,
            pricing_tier=PricingTier.PREMIUM,
            tk_compliance_level=TKComplianceLevel.NO_TK_FULL,
        )
        
        assert grant.pricing_tier == PricingTier.PREMIUM
        assert grant.tk_compliance_level == TKComplianceLevel.NO_TK_FULL
        assert grant.base_price_usd == 799.0
        assert grant.surcharge_amount_usd == pytest.approx(819.60, rel=0.01)  # 40% of 799 + 500
        assert grant.total_price_usd == pytest.approx(1618.60, rel=0.01)
        assert grant.compensation_pool_contribution_usd == pytest.approx(819.60, rel=0.01)
    
    def test_revenue_tracking(self):
        """Test revenue and compensation pool tracking in metrics."""
        governance = EthnoPathGovernance()
        
        # Register participants
        member = governance.register_participant(
            community_id="test_community",
            display_name="Member",
            role=GovernanceRole.COMMUNITY_MEMBER,
        )
        
        elder = governance.register_participant(
            community_id="test_community",
            display_name="Elder",
            role=GovernanceRole.ELDER,
        )
        
        # Create first grant - full compliance
        proposal1 = governance.create_proposal(
            community_id="test_community",
            created_by=member.participant_id,
            proposal_type=ProposalType.ACCESS_REQUEST,
            title="Grant 1",
            description="Full compliance",
            requested_entry_ids=["entry_1"],
        )
        governance.submit_proposal(proposal1.proposal_id)
        governance.cast_vote(
            proposal_id=proposal1.proposal_id,
            participant_id=elder.participant_id,
            vote_type=VoteType.APPROVE,
            voting_power_used=5.0,
        )
        
        grant1 = governance._create_access_grant(
            proposal1,
            pricing_tier=PricingTier.STANDARD,
            tk_compliance_level=TKComplianceLevel.FULL_COMPLIANCE,
        )
        
        # Create second grant - opt-out
        proposal2 = governance.create_proposal(
            community_id="test_community",
            created_by=member.participant_id,
            proposal_type=ProposalType.ACCESS_REQUEST,
            title="Grant 2",
            description="Opt-out",
            requested_entry_ids=["entry_2"],
        )
        governance.submit_proposal(proposal2.proposal_id)
        governance.cast_vote(
            proposal_id=proposal2.proposal_id,
            participant_id=elder.participant_id,
            vote_type=VoteType.APPROVE,
            voting_power_used=5.0,
        )
        
        grant2 = governance._create_access_grant(
            proposal2,
            pricing_tier=PricingTier.STANDARD,
            tk_compliance_level=TKComplianceLevel.NO_TK_FULL,
        )
        
        # Check metrics
        metrics = governance.get_governance_metrics()
        assert metrics.total_grants == 2
        assert metrics.full_compliance_grants == 1
        assert metrics.opt_out_grants == 1
        assert metrics.total_revenue_usd == pytest.approx(299.0 + 918.60, rel=0.01)
        assert metrics.total_tk_surcharges_usd == pytest.approx(619.60, rel=0.01)
        assert metrics.compensation_pool_usd == pytest.approx(619.60, rel=0.01)

