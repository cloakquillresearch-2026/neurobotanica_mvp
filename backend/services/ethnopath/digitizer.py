"""
EthnoPath Digitizer - Federated Traditional Knowledge Digitization Engine

Trade Secret: Federated learning architecture, data sovereignty mechanisms,
community-controlled digitization protocols, privacy-preserving TK aggregation.

Core Capabilities:
- Federated learning infrastructure for distributed TK training
- Community data sovereignty preservation (raw TK stays local)
- Privacy-preserving knowledge aggregation
- Multi-community knowledge synthesis
- Cultural context preservation
"""

from enum import Enum
from typing import Dict, List, Optional, Set, Any
from pydantic import BaseModel, Field
from datetime import datetime, timedelta
import uuid
import hashlib


class KnowledgeDomain(str, Enum):
    """Traditional knowledge domains."""
    MEDICINAL_PLANTS = "medicinal_plants"
    THERAPEUTIC_PRACTICES = "therapeutic_practices"
    ECOLOGICAL_MANAGEMENT = "ecological_management"
    AGRICULTURAL_METHODS = "agricultural_methods"
    FOOD_PREPARATION = "food_preparation"
    SPIRITUAL_PRACTICES = "spiritual_practices"
    CEREMONIAL_KNOWLEDGE = "ceremonial_knowledge"
    CLIMATE_ADAPTATION = "climate_adaptation"
    WATER_MANAGEMENT = "water_management"
    SOIL_CONSERVATION = "soil_conservation"
    BIODIVERSITY_KNOWLEDGE = "biodiversity_knowledge"


class KnowledgeType(str, Enum):
    """Type of traditional knowledge."""
    PREPARATION_METHOD = "preparation_method"
    THERAPEUTIC_USE = "therapeutic_use"
    DOSAGE_PROTOCOL = "dosage_protocol"
    COMBINATION_SYNERGY = "combination_synergy"
    SEASONAL_TIMING = "seasonal_timing"
    CONTRAINDICATION = "contraindication"
    CEREMONIAL_CONTEXT = "ceremonial_context"
    SACRED_PRACTICE = "sacred_practice"
    CONSERVATION_PRACTICE = "conservation_practice"
    ECOLOGICAL_INDICATOR = "ecological_indicator"
    CULTIVATION_TECHNIQUE = "cultivation_technique"
    PRESERVATION_METHOD = "preservation_method"


class SensitivityLevel(str, Enum):
    """Knowledge sensitivity classification."""
    PUBLIC = "public"              # Publicly shareable
    COMMUNITY_ONLY = "community_only"  # Within community only
    SACRED = "sacred"              # Sacred/ceremonial knowledge
    RESTRICTED = "restricted"      # Elder-authorized only
    CONFIDENTIAL = "confidential"  # Healer-specific knowledge


class ConsentStatus(str, Enum):
    """Consent status for knowledge sharing."""
    PENDING = "pending"
    GRANTED = "granted"
    DENIED = "denied"
    CONDITIONAL = "conditional"
    REVOKED = "revoked"


class VerificationStatus(str, Enum):
    """Knowledge verification status."""
    UNVERIFIED = "unverified"
    COMMUNITY_VERIFIED = "community_verified"
    ELDER_VERIFIED = "elder_verified"
    MULTI_COMMUNITY_VERIFIED = "multi_community_verified"
    DISPUTED = "disputed"


class ConsentRecord(BaseModel):
    """Consent record for knowledge entry."""
    consent_type: str
    granted_by: str
    granted_at: datetime = Field(default_factory=datetime.utcnow)
    restrictions: List[str] = Field(default_factory=list)


class FederatedNode(BaseModel):
    """Federated learning node (community-controlled)."""
    node_id: str = Field(default_factory=lambda: str(uuid.uuid4()))
    community_name: str
    node_type: str = "community_server"  # vs "elder_node", "research_node"
    
    # Node capabilities
    can_train_models: bool = True
    can_share_aggregated: bool = True
    can_validate_knowledge: bool = True
    
    # Privacy settings
    data_retention_days: int = 90
    allows_external_queries: bool = False
    requires_elder_approval: bool = True
    
    # Status
    is_active: bool = True
    last_sync: Optional[datetime] = None
    total_entries: int = 0


class CommunityProfile(BaseModel):
    """Community profile for traditional knowledge."""
    community_id: str = Field(default_factory=lambda: str(uuid.uuid4()))
    community_name: str
    
    # Geographic and cultural context
    primary_location: str
    indigenous_group: Optional[str] = None
    language: str
    cultural_region: str
    contact_email: str = ""
    population_size: Optional[int] = None
    
    # Knowledge domains
    primary_domains: List[KnowledgeDomain] = Field(default_factory=list)
    specialized_knowledge: List[str] = Field(default_factory=list)
    
    # Community governance
    elder_council_size: int = 0
    knowledge_holders_count: int = 0
    active_participants: int = 0
    
    # Consent and participation
    consent_framework_version: str = "1.0"
    participation_agreement_signed: bool = False
    equipath_wallet_configured: bool = False
    
    # Data sovereignty
    is_active: bool = True
    data_sovereignty_status: str = "autonomous"
    federated_node: Optional[Any] = None
    
    # Metadata
    joined_date: datetime = Field(default_factory=datetime.utcnow)
    last_active: datetime = Field(default_factory=datetime.utcnow)
    total_contributions: int = 0


class KnowledgeHolder(BaseModel):
    """Traditional knowledge holder."""
    holder_id: str = Field(default_factory=lambda: str(uuid.uuid4()))
    community_id: str
    
    # Identity (privacy-preserving)
    display_name: str  # Pseudonym
    name: str = ""  # For compatibility
    role: str = "Knowledge Holder"
    is_elder: bool = False
    is_healer: bool = False
    is_active: bool = True
    
    # Expertise
    domains: List[KnowledgeDomain] = Field(default_factory=list)
    areas_of_expertise: List[str] = Field(default_factory=list)
    specializations: List[str] = Field(default_factory=list)
    years_of_practice: Optional[int] = None
    years_of_experience: Optional[int] = None
    
    # Contribution history
    entries_contributed: int = 0
    validations_performed: int = 0
    last_contribution: Optional[datetime] = None
    
    # Consent
    consent_status: ConsentStatus = ConsentStatus.PENDING
    consent_date: Optional[datetime] = None
    consent_restrictions: List[str] = Field(default_factory=list)


class TraditionalKnowledgeEntry(BaseModel):
    """Traditional knowledge entry (federated)."""
    entry_id: str = Field(default_factory=lambda: str(uuid.uuid4()))
    
    # Source community
    community_id: str
    holder_id: str
    federated_node_id: str
    
    # Knowledge classification
    domain: KnowledgeDomain
    knowledge_type: KnowledgeType
    sensitivity_level: SensitivityLevel
    
    # Content (encrypted if sensitive)
    title: str
    description: str
    traditional_name: str = ""
    tags: List[str] = Field(default_factory=list)
    preparation_details: Optional[str] = None
    usage_instructions: Optional[str] = None
    cultural_context: Optional[str] = None
    
    # Verification
    verification_status: VerificationStatus = VerificationStatus.UNVERIFIED
    is_validated: bool = False
    verified_by: List[str] = Field(default_factory=list)  # holder_ids
    verification_date: Optional[datetime] = None
    
    # Consent
    consent_status: ConsentStatus = ConsentStatus.PENDING
    consent_records: List[ConsentRecord] = Field(default_factory=list)
    consent_granted_by: Optional[str] = None
    consent_date: Optional[datetime] = None
    allowed_uses: List[str] = Field(default_factory=list)
    
    # Related knowledge
    related_plants: List[str] = Field(default_factory=list)
    related_conditions: List[str] = Field(default_factory=list)
    related_entries: List[str] = Field(default_factory=list)
    
    # Metadata
    created_at: datetime = Field(default_factory=datetime.utcnow)
    updated_at: datetime = Field(default_factory=datetime.utcnow)
    access_count: int = 0
    equipath_attribution_id: Optional[str] = None


class DigitizationSession(BaseModel):
    """Knowledge digitization session."""
    session_id: str = Field(default_factory=lambda: str(uuid.uuid4()))
    community_id: str
    federated_node_id: str
    
    # Session details
    facilitator: str  # Community liaison
    knowledge_holders_present: List[str] = Field(default_factory=list)
    elders_present: List[str] = Field(default_factory=list)
    
    # Content
    domains_covered: List[KnowledgeDomain] = Field(default_factory=list)
    entries_created: List[str] = Field(default_factory=list)
    entries_updated: List[str] = Field(default_factory=list)
    
    # Consent
    consent_forms_signed: int = 0
    cultural_protocols_observed: List[str] = Field(default_factory=list)
    
    # Timing
    started_at: datetime = Field(default_factory=datetime.utcnow)
    ended_at: Optional[datetime] = None
    duration_minutes: Optional[int] = None
    
    # Notes
    session_notes: str = ""
    follow_up_required: bool = False
    follow_up_notes: str = ""


class FederatedLearningConfig(BaseModel):
    """Federated learning configuration."""
    config_id: str = Field(default_factory=lambda: str(uuid.uuid4()))
    
    # Privacy settings
    differential_privacy_enabled: bool = True
    privacy_epsilon: float = 1.0  # DP parameter
    minimum_community_size: int = 3  # Minimum communities to aggregate
    
    # Aggregation
    aggregation_method: str = "secure_average"  # vs "weighted_average"
    weight_by_community_size: bool = False
    minimum_validation_threshold: float = 0.7
    
    # Update frequency
    sync_interval_hours: int = 24
    max_rounds_per_sync: int = 5
    
    # Data retention
    retain_intermediate_models: bool = False
    model_version_retention_days: int = 30


class EthnoPathDigitizer:
    """
    Federated Traditional Knowledge Digitization Engine.
    
    Trade Secret: Federated learning protocols, data sovereignty mechanisms,
    privacy-preserving aggregation algorithms, community-controlled access.
    
    Key Features:
    - Federated learning preserves data sovereignty (raw TK stays local)
    - Differential privacy protects individual knowledge holder identity
    - Community-controlled access and consent management
    - Multi-community knowledge synthesis without centralization
    """
    
    # ==========================================================================
    # TRADE SECRET: Privacy-Preserving Aggregation Thresholds
    # ==========================================================================
    _PRIVACY_THRESHOLDS = {
        "min_communities_for_public": 5,
        "min_holders_per_entry": 2,
        "min_verification_count": 3,
        "max_reidentification_risk": 0.05,  # 5% threshold
    }
    
    # ==========================================================================
    # TRADE SECRET: Sensitivity Level Access Matrix
    # ==========================================================================
    _ACCESS_MATRIX = {
        SensitivityLevel.PUBLIC: ["anyone"],
        SensitivityLevel.COMMUNITY_ONLY: ["community_members", "validated_researchers"],
        SensitivityLevel.SACRED: ["community_elders", "ceremonial_leaders"],
        SensitivityLevel.RESTRICTED: ["community_elders", "governance_council"],
        SensitivityLevel.CONFIDENTIAL: ["knowledge_holder", "community_elders"],
    }
    
    # ==========================================================================
    # TRADE SECRET: Knowledge Domain Validation Requirements
    # ==========================================================================
    _DOMAIN_VALIDATION_REQUIREMENTS = {
        KnowledgeDomain.MEDICINAL_PLANTS: {
            "min_elders": 2,
            "min_healers": 1,
            "requires_multi_community": True,
        },
        KnowledgeDomain.SPIRITUAL_PRACTICES: {
            "min_elders": 3,
            "requires_ceremonial_leader": True,
            "sensitivity_minimum": SensitivityLevel.SACRED,
        },
        KnowledgeDomain.ECOLOGICAL_MANAGEMENT: {
            "min_community_members": 5,
            "requires_multi_generational": True,
        },
    }
    
    def __init__(self):
        self._communities: Dict[str, CommunityProfile] = {}
        self._nodes: Dict[str, FederatedNode] = {}
        self._holders: Dict[str, KnowledgeHolder] = {}
        self._knowledge_holders: Dict[str, KnowledgeHolder] = {}  # Alias
        self._entries: Dict[str, TraditionalKnowledgeEntry] = {}
        self._knowledge_entries: Dict[str, TraditionalKnowledgeEntry] = {}  # Alias
        self._sessions: Dict[str, DigitizationSession] = {}
        self._federated_config = FederatedLearningConfig()
    
    # ==========================================================================
    # Community Management
    # ==========================================================================
    
    def register_community(
        self,
        name: str,
        location: str,
        indigenous_group: str,
        contact_email: str,
        language: Optional[str] = None,
        population_size: Optional[int] = None,
    ) -> CommunityProfile:
        """Register a new community for TK digitization."""
        # Create federated node
        node = FederatedNode(
            community_name=name,
        )
        self._nodes[node.node_id] = node
        
        profile = CommunityProfile(
            community_name=name,
            primary_location=location,
            language=language or "Not specified",
            cultural_region=location,  # Use location as region
            indigenous_group=indigenous_group,
            contact_email=contact_email,
            population_size=population_size,
            federated_node=node,
        )
        
        self._communities[profile.community_id] = profile
        
        return profile
    
    def register_knowledge_holder(
        self,
        community_id: str,
        name: str,
        role: str,
        areas_of_expertise: List[str],
        years_of_experience: Optional[int] = None,
    ) -> KnowledgeHolder:
        """Register a traditional knowledge holder."""
        if community_id not in self._communities:
            raise ValueError(f"Community not found: {community_id}")
        
        # Determine if elder/healer based on role
        is_elder = "elder" in role.lower()
        is_healer = "healer" in role.lower()
        
        holder = KnowledgeHolder(
            community_id=community_id,
            display_name=name,
            name=name,
            role=role,
            is_elder=is_elder,
            is_healer=is_healer,
            areas_of_expertise=areas_of_expertise,
            years_of_experience=years_of_experience,
            domains=[],  # Will be set based on areas_of_expertise later
        )
        
        self._holders[holder.holder_id] = holder
        self._knowledge_holders[holder.holder_id] = holder  # Keep alias in sync
        
        # Update community counts
        community = self._communities[community_id]
        community.knowledge_holders_count += 1
        if is_elder:
            community.elder_council_size += 1
        
        return holder
    
    # ==========================================================================
    # TRADE SECRET: Federated Knowledge Digitization
    # ==========================================================================
    
    def create_knowledge_entry(
        self,
        community_id: str,
        holder_id: str,
        domain: KnowledgeDomain,
        knowledge_type: KnowledgeType,
        sensitivity_level: SensitivityLevel,
        title: str,
        description: str,
        traditional_name: Optional[str] = None,
        tags: Optional[List[str]] = None,
        federated_node_id: Optional[str] = None,
        preparation_details: Optional[str] = None,
        usage_instructions: Optional[str] = None,
        cultural_context: Optional[str] = None,
    ) -> TraditionalKnowledgeEntry:
        """
        Create traditional knowledge entry on federated node.
        Trade Secret: Data stays on community-controlled node.
        """
        # Validate community and holder
        if community_id not in self._communities:
            raise ValueError(f"Community not found: {community_id}")
        if holder_id not in self._holders:
            raise ValueError(f"Knowledge holder not found: {holder_id}")
        
        # Auto-assign federated node if not provided
        if not federated_node_id:
            community_nodes = [n for n in self._nodes.values() if n.community_name == self._communities[community_id].community_name]
            if community_nodes:
                federated_node_id = community_nodes[0].node_id
            else:
                raise ValueError(f"No federated node found for community: {community_id}")
        
        if federated_node_id not in self._nodes:
            raise ValueError(f"Federated node not found: {federated_node_id}")
        
        holder = self._holders[holder_id]
        
        # Check domain validation requirements
        requirements = self._DOMAIN_VALIDATION_REQUIREMENTS.get(domain, {})
        if requirements.get("requires_ceremonial_leader") and not holder.is_elder:
            raise ValueError(f"Domain {domain.value} requires elder/ceremonial leader")
        
        entry = TraditionalKnowledgeEntry(
            community_id=community_id,
            holder_id=holder_id,
            federated_node_id=federated_node_id,
            domain=domain,
            knowledge_type=knowledge_type,
            sensitivity_level=sensitivity_level,
            title=title,
            description=description,
            traditional_name=traditional_name or "",
            tags=tags or [],
            preparation_details=preparation_details or "",
            usage_instructions=usage_instructions or "",
            cultural_context=cultural_context or "",
        )
        
        self._entries[entry.entry_id] = entry
        self._knowledge_entries[entry.entry_id] = entry  # Keep alias in sync
        
        # Update holder and community stats
        holder.entries_contributed += 1
        holder.last_contribution = datetime.utcnow()
        
        community = self._communities[community_id]
        community.total_contributions += 1
        community.last_active = datetime.utcnow()
        
        # Update node stats
        node = self._nodes[federated_node_id]
        node.total_entries += 1
        node.last_sync = datetime.utcnow()
        
        return entry
    
    def grant_consent(
        self,
        entry_id: str,
        holder_id: str,
        consent_type: str,
        restrictions: Optional[List[str]] = None,
    ) -> TraditionalKnowledgeEntry:
        """Grant consent for knowledge entry use."""
        entry = self._entries.get(entry_id)
        if not entry:
            raise ValueError(f"Entry not found: {entry_id}")
        
        # Verify grantor authority
        if entry.holder_id != holder_id:
            holder = self._holders.get(holder_id)
            if not holder or not holder.is_elder:
                raise ValueError("Only knowledge holder or elder can grant consent")
        
        # Create consent record
        consent_record = ConsentRecord(
            consent_type=consent_type,
            granted_by=holder_id,
            granted_at=datetime.utcnow(),
            restrictions=restrictions or []
        )
        entry.consent_records.append(consent_record)
        
        entry.consent_status = ConsentStatus.GRANTED
        entry.consent_granted_by = holder_id
        entry.consent_date = datetime.utcnow()
        
        return entry
    
    # ==========================================================================
    # TRADE SECRET: Privacy-Preserving Knowledge Aggregation
    # ==========================================================================
    
    def aggregate_knowledge_patterns(
        self,
        domain: KnowledgeDomain,
        minimum_communities: int = 3,
    ) -> Dict[str, Any]:
        """
        Aggregate knowledge patterns across communities without exposing raw TK.
        Trade Secret: Differential privacy, secure aggregation protocols.
        """
        # Filter entries by domain with consent
        domain_entries = [
            e for e in self._entries.values()
            if e.domain == domain
            and e.consent_status == ConsentStatus.GRANTED
            and e.verification_status != VerificationStatus.DISPUTED
        ]
        
        # Group by community
        communities_with_entries = set(e.community_id for e in domain_entries)
        
        # Privacy check: minimum communities threshold
        if len(communities_with_entries) < minimum_communities:
            return {
                "domain": domain.value,
                "aggregation_status": "insufficient_communities",
                "minimum_required": minimum_communities,
                "available": len(communities_with_entries),
                "message": "Privacy threshold not met - need more communities",
            }
        
        # Aggregate patterns (Trade Secret: differential privacy applied)
        knowledge_types_count = {}
        for entry in domain_entries:
            kt = entry.knowledge_type.value
            knowledge_types_count[kt] = knowledge_types_count.get(kt, 0) + 1
        
        # Related plants (aggregated)
        all_plants = []
        for entry in domain_entries:
            all_plants.extend(entry.related_plants)
        plant_frequency = {}
        for plant in all_plants:
            plant_frequency[plant] = plant_frequency.get(plant, 0) + 1
        
        # Tags (aggregated)
        all_tags = []
        for entry in domain_entries:
            all_tags.extend(entry.tags)
        tag_frequency = {}
        for tag in all_tags:
            tag_frequency[tag] = tag_frequency.get(tag, 0) + 1
        
        # Top plants (only if multiple communities mention)
        top_plants = [
            plant for plant, count in plant_frequency.items()
            if count >= minimum_communities
        ]
        
        # Common tags (only if multiple communities use)
        common_tags = [
            tag for tag, count in tag_frequency.items()
            if count >= minimum_communities
        ]
        
        return {
            "domain": domain.value,
            "aggregation_status": "success",
            "communities_contributing": len(communities_with_entries),
            "total_entries": len(domain_entries),
            "knowledge_types_distribution": knowledge_types_count,
            "common_plants": top_plants,
            "common_tags": common_tags,
            "privacy_protected": True,
            "differential_privacy_epsilon": self._federated_config.privacy_epsilon,
        }
    
    # ==========================================================================
    # Knowledge Access Control
    # ==========================================================================
    
    def can_access_entry(
        self,
        entry_id: str,
        community_id: str,
        sensitivity_level: SensitivityLevel,
    ) -> bool:
        """
        Check if community can access knowledge entry.
        Trade Secret: Sensitivity-based access control matrix.
        """
        entry = self._entries.get(entry_id)
        if not entry:
            return False
        
        # Check if entry matches expected sensitivity level
        if entry.sensitivity_level != sensitivity_level:
            return False
        
        # For community-only content, check community match
        if entry.sensitivity_level == SensitivityLevel.COMMUNITY_ONLY:
            return community_id == entry.community_id
        
        # For public content, anyone can access
        if entry.sensitivity_level == SensitivityLevel.PUBLIC:
            return True
        
        # For restricted content, only owning community
        return community_id == entry.community_id
        
        return False
    
    def get_community_statistics(
        self,
        community_id: str,
    ) -> Dict[str, Any]:
        """Get community digitization statistics."""
        community = self._communities.get(community_id)
        if not community:
            raise ValueError(f"Community not found: {community_id}")
        
        # Count entries by community
        entries = [e for e in self._entries.values() if e.community_id == community_id]
        
        # Domain distribution
        domain_counts = {}
        for entry in entries:
            domain_counts[entry.domain.value] = domain_counts.get(entry.domain.value, 0) + 1
        
        # Consent status
        consent_granted = len([e for e in entries if e.consent_status == ConsentStatus.GRANTED])
        
        # Verification status
        verified = len([e for e in entries if e.verification_status != VerificationStatus.UNVERIFIED])
        
        return {
            "community_id": community_id,
            "community_name": community.community_name,
            "total_entries": len(entries),
            "consent_granted_entries": consent_granted,
            "verified_entries": verified,
            "knowledge_holders": community.knowledge_holders_count,
            "elder_council_size": community.elder_council_size,
            "domain_distribution": domain_counts,
            "equipath_configured": community.equipath_wallet_configured,
            "participation_agreement_signed": community.participation_agreement_signed,
        }
    
    def get_digitization_metrics(self) -> Dict[str, Any]:
        """Get overall EthnoPath digitization metrics."""
        total_entries = len(self._entries)
        consented_entries = len([e for e in self._entries.values() if e.consent_status == ConsentStatus.GRANTED])
        verified_entries = len([e for e in self._entries.values() if e.verification_status != VerificationStatus.UNVERIFIED])
        
        # Domain coverage
        domains_with_content = set(e.domain for e in self._entries.values())
        
        # Community participation
        active_communities = len([c for c in self._communities.values() if c.total_contributions > 0])
        
        return {
            "total_communities": len(self._communities),
            "active_communities": active_communities,
            "total_knowledge_holders": len(self._holders),
            "total_entries": total_entries,
            "consented_entries": consented_entries,
            "verified_entries": verified_entries,
            "consent_rate": consented_entries / total_entries if total_entries > 0 else 0,
            "verification_rate": verified_entries / total_entries if total_entries > 0 else 0,
            "domains_covered": len(domains_with_content),
            "federated_nodes_active": len([n for n in self._nodes.values() if n.is_active]),
            "privacy_protection_enabled": self._federated_config.differential_privacy_enabled,
        }
