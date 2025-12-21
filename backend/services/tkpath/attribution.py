"""
TKPath Attribution Engine

Determines which indigenous communities contributed to formulations and
calculates attribution percentages based on knowledge elements used.

Trade Secret: Community weighting algorithms, contribution scoring formulas,
traditional knowledge element classification system.

Nagoya Protocol Compliance: Implements prior informed consent tracking and
fair/equitable benefit-sharing calculations per CBD Article 8(j).
"""

from enum import Enum
from typing import Dict, List, Optional, Set
from pydantic import BaseModel, Field
from datetime import datetime
import uuid
import hashlib


class KnowledgeElementType(str, Enum):
    """Types of traditional knowledge contributions."""
    PREPARATION_METHOD = "preparation_method"      # Traditional extraction/preparation
    TIMING_KNOWLEDGE = "timing_knowledge"          # Harvest timing (lunar, seasonal)
    SYNERGY_KNOWLEDGE = "synergy_knowledge"        # Compound combination wisdom
    DOSING_PROTOCOL = "dosing_protocol"            # Traditional dosage guidance
    CONDITION_APPLICATION = "condition_application"  # Historical medical use
    CULTIVATION_METHOD = "cultivation_method"      # Traditional growing practices
    STORAGE_PRESERVATION = "storage_preservation"  # Traditional preservation methods
    CEREMONIAL_CONTEXT = "ceremonial_context"      # Spiritual/ceremonial significance


class KnowledgeElement(BaseModel):
    """Individual element of traditional knowledge."""
    element_id: str = Field(default_factory=lambda: str(uuid.uuid4()))
    element_type: KnowledgeElementType
    description: str
    source_community: str
    documentation_level: str  # "oral_tradition", "documented", "verified"
    contribution_weight: float = Field(ge=0.0, le=1.0)  # 0-1 weight factor
    verification_status: str = "pending"  # "pending", "verified", "disputed"
    nagoya_consent_ref: Optional[str] = None  # Reference to consent documentation
    
    class Config:
        use_enum_values = True


class CommunityAttribution(BaseModel):
    """Attribution record for a specific indigenous community."""
    community_id: str = Field(default_factory=lambda: str(uuid.uuid4()))
    community_name: str
    region: str
    country: str
    knowledge_elements: List[KnowledgeElement] = Field(default_factory=list)
    total_contribution_percentage: float = 0.0
    compensation_address: Optional[str] = None  # Blockchain/payment address
    nagoya_pico_reference: Optional[str] = None  # Prior Informed Consent ref
    undrip_acknowledgment: bool = False  # UN Declaration compliance
    
    def calculate_contribution(self) -> float:
        """Calculate total contribution weight from all knowledge elements."""
        if not self.knowledge_elements:
            return 0.0
        total = sum(elem.contribution_weight for elem in self.knowledge_elements)
        return min(total, 1.0)  # Cap at 100%


class TKContribution(BaseModel):
    """Complete TK contribution analysis for a formulation."""
    contribution_id: str = Field(default_factory=lambda: str(uuid.uuid4()))
    formulation_id: str
    formulation_name: str
    analysis_timestamp: datetime = Field(default_factory=datetime.utcnow)
    
    # Community attributions
    communities: List[CommunityAttribution] = Field(default_factory=list)
    
    # Summary metrics
    total_tk_contribution_percentage: float = 0.0
    primary_community: Optional[str] = None
    knowledge_element_count: int = 0
    
    # Compensation suggestions
    suggested_compensation_percentage: float = 0.0  # % of revenue to TK holders
    suggested_minimum_payment: float = 0.0  # Minimum monthly payment
    
    # Compliance
    nagoya_compliant: bool = False
    undrip_compliant: bool = False
    verification_hash: Optional[str] = None  # Immutable record hash
    
    def generate_verification_hash(self) -> str:
        """Generate tamper-evident hash of attribution data."""
        data = f"{self.formulation_id}:{self.analysis_timestamp.isoformat()}"
        for comm in self.communities:
            data += f":{comm.community_name}:{comm.total_contribution_percentage}"
        self.verification_hash = hashlib.sha256(data.encode()).hexdigest()
        return self.verification_hash


class TKPathAttribution:
    """
    Traditional Knowledge attribution and compensation engine.
    
    Trade Secret: Community weighting algorithms, contribution scoring,
    knowledge element classification, compensation distribution formulas.
    
    Key Innovation: Economic rationality aligned with ethics - TK features
    are FREE, opting out costs MORE, creating sustainable compensation model.
    """
    
    # Trade Secret: Knowledge element weights by type
    _ELEMENT_TYPE_WEIGHTS: Dict[KnowledgeElementType, float] = {
        KnowledgeElementType.PREPARATION_METHOD: 0.25,
        KnowledgeElementType.SYNERGY_KNOWLEDGE: 0.20,
        KnowledgeElementType.CONDITION_APPLICATION: 0.15,
        KnowledgeElementType.DOSING_PROTOCOL: 0.15,
        KnowledgeElementType.TIMING_KNOWLEDGE: 0.10,
        KnowledgeElementType.CULTIVATION_METHOD: 0.08,
        KnowledgeElementType.STORAGE_PRESERVATION: 0.04,
        KnowledgeElementType.CEREMONIAL_CONTEXT: 0.03,
    }
    
    # Trade Secret: Compensation tier thresholds
    _COMPENSATION_TIERS = {
        "minimal": {"threshold": 0.10, "percentage": 2.0, "minimum": 50.0},
        "moderate": {"threshold": 0.30, "percentage": 5.0, "minimum": 150.0},
        "significant": {"threshold": 0.50, "percentage": 8.0, "minimum": 300.0},
        "primary": {"threshold": 0.70, "percentage": 12.0, "minimum": 500.0},
        "foundational": {"threshold": 0.90, "percentage": 18.0, "minimum": 1000.0},
    }
    
    # Known indigenous communities with cannabis TK
    _KNOWN_COMMUNITIES: Dict[str, Dict] = {
        "hindu_kush_region": {
            "name": "Hindu Kush Traditional Cultivators",
            "region": "Hindu Kush Mountains",
            "country": "Afghanistan/Pakistan",
            "primary_knowledge": [
                KnowledgeElementType.CULTIVATION_METHOD,
                KnowledgeElementType.PREPARATION_METHOD,
            ],
        },
        "jamaican_rastafari": {
            "name": "Jamaican Rastafari Communities",
            "region": "Jamaica",
            "country": "Jamaica",
            "primary_knowledge": [
                KnowledgeElementType.CEREMONIAL_CONTEXT,
                KnowledgeElementType.PREPARATION_METHOD,
                KnowledgeElementType.DOSING_PROTOCOL,
            ],
        },
        "moroccan_rif": {
            "name": "Moroccan Rif Mountain Cultivators",
            "region": "Rif Mountains",
            "country": "Morocco",
            "primary_knowledge": [
                KnowledgeElementType.CULTIVATION_METHOD,
                KnowledgeElementType.PREPARATION_METHOD,
            ],
        },
        "indian_ayurvedic": {
            "name": "Indian Ayurvedic Practitioners",
            "region": "Northern India",
            "country": "India",
            "primary_knowledge": [
                KnowledgeElementType.CONDITION_APPLICATION,
                KnowledgeElementType.DOSING_PROTOCOL,
                KnowledgeElementType.SYNERGY_KNOWLEDGE,
            ],
        },
        "chinese_tcm": {
            "name": "Chinese Traditional Medicine Practitioners",
            "region": "Various",
            "country": "China",
            "primary_knowledge": [
                KnowledgeElementType.CONDITION_APPLICATION,
                KnowledgeElementType.SYNERGY_KNOWLEDGE,
                KnowledgeElementType.TIMING_KNOWLEDGE,
            ],
        },
        "mexican_indigenous": {
            "name": "Mexican Indigenous Communities",
            "region": "Oaxaca/Guerrero",
            "country": "Mexico",
            "primary_knowledge": [
                KnowledgeElementType.CULTIVATION_METHOD,
                KnowledgeElementType.CEREMONIAL_CONTEXT,
            ],
        },
        "thai_traditional": {
            "name": "Thai Traditional Medicine Practitioners",
            "region": "Northern Thailand",
            "country": "Thailand",
            "primary_knowledge": [
                KnowledgeElementType.CONDITION_APPLICATION,
                KnowledgeElementType.PREPARATION_METHOD,
            ],
        },
        "south_african_dagga": {
            "name": "South African Dagga Cultivators",
            "region": "Eastern Cape",
            "country": "South Africa",
            "primary_knowledge": [
                KnowledgeElementType.CULTIVATION_METHOD,
                KnowledgeElementType.CEREMONIAL_CONTEXT,
            ],
        },
    }
    
    def __init__(self):
        self._attribution_cache: Dict[str, TKContribution] = {}
    
    def calculate_tk_contribution(
        self,
        formulation: str,
        formulation_id: Optional[str] = None,
        compounds: Optional[List[str]] = None,
        preparation_method: Optional[str] = None,
        claimed_communities: Optional[List[str]] = None,
    ) -> TKContribution:
        """
        Determine which indigenous communities contributed to formulation.
        
        Returns:
        - Community attribution percentages
        - Specific knowledge elements used (preparation, timing, synergies)
        - Suggested compensation amounts
        
        Trade Secret: Attribution scoring algorithm, community matching logic.
        """
        contribution = TKContribution(
            formulation_id=formulation_id or str(uuid.uuid4()),
            formulation_name=formulation,
        )
        
        # Analyze formulation for TK elements
        detected_elements = self._detect_knowledge_elements(
            formulation, compounds, preparation_method
        )
        
        # Match elements to communities
        community_matches = self._match_communities(
            detected_elements, claimed_communities
        )
        
        # Build attribution records
        for community_id, elements in community_matches.items():
            comm_info = self._KNOWN_COMMUNITIES.get(community_id, {})
            
            community_attribution = CommunityAttribution(
                community_name=comm_info.get("name", community_id),
                region=comm_info.get("region", "Unknown"),
                country=comm_info.get("country", "Unknown"),
                knowledge_elements=elements,
            )
            
            # Calculate contribution percentage
            community_attribution.total_contribution_percentage = (
                community_attribution.calculate_contribution() * 100
            )
            
            contribution.communities.append(community_attribution)
        
        # Calculate summary metrics
        contribution.total_tk_contribution_percentage = sum(
            c.total_contribution_percentage for c in contribution.communities
        )
        contribution.knowledge_element_count = sum(
            len(c.knowledge_elements) for c in contribution.communities
        )
        
        # Identify primary community
        if contribution.communities:
            primary = max(
                contribution.communities,
                key=lambda c: c.total_contribution_percentage
            )
            contribution.primary_community = primary.community_name
        
        # Calculate compensation suggestions
        self._calculate_compensation_suggestions(contribution)
        
        # Check compliance
        contribution.nagoya_compliant = self._check_nagoya_compliance(contribution)
        contribution.undrip_compliant = self._check_undrip_compliance(contribution)
        
        # Generate verification hash
        contribution.generate_verification_hash()
        
        # Cache for audit
        self._attribution_cache[contribution.contribution_id] = contribution
        
        return contribution
    
    def _detect_knowledge_elements(
        self,
        formulation: str,
        compounds: Optional[List[str]],
        preparation_method: Optional[str],
    ) -> List[KnowledgeElement]:
        """
        Detect TK elements present in formulation.
        Trade Secret: Detection algorithms and keyword matching.
        """
        elements = []
        formulation_lower = formulation.lower()
        
        # Detection patterns (simplified - real implementation more complex)
        detection_patterns = {
            KnowledgeElementType.PREPARATION_METHOD: [
                "traditional", "extract", "infusion", "decoction", "tincture",
                "cold-press", "sun-dried", "fermented", "charras", "hashish"
            ],
            KnowledgeElementType.SYNERGY_KNOWLEDGE: [
                "entourage", "synergy", "full-spectrum", "whole-plant",
                "combined", "ratio", "blend"
            ],
            KnowledgeElementType.CONDITION_APPLICATION: [
                "pain", "anxiety", "sleep", "inflammation", "nausea",
                "appetite", "seizure", "epilepsy", "ptsd", "cancer"
            ],
            KnowledgeElementType.DOSING_PROTOCOL: [
                "microdose", "titration", "sublingual", "ceremonial dose",
                "therapeutic dose", "maintenance"
            ],
            KnowledgeElementType.TIMING_KNOWLEDGE: [
                "harvest moon", "morning harvest", "sunset", "seasonal",
                "flowering time", "trichome maturity"
            ],
            KnowledgeElementType.CULTIVATION_METHOD: [
                "landrace", "heirloom", "organic", "sun-grown", "outdoor",
                "mountain-grown", "high-altitude"
            ],
        }
        
        for elem_type, patterns in detection_patterns.items():
            for pattern in patterns:
                if pattern in formulation_lower:
                    elements.append(KnowledgeElement(
                        element_type=elem_type,
                        description=f"Detected '{pattern}' in formulation",
                        source_community="pending_attribution",
                        documentation_level="detected",
                        contribution_weight=self._ELEMENT_TYPE_WEIGHTS[elem_type],
                    ))
                    break  # One element per type max
        
        # Check compounds for landrace/traditional strain origins
        if compounds:
            for compound in compounds:
                if any(x in compound.lower() for x in ["cbd", "thc", "cbg", "cbn"]):
                    elements.append(KnowledgeElement(
                        element_type=KnowledgeElementType.SYNERGY_KNOWLEDGE,
                        description=f"Cannabinoid {compound} - traditional knowledge of effects",
                        source_community="pending_attribution",
                        documentation_level="established",
                        contribution_weight=0.10,
                    ))
                    break
        
        return elements
    
    def _match_communities(
        self,
        elements: List[KnowledgeElement],
        claimed_communities: Optional[List[str]],
    ) -> Dict[str, List[KnowledgeElement]]:
        """
        Match knowledge elements to contributing communities.
        Trade Secret: Community matching algorithm.
        """
        community_elements: Dict[str, List[KnowledgeElement]] = {}
        
        # If communities are claimed, use those
        if claimed_communities:
            for comm_id in claimed_communities:
                if comm_id in self._KNOWN_COMMUNITIES:
                    community_elements[comm_id] = []
                    comm_info = self._KNOWN_COMMUNITIES[comm_id]
                    
                    # Assign elements matching community's primary knowledge
                    for elem in elements:
                        if elem.element_type in comm_info.get("primary_knowledge", []):
                            elem.source_community = comm_info["name"]
                            community_elements[comm_id].append(elem)
            return community_elements
        
        # Auto-match based on element types
        for elem in elements:
            for comm_id, comm_info in self._KNOWN_COMMUNITIES.items():
                if elem.element_type in comm_info.get("primary_knowledge", []):
                    if comm_id not in community_elements:
                        community_elements[comm_id] = []
                    elem.source_community = comm_info["name"]
                    community_elements[comm_id].append(elem)
                    break  # First match wins
        
        return community_elements
    
    def _calculate_compensation_suggestions(self, contribution: TKContribution) -> None:
        """
        Calculate suggested compensation based on TK contribution level.
        Trade Secret: Compensation tier calculation algorithm.
        """
        total_contribution = contribution.total_tk_contribution_percentage / 100
        
        # Find appropriate tier
        selected_tier = None
        for tier_name, tier_config in sorted(
            self._COMPENSATION_TIERS.items(),
            key=lambda x: x[1]["threshold"],
            reverse=True
        ):
            if total_contribution >= tier_config["threshold"]:
                selected_tier = tier_config
                break
        
        if selected_tier:
            contribution.suggested_compensation_percentage = selected_tier["percentage"]
            contribution.suggested_minimum_payment = selected_tier["minimum"]
        else:
            # Below minimal threshold - still suggest token compensation
            contribution.suggested_compensation_percentage = 1.0
            contribution.suggested_minimum_payment = 25.0
    
    def _check_nagoya_compliance(self, contribution: TKContribution) -> bool:
        """
        Check Nagoya Protocol compliance.
        Requires: PIC (Prior Informed Consent), MAT (Mutually Agreed Terms).
        """
        # Check each community has consent documentation
        for comm in contribution.communities:
            if comm.nagoya_pico_reference is None:
                return False
        return len(contribution.communities) > 0
    
    def _check_undrip_compliance(self, contribution: TKContribution) -> bool:
        """
        Check UN Declaration on Rights of Indigenous Peoples compliance.
        Articles 31 & 32 require FPIC (Free, Prior, Informed Consent).
        """
        for comm in contribution.communities:
            if not comm.undrip_acknowledgment:
                return False
        return len(contribution.communities) > 0
    
    def get_attribution_audit(self, contribution_id: str) -> Optional[TKContribution]:
        """Retrieve cached attribution for audit purposes."""
        return self._attribution_cache.get(contribution_id)
    
    def list_known_communities(self) -> List[Dict]:
        """List all known communities for UI selection."""
        return [
            {
                "id": comm_id,
                "name": info["name"],
                "region": info["region"],
                "country": info["country"],
                "primary_knowledge_types": [kt.value for kt in info["primary_knowledge"]],
            }
            for comm_id, info in self._KNOWN_COMMUNITIES.items()
        ]
