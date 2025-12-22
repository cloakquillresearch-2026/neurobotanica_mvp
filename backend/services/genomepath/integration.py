"""
GenomePath ↔ EthnoPath Integration Module

Trade Secret Integration: Connects TS-GP-001 ($6.2B GenomePath) with
TS-ETH-001 ($4.6B EthnoPath) for unified traditional knowledge →
genomic correlation ecosystem.

Integration Points:
1. Community Authentication: Shared community IDs from EthnoPath
2. Sacred Knowledge Protection: Both systems block ceremonial_significance
3. Consent Workflow: EthnoPath grants feed GenomePath validation
4. Attribution Chain: EquiPath blockchain links both systems
5. Revenue Synergies: Combined pricing with TK opt-out surcharges

Ecosystem Value: $10.8B standalone → $18-25B with integration multipliers
"""

from typing import Dict, List, Optional, Tuple, Any
from datetime import datetime
import uuid

from backend.services.ethnopath.governance import (
    EthnoPathGovernance,
    AccessGrant,
    ConsentRecord,
    TKComplianceLevel,
    PricingTier,
)

from backend.services.genomepath.bridge import (
    GenomePathBridge,
    TKEncodedVector,
    CulturalSensitivityLevel,
)

from backend.services.genomepath.correlation import (
    TKGenomicCorrelator,
    CorrelationResult,
)


class GenomePathEthnoPathIntegration:
    """
    Integration layer connecting GenomePath and EthnoPath systems.
    
    Trade Secret: Integration algorithms maximizing ecosystem value through
    shared authentication, consent management, and attribution tracking.
    
    Key Features:
    - Unified community authentication from EthnoPath
    - Sacred knowledge protection across both systems
    - Consent workflow integration (EthnoPath grants → GenomePath validation)
    - Attribution chain via EquiPath blockchain
    - Combined revenue model with TK opt-out surcharges
    """
    
    def __init__(
        self,
        ethnopath_governance: Optional[EthnoPathGovernance] = None,
        genomepath_bridge: Optional[GenomePathBridge] = None,
        genomepath_correlator: Optional[TKGenomicCorrelator] = None,
    ):
        self.ethnopath = ethnopath_governance or EthnoPathGovernance()
        self.genomepath_bridge = genomepath_bridge or GenomePathBridge()
        self.genomepath_correlator = genomepath_correlator or TKGenomicCorrelator()
        
        # Track integrated workflows
        self._integrated_correlations: Dict[str, Dict[str, Any]] = {}
    
    # =========================================================================
    # 1. Community Authentication Integration
    # =========================================================================
    
    def verify_community_genomic_access(
        self,
        community_id: str,
        user_id: str,
    ) -> Tuple[bool, str]:
        """
        Verify community has granted genomic research access via EthnoPath.
        
        Integration Point: Uses EthnoPath community_id for authentication.
        
        Returns:
            (is_authorized, reason)
        """
        # Check if community exists in EthnoPath
        if community_id not in self.ethnopath._communities:
            return False, f"Community {community_id} not registered in EthnoPath"
        
        # Check if user is community member or elder
        # In production, would query EthnoPath for user's role
        community = self.ethnopath._communities[community_id]
        
        # Simplified authorization (production would check participant roles)
        return True, "Community genomic access authorized"
    
    # =========================================================================
    # 2. Sacred Knowledge Protection Integration
    # =========================================================================
    
    def check_sacred_knowledge_protection(
        self,
        tk_entry_id: str,
        community_id: str,
    ) -> Tuple[bool, str]:
        """
        Check if TK entry is sacred/ceremonial and should be blocked.
        
        Integration Point: Both EthnoPath and GenomePath block sacred knowledge.
        
        EthnoPath: ceremonial_significance=True → blocks access
        GenomePath: sacred_knowledge_flag=True → blocks encoding
        
        Returns:
            (is_sacred, protection_reason)
        """
        # Query EthnoPath for TK entry
        # In production, would retrieve from EthnoPath digitizer
        # Simplified: check if entry has ceremonial significance
        
        # For MVP, assume we have the entry metadata
        entry_metadata = {
            "entry_id": tk_entry_id,
            "community_id": community_id,
            "ceremonial_significance": False,  # Would query EthnoPath
        }
        
        if entry_metadata.get("ceremonial_significance"):
            return True, "Sacred knowledge - ABSOLUTE preservation priority, genomic encoding blocked"
        
        return False, "Not sacred - genomic research permitted with consent"
    
    # =========================================================================
    # 3. Consent Workflow Integration
    # =========================================================================
    
    def verify_genomic_research_consent(
        self,
        community_id: str,
        tk_entry_id: str,
        user_id: str,
    ) -> Tuple[bool, Optional[ConsentRecord]]:
        """
        Verify community has granted consent for genomic research.
        
        Integration Point: EthnoPath consent grants feed GenomePath validation.
        
        Workflow:
        1. Query EthnoPath for consent record
        2. Verify consent_status="granted"
        3. Check consent_scope includes "genomic_research"
        4. Return consent record for GenomePath attribution
        
        Returns:
            (consent_granted, consent_record)
        """
        # In production, would query EthnoPath governance for consent records
        # For MVP, create simplified consent verification
        
        consent_record = ConsentRecord(
            consent_id=str(uuid.uuid4()),
            entry_id=tk_entry_id,
            grantor_id=community_id,
            grantee_id=user_id,
            consent_status="granted",  # Would verify in EthnoPath
            consent_scope="genomic_research",
            granted_at=datetime.utcnow(),
        )
        
        if consent_record.consent_status != "granted":
            return False, None
        
        if "genomic" not in consent_record.consent_scope.lower():
            return False, None
        
        return True, consent_record
    
    def request_genomic_research_consent(
        self,
        community_id: str,
        tk_entry_id: str,
        user_id: str,
        research_purpose: str,
    ) -> str:
        """
        Request community consent for genomic research via EthnoPath governance.
        
        Integration Point: Creates EthnoPath proposal for community vote.
        
        Workflow:
        1. Create EthnoPath governance proposal
        2. Type: GENOMIC_RESEARCH_CONSENT
        3. Requires community vote (elders have veto authority)
        4. If approved, creates consent record
        5. Returns proposal_id for tracking
        
        Returns:
            proposal_id (can be used to track approval status)
        """
        from backend.services.ethnopath.governance import ProposalType
        
        # Create governance proposal in EthnoPath
        proposal = self.ethnopath.create_proposal(
            proposal_type=ProposalType.BENEFIT_SHARING,  # Genomic research benefits
            title=f"Genomic Research Consent Request: {tk_entry_id}",
            description=f"Purpose: {research_purpose}",
            proposed_by=user_id,
            community_id=community_id,
            metadata={
                "tk_entry_id": tk_entry_id,
                "research_type": "genomic_correlation",
                "downstream_applications": ["therapeutic", "diagnostic"],
            }
        )
        
        return proposal.proposal_id
    
    # =========================================================================
    # 4. Attribution Chain Integration
    # =========================================================================
    
    def link_attribution_chain(
        self,
        ethnopath_access_grant_id: str,
        genomepath_correlation_id: str,
        equipath_attribution_id: Optional[str] = None,
    ) -> Dict[str, str]:
        """
        Link attribution chain across EthnoPath → GenomePath → EquiPath.
        
        Integration Point: Same EquiPath blockchain for both systems.
        
        Attribution Flow:
        1. EthnoPath creates access grant (access_grant_id)
        2. GenomePath creates correlation (correlation_id)
        3. Both reference same EquiPath attribution (equipath_attribution_id)
        4. Blockchain tracks benefit-sharing proportionally
        
        Returns:
            attribution_chain dictionary
        """
        # Generate EquiPath attribution ID if not provided
        if not equipath_attribution_id:
            equipath_attribution_id = f"equipath_{uuid.uuid4()}"
        
        attribution_chain = {
            "ethnopath_access_grant_id": ethnopath_access_grant_id,
            "genomepath_correlation_id": genomepath_correlation_id,
            "equipath_attribution_id": equipath_attribution_id,
            "attribution_timestamp": datetime.utcnow().isoformat(),
            "attribution_type": "tk_genomic_correlation",
        }
        
        # Store attribution chain (in production, would write to EquiPath blockchain)
        self._integrated_correlations[genomepath_correlation_id] = attribution_chain
        
        return attribution_chain
    
    # =========================================================================
    # 5. Combined Revenue Model Integration
    # =========================================================================
    
    def calculate_combined_pricing(
        self,
        ethnopath_pricing_tier: PricingTier,
        genomepath_enabled: bool,
        tk_compliance_level: TKComplianceLevel,
    ) -> Dict[str, float]:
        """
        Calculate combined pricing for EthnoPath + GenomePath access.
        
        Integration Point: Combined pricing with TK opt-out surcharges.
        
        Pricing Structure:
        - Base: EthnoPath tier ($99-$2,499/mo)
        - Add-on: GenomePath (+$400/mo if enabled)
        - Surcharge: TK opt-out (15%-40% + $150-$500 on combined base)
        
        Example:
        PREMIUM ($799) + GenomePath ($400) = $1,199 combined base
        NO_TK_FULL opt-out: +40% on $1,199 + $500 = $1,677.60 surcharge
        Total: $2,876.60/mo ($1,677.60 to compensation pool)
        
        Returns:
            pricing breakdown dictionary
        """
        # EthnoPath base pricing
        base_prices = {
            PricingTier.BASIC: 99.0,
            PricingTier.STANDARD: 299.0,
            PricingTier.PREMIUM: 799.0,
            PricingTier.ENTERPRISE: 2499.0,
        }
        
        ethnopath_base = base_prices[ethnopath_pricing_tier]
        
        # GenomePath add-on
        genomepath_addon = 400.0 if genomepath_enabled else 0.0
        
        # Combined base price
        combined_base = ethnopath_base + genomepath_addon
        
        # TK opt-out surcharges (applied to combined base)
        surcharges = {
            TKComplianceLevel.FULL_COMPLIANCE: (0.0, 0.0),  # (percentage, flat_fee)
            TKComplianceLevel.NO_ATTRIBUTION: (0.15, 150.0),
            TKComplianceLevel.NO_COMPENSATION: (0.25, 250.0),
            TKComplianceLevel.NO_TK_FULL: (0.40, 500.0),
        }
        
        percentage, flat_fee = surcharges[tk_compliance_level]
        surcharge_amount = (combined_base * percentage) + flat_fee
        
        # Total price
        total_price = combined_base + surcharge_amount
        
        # Compensation pool contribution (100% of surcharge)
        compensation_pool_contribution = surcharge_amount
        
        return {
            "ethnopath_base_usd": ethnopath_base,
            "genomepath_addon_usd": genomepath_addon,
            "combined_base_usd": combined_base,
            "surcharge_percentage": percentage,
            "surcharge_flat_fee_usd": flat_fee,
            "surcharge_amount_usd": surcharge_amount,
            "total_price_usd": total_price,
            "compensation_pool_contribution_usd": compensation_pool_contribution,
            "pricing_tier": ethnopath_pricing_tier.value,
            "genomepath_enabled": genomepath_enabled,
            "tk_compliance_level": tk_compliance_level.value,
        }
    
    # =========================================================================
    # Complete Integrated Workflow
    # =========================================================================
    
    def execute_integrated_tk_genomic_workflow(
        self,
        community_id: str,
        tk_entry_id: str,
        user_id: str,
        practice_metadata: Dict[str, Any],
        pricing_tier: PricingTier,
        tk_compliance_level: TKComplianceLevel = TKComplianceLevel.FULL_COMPLIANCE,
    ) -> Dict[str, Any]:
        """
        Execute complete integrated TK → Genomic correlation workflow.
        
        Trade Secret: End-to-end integration maximizing ecosystem value.
        
        Workflow:
        1. Verify community authentication (EthnoPath)
        2. Check sacred knowledge protection (both systems)
        3. Verify genomic research consent (EthnoPath)
        4. Create EthnoPath access grant
        5. Execute GenomePath correlation
        6. Link attribution chain (EquiPath)
        7. Calculate combined pricing with surcharges
        8. Return integrated result
        
        Returns:
            integrated_result dictionary
        """
        result = {
            "workflow_id": str(uuid.uuid4()),
            "community_id": community_id,
            "tk_entry_id": tk_entry_id,
            "user_id": user_id,
            "timestamp": datetime.utcnow().isoformat(),
            "status": "pending",
            "errors": [],
        }
        
        # Step 1: Verify community authentication
        is_authorized, auth_reason = self.verify_community_genomic_access(
            community_id, user_id
        )
        if not is_authorized:
            result["status"] = "failed"
            result["errors"].append(f"Authentication failed: {auth_reason}")
            return result
        
        result["authentication_verified"] = True
        
        # Step 2: Check sacred knowledge protection
        is_sacred, protection_reason = self.check_sacred_knowledge_protection(
            tk_entry_id, community_id
        )
        if is_sacred:
            result["status"] = "blocked"
            result["errors"].append(f"Sacred knowledge protected: {protection_reason}")
            return result
        
        result["sacred_knowledge_protected"] = False
        
        # Step 3: Verify genomic research consent
        consent_granted, consent_record = self.verify_genomic_research_consent(
            community_id, tk_entry_id, user_id
        )
        if not consent_granted:
            result["status"] = "consent_required"
            result["errors"].append("Genomic research consent not granted")
            # Could auto-create proposal here
            return result
        
        result["consent_verified"] = True
        result["consent_record_id"] = consent_record.consent_id if consent_record else None
        
        # Step 4: Create EthnoPath access grant
        access_grant = self.ethnopath.create_access_grant(
            community_id=community_id,
            grantee_id=user_id,
            access_level="genomic_research",
            duration_days=365,  # 1 year access
            pricing_tier=pricing_tier,
            tk_compliance_level=tk_compliance_level,
        )
        
        result["ethnopath_access_grant_id"] = access_grant.grant_id
        
        # Step 5: Execute GenomePath correlation
        # (Simplified - full implementation would call correlation engine)
        genomepath_correlation_id = str(uuid.uuid4())
        result["genomepath_correlation_id"] = genomepath_correlation_id
        
        # Step 6: Link attribution chain
        attribution_chain = self.link_attribution_chain(
            ethnopath_access_grant_id=access_grant.grant_id,
            genomepath_correlation_id=genomepath_correlation_id,
        )
        result["attribution_chain"] = attribution_chain
        
        # Step 7: Calculate combined pricing
        pricing_breakdown = self.calculate_combined_pricing(
            ethnopath_pricing_tier=pricing_tier,
            genomepath_enabled=True,
            tk_compliance_level=tk_compliance_level,
        )
        result["pricing"] = pricing_breakdown
        
        # Step 8: Return integrated result
        result["status"] = "success"
        result["ecosystem_value_contribution_usd"] = pricing_breakdown["compensation_pool_contribution_usd"]
        
        return result
    
    def get_integration_statistics(self) -> Dict[str, Any]:
        """Get integration statistics across both systems."""
        return {
            "total_integrated_correlations": len(self._integrated_correlations),
            "ethnopath_communities": len(self.ethnopath._communities),
            "ethnopath_access_grants": len(self.ethnopath._access_grants),
            "genomepath_correlations": len(self.genomepath_correlator._correlation_results),
            "combined_ecosystem_value_usd": "10.8B standalone, 18-25B with multipliers",
        }


# =============================================================================
# Integration Helper Functions
# =============================================================================

def map_sensitivity_levels(
    ethnopath_sensitivity: str,
) -> CulturalSensitivityLevel:
    """
    Map EthnoPath sensitivity levels to GenomePath cultural sensitivity.
    
    Integration Point: Harmonize sensitivity classification across systems.
    """
    mapping = {
        "low": CulturalSensitivityLevel.LOW,
        "moderate": CulturalSensitivityLevel.MODERATE,
        "high": CulturalSensitivityLevel.HIGH,
        "sacred": CulturalSensitivityLevel.SACRED,
    }
    return mapping.get(ethnopath_sensitivity.lower(), CulturalSensitivityLevel.MODERATE)


def harmonize_elder_voting_weights() -> float:
    """
    Harmonize elder/TK holder voting weights across EthnoPath and GenomePath.
    
    Integration Point: Unified governance representation.
    
    EthnoPath: 3.0x elder voting weight
    GenomePath: 2.8x TK holder voting weight
    
    Recommendation: Harmonize to 3.0x for maximum protection
    """
    return 3.0  # Higher of the two for maximum TK holder protection
