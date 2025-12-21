"""
NeuroBotanica Feature Flags & Ethical Pricing System

Enables/disables features per customer with ethical tiered pricing.
TK (Traditional Knowledge) features provide DISCOUNTS - customers who
opt OUT of TK attribution/compensation pay a PREMIUM surcharge.

This incentivizes ethical behavior by making it financially advantageous
to properly attribute and compensate traditional knowledge holders.

Trade Secret: Pricing tiers, ethical surcharge calculations, discount logic.
"""

from enum import Enum
from typing import Dict, List, Optional, Set
from pydantic import BaseModel, Field
from datetime import datetime
import uuid


class SubscriptionTier(str, Enum):
    """Base subscription tiers."""
    BASIC = "basic"
    STANDARD = "standard"
    PREMIUM = "premium"
    ENTERPRISE = "enterprise"


class FeatureCategory(str, Enum):
    """Feature categories for organization."""
    CORE = "core"                    # Always included
    ANALYSIS = "analysis"            # ChemPath, ToxPath, etc.
    REGULATORY = "regulatory"        # RegPath, ClinPath
    IP_PROTECTION = "ip_protection"  # PatentPath
    GENOMICS = "genomics"            # GenomePath, BioPath
    TK_OPTIONAL = "tk_optional"      # Traditional Knowledge (premium add-on)
    INFRASTRUCTURE = "infrastructure"  # MetaPath, TechPath


# =============================================================================
# Feature Definitions
# =============================================================================

class FeatureDefinition(BaseModel):
    """Definition of a feature module."""
    code: str
    name: str
    description: str
    category: FeatureCategory
    min_tier: SubscriptionTier
    is_addon: bool = False
    addon_price_monthly: float = 0.0
    requires: List[str] = Field(default_factory=list)  # Dependencies


# Core features (included in all tiers)
CORE_FEATURES = {
    "api_access": FeatureDefinition(
        code="api_access",
        name="API Access",
        description="Base API access with rate limiting",
        category=FeatureCategory.CORE,
        min_tier=SubscriptionTier.BASIC,
    ),
    "health_endpoints": FeatureDefinition(
        code="health_endpoints",
        name="Health & Status Endpoints",
        description="System health and status monitoring",
        category=FeatureCategory.CORE,
        min_tier=SubscriptionTier.BASIC,
    ),
}

# Analysis features
ANALYSIS_FEATURES = {
    "chempath_basic": FeatureDefinition(
        code="chempath_basic",
        name="ChemPath Basic",
        description="Basic molecular characterization (2D descriptors)",
        category=FeatureCategory.ANALYSIS,
        min_tier=SubscriptionTier.BASIC,
    ),
    "chempath_full": FeatureDefinition(
        code="chempath_full",
        name="ChemPath Full",
        description="Full characterization with 3D conformers and COA analysis",
        category=FeatureCategory.ANALYSIS,
        min_tier=SubscriptionTier.STANDARD,
        requires=["chempath_basic"],
    ),
    "toxpath_basic": FeatureDefinition(
        code="toxpath_basic",
        name="ToxPath Basic",
        description="Basic toxicity risk assessment",
        category=FeatureCategory.ANALYSIS,
        min_tier=SubscriptionTier.BASIC,
    ),
    "toxpath_full": FeatureDefinition(
        code="toxpath_full",
        name="ToxPath Full",
        description="Full toxicity assessment with testing plans and memos",
        category=FeatureCategory.ANALYSIS,
        min_tier=SubscriptionTier.STANDARD,
        requires=["toxpath_basic"],
    ),
    "terpene_analysis": FeatureDefinition(
        code="terpene_analysis",
        name="Terpene Analysis",
        description="Terpene profiling and effect predictions",
        category=FeatureCategory.ANALYSIS,
        min_tier=SubscriptionTier.STANDARD,
    ),
}

# Regulatory features
REGULATORY_FEATURES = {
    "regpath": FeatureDefinition(
        code="regpath",
        name="RegPath",
        description="FDA regulatory pathway guidance and memos",
        category=FeatureCategory.REGULATORY,
        min_tier=SubscriptionTier.STANDARD,
    ),
    "clinpath": FeatureDefinition(
        code="clinpath",
        name="ClinPath",
        description="Adaptive clinical trial design optimization",
        category=FeatureCategory.REGULATORY,
        min_tier=SubscriptionTier.PREMIUM,
        requires=["regpath"],
    ),
    "fda_docs": FeatureDefinition(
        code="fda_docs",
        name="FDA Document Generation",
        description="Generate FDA submission document templates",
        category=FeatureCategory.REGULATORY,
        min_tier=SubscriptionTier.PREMIUM,
        requires=["regpath"],
    ),
}

# IP Protection features
IP_FEATURES = {
    "patentpath_basic": FeatureDefinition(
        code="patentpath_basic",
        name="PatentPath Basic",
        description="Prior art search and basic novelty assessment",
        category=FeatureCategory.IP_PROTECTION,
        min_tier=SubscriptionTier.STANDARD,
    ),
    "patentpath_full": FeatureDefinition(
        code="patentpath_full",
        name="PatentPath Full",
        description="Full IP suite: FTO, claims, cost estimation",
        category=FeatureCategory.IP_PROTECTION,
        min_tier=SubscriptionTier.PREMIUM,
        requires=["patentpath_basic"],
    ),
}

# Genomics features
GENOMICS_FEATURES = {
    "genomepath": FeatureDefinition(
        code="genomepath",
        name="GenomePath",
        description="Genomic target prediction and pathway analysis",
        category=FeatureCategory.GENOMICS,
        min_tier=SubscriptionTier.PREMIUM,
    ),
    "biopath": FeatureDefinition(
        code="biopath",
        name="BioPath",
        description="Bias-corrected efficacy validation",
        category=FeatureCategory.GENOMICS,
        min_tier=SubscriptionTier.PREMIUM,
        requires=["genomepath"],
    ),
}

# Traditional Knowledge features (ETHICAL DISCOUNT - opt-out pays MORE)
TK_FEATURES = {
    "tk_attribution": FeatureDefinition(
        code="tk_attribution",
        name="TK Attribution",
        description="Traditional knowledge source attribution in patents - INCLUDED (opt-out adds surcharge)",
        category=FeatureCategory.TK_OPTIONAL,
        min_tier=SubscriptionTier.STANDARD,
        is_addon=False,  # Included by default
        addon_price_monthly=0.0,  # Free when enabled
    ),
    "tk_sourcing": FeatureDefinition(
        code="tk_sourcing",
        name="TK Sourcing",
        description="Access traditional knowledge databases with proper attribution - INCLUDED",
        category=FeatureCategory.TK_OPTIONAL,
        min_tier=SubscriptionTier.STANDARD,
        is_addon=False,
        addon_price_monthly=0.0,
    ),
    "equipath": FeatureDefinition(
        code="equipath",
        name="EquiPath",
        description="Zero-knowledge proof compensation for TK holders - INCLUDED",
        category=FeatureCategory.TK_OPTIONAL,
        min_tier=SubscriptionTier.PREMIUM,
        is_addon=False,
        addon_price_monthly=0.0,
        requires=["tk_attribution"],
    ),
    "federated_tk": FeatureDefinition(
        code="federated_tk",
        name="Federated TK",
        description="Cross-community TK data federation via TechPath - INCLUDED",
        category=FeatureCategory.TK_OPTIONAL,
        min_tier=SubscriptionTier.PREMIUM,
        is_addon=False,
        addon_price_monthly=0.0,
        requires=["equipath"],
    ),
    "tk_full_suite": FeatureDefinition(
        code="tk_full_suite",
        name="TK Full Suite",
        description="All TK features - DEFAULT ETHICAL PACKAGE",
        category=FeatureCategory.TK_OPTIONAL,
        min_tier=SubscriptionTier.STANDARD,
        is_addon=False,
        addon_price_monthly=0.0,
    ),
}

# =============================================================================
# Ethical Surcharges (for opting OUT of TK features)
# =============================================================================

class TKOptOutSurcharge(BaseModel):
    """Surcharge applied when customer opts out of TK features."""
    code: str
    name: str
    description: str
    surcharge_percentage: float  # Percentage added to base price
    surcharge_flat: float  # Flat monthly fee added


TK_OPT_OUT_SURCHARGES = {
    "no_tk_attribution": TKOptOutSurcharge(
        code="no_tk_attribution",
        name="No TK Attribution Surcharge",
        description="Surcharge for not attributing traditional knowledge sources",
        surcharge_percentage=15.0,  # +15% of base price
        surcharge_flat=150.0,
    ),
    "no_tk_compensation": TKOptOutSurcharge(
        code="no_tk_compensation",
        name="No TK Compensation Surcharge",
        description="Surcharge for not compensating traditional knowledge holders",
        surcharge_percentage=25.0,  # +25% of base price
        surcharge_flat=250.0,
    ),
    "no_tk_full": TKOptOutSurcharge(
        code="no_tk_full",
        name="Full TK Opt-Out Surcharge",
        description="Maximum surcharge for opting out of all TK features",
        surcharge_percentage=40.0,  # +40% of base price
        surcharge_flat=500.0,
    ),
}

# Infrastructure features
INFRASTRUCTURE_FEATURES = {
    "metapath": FeatureDefinition(
        code="metapath",
        name="MetaPath",
        description="Cross-pathway orchestration for complex analyses",
        category=FeatureCategory.INFRASTRUCTURE,
        min_tier=SubscriptionTier.PREMIUM,
    ),
    "techpath": FeatureDefinition(
        code="techpath",
        name="TechPath",
        description="Federated AI infrastructure",
        category=FeatureCategory.INFRASTRUCTURE,
        min_tier=SubscriptionTier.ENTERPRISE,
    ),
}

# Combine all features
ALL_FEATURES: Dict[str, FeatureDefinition] = {
    **CORE_FEATURES,
    **ANALYSIS_FEATURES,
    **REGULATORY_FEATURES,
    **IP_FEATURES,
    **GENOMICS_FEATURES,
    **TK_FEATURES,
    **INFRASTRUCTURE_FEATURES,
}


# =============================================================================
# Tier Definitions
# =============================================================================

class TierDefinition(BaseModel):
    """Definition of a subscription tier."""
    tier: SubscriptionTier
    name: str
    base_price_monthly: float
    included_features: List[str]
    rate_limit_per_minute: int
    rate_limit_per_day: int
    support_level: str
    tk_addons_available: bool


TIER_DEFINITIONS: Dict[SubscriptionTier, TierDefinition] = {
    SubscriptionTier.BASIC: TierDefinition(
        tier=SubscriptionTier.BASIC,
        name="Basic",
        base_price_monthly=99.0,
        included_features=[
            "api_access",
            "health_endpoints",
            "chempath_basic",
            "toxpath_basic",
        ],
        rate_limit_per_minute=30,
        rate_limit_per_day=1000,
        support_level="Community",
        tk_addons_available=False,
    ),
    SubscriptionTier.STANDARD: TierDefinition(
        tier=SubscriptionTier.STANDARD,
        name="Standard",
        base_price_monthly=299.0,
        included_features=[
            "api_access",
            "health_endpoints",
            "chempath_basic",
            "chempath_full",
            "toxpath_basic",
            "toxpath_full",
            "terpene_analysis",
            "regpath",
            "patentpath_basic",
        ],
        rate_limit_per_minute=60,
        rate_limit_per_day=5000,
        support_level="Email (48h)",
        tk_addons_available=True,
    ),
    SubscriptionTier.PREMIUM: TierDefinition(
        tier=SubscriptionTier.PREMIUM,
        name="Premium",
        base_price_monthly=799.0,
        included_features=[
            "api_access",
            "health_endpoints",
            "chempath_basic",
            "chempath_full",
            "toxpath_basic",
            "toxpath_full",
            "terpene_analysis",
            "regpath",
            "clinpath",
            "fda_docs",
            "patentpath_basic",
            "patentpath_full",
            "genomepath",
            "biopath",
            "metapath",
        ],
        rate_limit_per_minute=120,
        rate_limit_per_day=20000,
        support_level="Priority (24h)",
        tk_addons_available=True,
    ),
    SubscriptionTier.ENTERPRISE: TierDefinition(
        tier=SubscriptionTier.ENTERPRISE,
        name="Enterprise",
        base_price_monthly=0.0,  # Custom pricing
        included_features=list(ALL_FEATURES.keys()),  # All features
        rate_limit_per_minute=1000,
        rate_limit_per_day=100000,
        support_level="Dedicated Account Manager",
        tk_addons_available=True,
    ),
}


# =============================================================================
# Customer Subscription Model
# =============================================================================

class CustomerSubscription(BaseModel):
    """Customer subscription with tier and TK opt-out tracking."""
    subscription_id: str = Field(default_factory=lambda: str(uuid.uuid4()))
    customer_id: str
    tier: SubscriptionTier
    active_addons: List[str] = Field(default_factory=list)
    custom_features: List[str] = Field(default_factory=list)  # Enterprise only
    
    # TK Opt-Out Tracking (opt-out = pay more)
    tk_opted_out: bool = False  # If True, customer pays surcharge
    tk_opt_out_level: Optional[str] = None  # "no_tk_attribution", "no_tk_compensation", "no_tk_full"
    tk_opt_out_reason: Optional[str] = None  # Required justification
    tk_opt_out_acknowledged_at: Optional[datetime] = None
    
    created_at: datetime = Field(default_factory=datetime.utcnow)
    expires_at: Optional[datetime] = None
    
    def get_enabled_features(self) -> Set[str]:
        """Get all features enabled for this subscription."""
        tier_def = TIER_DEFINITIONS[self.tier]
        features = set(tier_def.included_features)
        features.update(self.active_addons)
        features.update(self.custom_features)
        
        # TK features included by default unless opted out
        if not self.tk_opted_out:
            features.update(TK_FEATURES.keys())
        
        return features
    
    def has_feature(self, feature_code: str) -> bool:
        """Check if subscription has access to a feature."""
        return feature_code in self.get_enabled_features()
    
    def has_tk_features(self) -> bool:
        """Check if TK features are enabled (default: yes)."""
        return not self.tk_opted_out
    
    def calculate_monthly_price(self) -> float:
        """Calculate total monthly price including ethical surcharges."""
        tier_def = TIER_DEFINITIONS[self.tier]
        base = tier_def.base_price_monthly
        
        # Add any custom add-ons
        addon_total = 0.0
        for addon_code in self.active_addons:
            if addon_code in ALL_FEATURES:
                addon_total += ALL_FEATURES[addon_code].addon_price_monthly
        
        subtotal = base + addon_total
        
        # Apply TK opt-out surcharge if applicable
        if self.tk_opted_out and self.tk_opt_out_level:
            surcharge = TK_OPT_OUT_SURCHARGES.get(self.tk_opt_out_level)
            if surcharge:
                percentage_surcharge = subtotal * (surcharge.surcharge_percentage / 100)
                subtotal += percentage_surcharge + surcharge.surcharge_flat
        
        return subtotal
    
    def get_tk_discount_amount(self) -> float:
        """
        Calculate potential savings from enabling TK features.
        
        If opted OUT: Returns how much they could save by re-enabling TK
        If opted IN: Returns 0 (they're already saving)
        """
        if not self.tk_opted_out or not self.tk_opt_out_level:
            return 0.0
        
        # Calculate the surcharge they're currently paying
        tier_def = TIER_DEFINITIONS[self.tier]
        base = tier_def.base_price_monthly
        
        surcharge = TK_OPT_OUT_SURCHARGES.get(self.tk_opt_out_level)
        if surcharge:
            return base * (surcharge.surcharge_percentage / 100) + surcharge.surcharge_flat
        
        return 0.0


# =============================================================================
# Feature Access Manager
# =============================================================================

class FeatureAccessManager:
    """
    Manages feature access for customers with ethical pricing.
    
    TK features are INCLUDED by default. Customers who opt-out
    pay a surcharge, incentivizing ethical behavior.
    
    Trade Secret: Access logic, surcharge calculations, opt-out workflows.
    """
    
    def __init__(self):
        self._subscriptions: Dict[str, CustomerSubscription] = {}
    
    def register_subscription(self, subscription: CustomerSubscription) -> None:
        """Register a customer subscription."""
        self._subscriptions[subscription.customer_id] = subscription
    
    def get_subscription(self, customer_id: str) -> Optional[CustomerSubscription]:
        """Get customer subscription."""
        return self._subscriptions.get(customer_id)
    
    def check_access(self, customer_id: str, feature_code: str) -> Dict:
        """
        Check if customer has access to a feature.
        
        Returns:
            {
                "allowed": bool,
                "reason": str,
                "upgrade_options": List[str] if not allowed
            }
        """
        subscription = self._subscriptions.get(customer_id)
        
        if not subscription:
            return {
                "allowed": False,
                "reason": "No active subscription",
                "upgrade_options": ["basic", "standard", "premium"],
            }
        
        if feature_code not in ALL_FEATURES:
            return {
                "allowed": False,
                "reason": f"Unknown feature: {feature_code}",
                "upgrade_options": [],
            }
        
        feature = ALL_FEATURES[feature_code]
        
        # Check if feature is enabled
        if subscription.has_feature(feature_code):
            # Check dependencies
            missing_deps = []
            for dep in feature.requires:
                if not subscription.has_feature(dep):
                    missing_deps.append(dep)
            
            if missing_deps:
                return {
                    "allowed": False,
                    "reason": f"Missing required features: {missing_deps}",
                    "upgrade_options": missing_deps,
                }
            
            return {
                "allowed": True,
                "reason": "Feature enabled",
                "upgrade_options": [],
            }
        
        # Feature not enabled - check if it's a TK feature they opted out of
        if feature_code in TK_FEATURES and subscription.tk_opted_out:
            return {
                "allowed": False,
                "reason": f"TK feature '{feature.name}' disabled due to opt-out. Re-enable TK features to access.",
                "upgrade_options": ["enable_tk_features"],
            }
        
        # Tier-locked feature
        upgrade_options = []
        for tier in SubscriptionTier:
            tier_def = TIER_DEFINITIONS[tier]
            if feature_code in tier_def.included_features:
                if list(SubscriptionTier).index(tier) > list(SubscriptionTier).index(subscription.tier):
                    upgrade_options.append(f"upgrade_tier:{tier.value}")
                break
        
        return {
            "allowed": False,
            "reason": f"Feature '{feature.name}' not included in subscription",
            "upgrade_options": upgrade_options,
        }
    
    def request_tk_opt_out(
        self, 
        customer_id: str, 
        opt_out_level: str,
        reason: str
    ) -> Dict:
        """
        Process request to opt-out of TK features (triggers surcharge).
        
        Requires explicit acknowledgment and reason.
        """
        subscription = self._subscriptions.get(customer_id)
        if not subscription:
            return {"success": False, "error": "No subscription found"}
        
        if opt_out_level not in TK_OPT_OUT_SURCHARGES:
            return {
                "success": False,
                "error": f"Invalid opt-out level. Options: {list(TK_OPT_OUT_SURCHARGES.keys())}",
            }
        
        if not reason or len(reason) < 20:
            return {
                "success": False,
                "error": "A detailed reason (minimum 20 characters) is required for TK opt-out",
            }
        
        surcharge = TK_OPT_OUT_SURCHARGES[opt_out_level]
        tier_def = TIER_DEFINITIONS[subscription.tier]
        base_price = tier_def.base_price_monthly
        
        # Calculate new price with surcharge
        percentage_add = base_price * (surcharge.surcharge_percentage / 100)
        new_price = base_price + percentage_add + surcharge.surcharge_flat
        price_increase = new_price - base_price
        
        return {
            "success": True,
            "warning": "TK opt-out will increase your monthly cost",
            "opt_out_level": opt_out_level,
            "surcharge_name": surcharge.name,
            "surcharge_description": surcharge.description,
            "current_price": base_price,
            "new_price": new_price,
            "price_increase": price_increase,
            "percentage_increase": surcharge.surcharge_percentage,
            "flat_surcharge": surcharge.surcharge_flat,
            "confirmation_required": True,
            "confirmation_message": (
                f"By opting out of TK features, you acknowledge that:\n"
                f"1. Your monthly cost will increase by ${price_increase:.2f} ({surcharge.surcharge_percentage}% + ${surcharge.surcharge_flat})\n"
                f"2. You will not have access to traditional knowledge attribution features\n"
                f"3. Your patents/products will not include TK compensation mechanisms\n"
                f"4. This may affect your ethical standing and regulatory compliance in some jurisdictions\n\n"
                f"Reason provided: {reason}"
            ),
        }
    
    def confirm_tk_opt_out(
        self,
        customer_id: str,
        opt_out_level: str,
        reason: str,
        confirmed: bool
    ) -> Dict:
        """Confirm and apply TK opt-out with surcharge."""
        if not confirmed:
            return {"success": False, "error": "Opt-out not confirmed"}
        
        subscription = self._subscriptions.get(customer_id)
        if not subscription:
            return {"success": False, "error": "No subscription found"}
        
        subscription.tk_opted_out = True
        subscription.tk_opt_out_level = opt_out_level
        subscription.tk_opt_out_reason = reason
        subscription.tk_opt_out_acknowledged_at = datetime.utcnow()
        
        return {
            "success": True,
            "message": "TK opt-out confirmed. Surcharge applied.",
            "new_monthly_price": subscription.calculate_monthly_price(),
            "tk_features_disabled": list(TK_FEATURES.keys()),
        }
    
    def enable_tk_features(self, customer_id: str) -> Dict:
        """Re-enable TK features (removes surcharge)."""
        subscription = self._subscriptions.get(customer_id)
        if not subscription:
            return {"success": False, "error": "No subscription found"}
        
        if not subscription.tk_opted_out:
            return {"success": False, "error": "TK features already enabled"}
        
        old_price = subscription.calculate_monthly_price()
        
        subscription.tk_opted_out = False
        subscription.tk_opt_out_level = None
        subscription.tk_opt_out_reason = None
        subscription.tk_opt_out_acknowledged_at = None
        
        new_price = subscription.calculate_monthly_price()
        
        return {
            "success": True,
            "message": "TK features re-enabled! Surcharge removed.",
            "old_price": old_price,
            "new_price": new_price,
            "savings": old_price - new_price,
            "tk_features_enabled": list(TK_FEATURES.keys()),
        }
    
    def get_pricing_summary(self, customer_id: str) -> Dict:
        """Get pricing summary for customer including ethical pricing details."""
        subscription = self._subscriptions.get(customer_id)
        if not subscription:
            return {"error": "No subscription found"}
        
        tier_def = TIER_DEFINITIONS[subscription.tier]
        
        addon_details = []
        for addon_code in subscription.active_addons:
            if addon_code in ALL_FEATURES:
                feature = ALL_FEATURES[addon_code]
                addon_details.append({
                    "code": addon_code,
                    "name": feature.name,
                    "price": feature.addon_price_monthly,
                })
        
        # Calculate surcharge details if applicable
        surcharge_details = None
        if subscription.tk_opted_out and subscription.tk_opt_out_level:
            surcharge = TK_OPT_OUT_SURCHARGES.get(subscription.tk_opt_out_level)
            if surcharge:
                surcharge_details = {
                    "level": subscription.tk_opt_out_level,
                    "name": surcharge.name,
                    "percentage": surcharge.surcharge_percentage,
                    "flat_fee": surcharge.surcharge_flat,
                    "reason": subscription.tk_opt_out_reason,
                    "opted_out_at": subscription.tk_opt_out_acknowledged_at.isoformat() if subscription.tk_opt_out_acknowledged_at else None,
                }
        
        return {
            "tier": subscription.tier.value,
            "tier_name": tier_def.name,
            "base_price": tier_def.base_price_monthly,
            "active_addons": addon_details,
            "addon_total": sum(a["price"] for a in addon_details),
            "tk_features_enabled": not subscription.tk_opted_out,
            "tk_opt_out_surcharge": surcharge_details,
            "total_monthly": subscription.calculate_monthly_price(),
            "ethical_discount_if_tk_enabled": subscription.get_tk_discount_amount() if subscription.tk_opted_out else 0,
        }


# =============================================================================
# Singleton Instance
# =============================================================================

_feature_manager: Optional[FeatureAccessManager] = None


def get_feature_manager() -> FeatureAccessManager:
    """Get singleton feature manager instance."""
    global _feature_manager
    if _feature_manager is None:
        _feature_manager = FeatureAccessManager()
    return _feature_manager
