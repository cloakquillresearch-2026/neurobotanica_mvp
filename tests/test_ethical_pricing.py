"""
Test Suite for NeuroBotanica Ethical Pricing System

Tests the ethical pricing model where TK (Traditional Knowledge) features
are INCLUDED by default and customers who opt-out pay a PREMIUM surcharge.

This incentivizes ethical behavior toward indigenous knowledge holders.
"""

import pytest
from datetime import datetime
from unittest.mock import patch

import sys
sys.path.insert(0, str(__file__).replace('/tests/test_ethical_pricing.py', '').replace('\\tests\\test_ethical_pricing.py', ''))

from backend.services.feature_flags import (
    SubscriptionTier,
    FeatureCategory,
    FeatureAccessManager,
    CustomerSubscription,
    TierDefinition,
    TIER_DEFINITIONS,
    TK_FEATURES,
    TK_OPT_OUT_SURCHARGES,
    ALL_FEATURES,
)


class TestEthicalPricingModel:
    """Tests for the ethical pricing model."""
    
    def test_tk_features_are_free_by_default(self):
        """TK features should be free (addon_price_monthly=0)."""
        for code, feature in TK_FEATURES.items():
            assert feature.addon_price_monthly == 0.0, f"{code} should be free"
            assert feature.is_addon is False, f"{code} should not be an addon"
    
    def test_tk_features_included_in_subscription(self):
        """TK features should be included by default in Standard+ tiers."""
        sub = CustomerSubscription(
            customer_id="test-001",
            tier=SubscriptionTier.STANDARD,
        )
        
        enabled = sub.get_enabled_features()
        
        # Basic TK features should be enabled
        assert "tk_attribution" in enabled
        assert "tk_sourcing" in enabled
        assert "tk_full_suite" in enabled
    
    def test_tk_opt_out_increases_price(self):
        """Opting out of TK should increase monthly price."""
        # With TK (default)
        sub_with_tk = CustomerSubscription(
            customer_id="test-with-tk",
            tier=SubscriptionTier.STANDARD,
            tk_opted_out=False,
        )
        
        # Without TK (opted out)
        sub_without_tk = CustomerSubscription(
            customer_id="test-without-tk",
            tier=SubscriptionTier.STANDARD,
            tk_opted_out=True,
            tk_opt_out_level="no_tk_full",
        )
        
        price_with_tk = sub_with_tk.calculate_monthly_price()
        price_without_tk = sub_without_tk.calculate_monthly_price()
        
        # Verify opt-out costs more
        assert price_without_tk > price_with_tk
        
        # Calculate expected surcharge
        surcharge = TK_OPT_OUT_SURCHARGES["no_tk_full"]
        base_price = TIER_DEFINITIONS[SubscriptionTier.STANDARD].base_price_monthly
        expected_surcharge = base_price * (surcharge.surcharge_percentage / 100) + surcharge.surcharge_flat
        
        assert price_without_tk == price_with_tk + expected_surcharge
    
    def test_surcharge_levels(self):
        """Test different surcharge levels for various opt-out types."""
        base_price = TIER_DEFINITIONS[SubscriptionTier.STANDARD].base_price_monthly
        
        # No attribution surcharge (15% + $150)
        no_attr = TK_OPT_OUT_SURCHARGES["no_tk_attribution"]
        attr_surcharge = base_price * 0.15 + 150.0
        expected_attr = base_price + attr_surcharge
        
        sub_no_attr = CustomerSubscription(
            customer_id="no-attr",
            tier=SubscriptionTier.STANDARD,
            tk_opted_out=True,
            tk_opt_out_level="no_tk_attribution",
        )
        assert sub_no_attr.calculate_monthly_price() == expected_attr
        
        # No compensation surcharge (25% + $250)
        no_comp = TK_OPT_OUT_SURCHARGES["no_tk_compensation"]
        comp_surcharge = base_price * 0.25 + 250.0
        expected_comp = base_price + comp_surcharge
        
        sub_no_comp = CustomerSubscription(
            customer_id="no-comp",
            tier=SubscriptionTier.STANDARD,
            tk_opted_out=True,
            tk_opt_out_level="no_tk_compensation",
        )
        assert sub_no_comp.calculate_monthly_price() == expected_comp
        
        # Full opt-out surcharge (40% + $500)
        no_full = TK_OPT_OUT_SURCHARGES["no_tk_full"]
        full_surcharge = base_price * 0.40 + 500.0
        expected_full = base_price + full_surcharge
        
        sub_no_full = CustomerSubscription(
            customer_id="no-full",
            tier=SubscriptionTier.STANDARD,
            tk_opted_out=True,
            tk_opt_out_level="no_tk_full",
        )
        assert sub_no_full.calculate_monthly_price() == expected_full


class TestTKOptOutWorkflow:
    """Tests for the TK opt-out workflow."""
    
    def setup_method(self):
        """Set up test fixtures."""
        self.manager = FeatureAccessManager()
        self.sub = CustomerSubscription(
            customer_id="workflow-test",
            tier=SubscriptionTier.STANDARD,
        )
        self.manager.register_subscription(self.sub)
    
    def test_request_tk_opt_out_requires_reason(self):
        """Opt-out should require a detailed reason."""
        result = self.manager.request_tk_opt_out(
            customer_id="workflow-test",
            opt_out_level="no_tk_full",
            reason="short"  # Too short
        )
        
        assert result["success"] is False
        assert "minimum 20 characters" in result["error"]
    
    def test_request_tk_opt_out_valid(self):
        """Valid opt-out request should return pricing info."""
        result = self.manager.request_tk_opt_out(
            customer_id="workflow-test",
            opt_out_level="no_tk_full",
            reason="Our company policy requires us to use proprietary-only knowledge sources.",
        )
        
        assert result["success"] is True
        assert result["confirmation_required"] is True
        assert result["price_increase"] > 0
        assert result["new_price"] > result["current_price"]
        assert "percentage_increase" in result  # The actual key name
        assert "flat_surcharge" in result
    
    def test_confirm_tk_opt_out(self):
        """Confirming opt-out should apply surcharge."""
        initial_price = self.sub.calculate_monthly_price()
        
        result = self.manager.confirm_tk_opt_out(
            customer_id="workflow-test",
            opt_out_level="no_tk_full",
            reason="Company policy requires proprietary-only sources.",
            confirmed=True,
        )
        
        assert result["success"] is True
        assert self.sub.tk_opted_out is True
        assert result["new_monthly_price"] > initial_price
        assert "tk_features_disabled" in result
    
    def test_confirm_tk_opt_out_not_confirmed(self):
        """Unconfirmed opt-out should fail."""
        result = self.manager.confirm_tk_opt_out(
            customer_id="workflow-test",
            opt_out_level="no_tk_full",
            reason="Some reason that is long enough for validation.",
            confirmed=False,
        )
        
        assert result["success"] is False
        assert self.sub.tk_opted_out is False


class TestTKReEnable:
    """Tests for re-enabling TK features after opt-out."""
    
    def setup_method(self):
        """Set up test fixtures with opted-out subscription."""
        self.manager = FeatureAccessManager()
        self.sub = CustomerSubscription(
            customer_id="re-enable-test",
            tier=SubscriptionTier.STANDARD,
            tk_opted_out=True,
            tk_opt_out_level="no_tk_full",
            tk_opt_out_reason="Test opt-out",
            tk_opt_out_acknowledged_at=datetime.utcnow(),
        )
        self.manager.register_subscription(self.sub)
    
    def test_enable_tk_features(self):
        """Re-enabling TK should remove surcharge."""
        opted_out_price = self.sub.calculate_monthly_price()
        
        result = self.manager.enable_tk_features("re-enable-test")
        
        assert result["success"] is True
        assert self.sub.tk_opted_out is False
        assert self.sub.tk_opt_out_level is None
        assert result["new_price"] < result["old_price"]
        assert result["savings"] == opted_out_price - result["new_price"]
        assert "tk_features_enabled" in result
    
    def test_enable_tk_features_already_enabled(self):
        """Re-enabling when already enabled should fail gracefully."""
        self.sub.tk_opted_out = False
        
        result = self.manager.enable_tk_features("re-enable-test")
        
        assert result["success"] is False
        assert "already enabled" in result["error"]


class TestFeatureAccessWithTK:
    """Tests for feature access with TK enabled/disabled."""
    
    def setup_method(self):
        """Set up test fixtures."""
        self.manager = FeatureAccessManager()
    
    def test_tk_features_accessible_by_default(self):
        """TK features should be accessible without opt-out."""
        sub = CustomerSubscription(
            customer_id="access-test",
            tier=SubscriptionTier.STANDARD,
        )
        self.manager.register_subscription(sub)
        
        result = self.manager.check_access("access-test", "tk_attribution")
        assert result["allowed"] is True
    
    def test_tk_features_blocked_after_opt_out(self):
        """TK features should be blocked after opt-out."""
        sub = CustomerSubscription(
            customer_id="blocked-test",
            tier=SubscriptionTier.STANDARD,
            tk_opted_out=True,
            tk_opt_out_level="no_tk_full",
        )
        self.manager.register_subscription(sub)
        
        result = self.manager.check_access("blocked-test", "tk_attribution")
        assert result["allowed"] is False
        assert "opt-out" in result["reason"].lower()
        assert "enable_tk_features" in result["upgrade_options"]


class TestPricingSummary:
    """Tests for the pricing summary with ethical pricing."""
    
    def setup_method(self):
        """Set up test fixtures."""
        self.manager = FeatureAccessManager()
    
    def test_pricing_summary_with_tk(self):
        """Pricing summary should show TK enabled."""
        sub = CustomerSubscription(
            customer_id="summary-test",
            tier=SubscriptionTier.STANDARD,
        )
        self.manager.register_subscription(sub)
        
        summary = self.manager.get_pricing_summary("summary-test")
        
        assert summary["tk_features_enabled"] is True
        assert summary["tk_opt_out_surcharge"] is None
        assert summary["ethical_discount_if_tk_enabled"] == 0
    
    def test_pricing_summary_without_tk(self):
        """Pricing summary should show surcharge details when opted out."""
        sub = CustomerSubscription(
            customer_id="surcharge-test",
            tier=SubscriptionTier.STANDARD,
            tk_opted_out=True,
            tk_opt_out_level="no_tk_full",
            tk_opt_out_reason="Test reason for compliance documentation",
            tk_opt_out_acknowledged_at=datetime.utcnow(),
        )
        self.manager.register_subscription(sub)
        
        summary = self.manager.get_pricing_summary("surcharge-test")
        
        assert summary["tk_features_enabled"] is False
        assert summary["tk_opt_out_surcharge"] is not None
        assert summary["tk_opt_out_surcharge"]["level"] == "no_tk_full"
        assert summary["tk_opt_out_surcharge"]["percentage"] == 40.0
        assert summary["ethical_discount_if_tk_enabled"] > 0


class TestTierBasedTKAvailability:
    """Tests for tier-based TK feature availability."""
    
    def test_basic_tier_features(self):
        """Basic tier should have core features."""
        sub = CustomerSubscription(
            customer_id="basic-test",
            tier=SubscriptionTier.BASIC,
        )
        
        enabled = sub.get_enabled_features()
        
        # Basic tier features should be present
        assert "api_access" in enabled
        assert "chempath_basic" in enabled
        
        # Note: TK features are included by default for ALL tiers (ethical model)
        # This is intentional - even basic users get TK features free
        assert "tk_attribution" in enabled
    
    def test_standard_tier_has_tk(self):
        """Standard tier should have basic TK features."""
        sub = CustomerSubscription(
            customer_id="standard-test",
            tier=SubscriptionTier.STANDARD,
        )
        
        enabled = sub.get_enabled_features()
        
        assert "tk_attribution" in enabled
        assert "tk_sourcing" in enabled
        assert "tk_full_suite" in enabled
    
    def test_premium_tier_has_equipath(self):
        """Premium tier should have EquiPath."""
        sub = CustomerSubscription(
            customer_id="premium-test",
            tier=SubscriptionTier.PREMIUM,
        )
        
        enabled = sub.get_enabled_features()
        
        assert "equipath" in enabled


class TestEthicalDiscountCalculation:
    """Tests for ethical discount/surcharge calculations."""
    
    def test_get_tk_discount_amount_when_opted_out(self):
        """Should return the potential savings from re-enabling TK."""
        sub = CustomerSubscription(
            customer_id="discount-test",
            tier=SubscriptionTier.STANDARD,
            tk_opted_out=True,
            tk_opt_out_level="no_tk_full",
        )
        
        discount = sub.get_tk_discount_amount()
        
        surcharge = TK_OPT_OUT_SURCHARGES["no_tk_full"]
        base_price = TIER_DEFINITIONS[SubscriptionTier.STANDARD].base_price_monthly
        expected_discount = base_price * (surcharge.surcharge_percentage / 100) + surcharge.surcharge_flat
        
        # Discount should equal the surcharge they're paying
        assert discount == expected_discount
        assert discount > 0
    
    def test_get_tk_discount_amount_when_enabled(self):
        """Should return 0 when TK is already enabled (no savings possible)."""
        sub = CustomerSubscription(
            customer_id="enabled-test",
            tier=SubscriptionTier.STANDARD,
            tk_opted_out=False,
        )
        
        discount = sub.get_tk_discount_amount()
        # Already saving - no additional discount available
        assert discount == 0.0


class TestPricingByTier:
    """Comprehensive pricing tests by tier."""
    
    @pytest.mark.parametrize("tier,base_price", [
        (SubscriptionTier.BASIC, 99.0),
        (SubscriptionTier.STANDARD, 299.0),
        (SubscriptionTier.PREMIUM, 799.0),
    ])
    def test_base_prices(self, tier, base_price):
        """Verify base prices match expected values."""
        tier_def = TIER_DEFINITIONS[tier]
        assert tier_def.base_price_monthly == base_price
    
    @pytest.mark.parametrize("tier,opt_out_level,expected_surcharge_pct,expected_surcharge_flat", [
        (SubscriptionTier.STANDARD, "no_tk_attribution", 15.0, 150.0),
        (SubscriptionTier.STANDARD, "no_tk_compensation", 25.0, 250.0),
        (SubscriptionTier.STANDARD, "no_tk_full", 40.0, 500.0),
        (SubscriptionTier.PREMIUM, "no_tk_full", 40.0, 500.0),
    ])
    def test_surcharge_calculations(self, tier, opt_out_level, expected_surcharge_pct, expected_surcharge_flat):
        """Verify surcharges are calculated correctly."""
        sub = CustomerSubscription(
            customer_id=f"calc-test-{tier.value}-{opt_out_level}",
            tier=tier,
            tk_opted_out=True,
            tk_opt_out_level=opt_out_level,
        )
        
        base_price = TIER_DEFINITIONS[tier].base_price_monthly
        expected_total = base_price + (base_price * expected_surcharge_pct / 100) + expected_surcharge_flat
        
        assert sub.calculate_monthly_price() == expected_total


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
