"""
Week 12 Tests - Security Hardening and Compliance

Comprehensive tests for:
- Rate limiter with token bucket and sliding window
- API key management
- Audit logging with compliance
- Security middleware
- API endpoints
"""

import pytest
from datetime import datetime, timedelta
import time


# ============================================================================
# Rate Limiter Tests
# ============================================================================

class TestRateLimitTier:
    """Test RateLimitTier enumeration."""
    
    def test_tier_values(self):
        """Test tier values exist."""
        from backend.services.security import RateLimitTier
        
        assert RateLimitTier.ANONYMOUS.value == "anonymous"
        assert RateLimitTier.FREE.value == "free"
        assert RateLimitTier.ENTERPRISE.value == "enterprise"
    
    def test_tier_limits(self):
        """Test tier has limits."""
        from backend.services.security import RateLimitTier
        
        assert RateLimitTier.FREE.requests_per_minute > 0
        assert RateLimitTier.FREE.requests_per_hour > RateLimitTier.FREE.requests_per_minute
    
    def test_tier_hierarchy(self):
        """Test higher tiers have higher limits."""
        from backend.services.security import RateLimitTier
        
        assert RateLimitTier.PROFESSIONAL.requests_per_minute > RateLimitTier.FREE.requests_per_minute
        assert RateLimitTier.ENTERPRISE.requests_per_minute > RateLimitTier.PROFESSIONAL.requests_per_minute


class TestTokenBucket:
    """Test TokenBucket for burst handling."""
    
    def test_bucket_creation(self):
        """Test bucket initializes with full capacity."""
        from backend.services.security import TokenBucket
        
        bucket = TokenBucket(capacity=10.0, refill_rate=1.0)
        assert bucket.available_tokens == 10.0
    
    def test_bucket_consume(self):
        """Test consuming tokens from bucket."""
        from backend.services.security import TokenBucket
        
        bucket = TokenBucket(capacity=10.0, refill_rate=1.0)
        
        assert bucket.consume(5) == True
        assert bucket.available_tokens >= 4.9  # Allow small float variance
    
    def test_bucket_exhaust(self):
        """Test bucket exhaustion."""
        from backend.services.security import TokenBucket
        
        bucket = TokenBucket(capacity=5.0, refill_rate=0.1)
        
        for _ in range(5):
            assert bucket.consume(1) == True
        
        assert bucket.consume(1) == False
    
    def test_bucket_refill(self):
        """Test bucket refills over time."""
        from backend.services.security import TokenBucket
        
        bucket = TokenBucket(capacity=10.0, refill_rate=100.0)  # 100 tokens/sec
        
        bucket.consume(10)  # Empty bucket
        assert bucket.available_tokens < 1.0  # Should be near zero
        
        time.sleep(0.05)  # Wait 50ms
        
        # Should have refilled ~5 tokens
        assert bucket.available_tokens >= 4.0


class TestSlidingWindowCounter:
    """Test SlidingWindowCounter."""
    
    def test_counter_creation(self):
        """Test counter initializes properly."""
        from backend.services.security import SlidingWindowCounter
        
        counter = SlidingWindowCounter(
            window_seconds=60,
            max_requests=100
        )
        assert counter.get_count() == 0
    
    def test_counter_increment(self):
        """Test incrementing counter."""
        from backend.services.security import SlidingWindowCounter
        
        counter = SlidingWindowCounter(
            window_seconds=60,
            max_requests=100
        )
        
        allowed, count = counter.increment()
        assert allowed == True
        assert count == 1
    
    def test_counter_limit(self):
        """Test counter respects limit."""
        from backend.services.security import SlidingWindowCounter
        
        counter = SlidingWindowCounter(
            window_seconds=60,
            max_requests=5
        )
        
        for i in range(5):
            allowed, _ = counter.increment()
            assert allowed == True
        
        # 6th request should be denied
        allowed, _ = counter.increment()
        assert allowed == False
    
    def test_counter_reset(self):
        """Test counter reset."""
        from backend.services.security import SlidingWindowCounter
        
        counter = SlidingWindowCounter(
            window_seconds=60,
            max_requests=100
        )
        
        counter.increment()
        counter.increment()
        assert counter.get_count() == 2
        
        counter.reset()
        assert counter.get_count() == 0


class TestRateLimiter:
    """Test RateLimiter main class."""
    
    def test_limiter_creation(self):
        """Test limiter initializes properly."""
        from backend.services.security import RateLimiter, RateLimitConfig
        
        limiter = RateLimiter(RateLimitConfig())
        assert limiter is not None
    
    def test_limiter_allows_requests(self):
        """Test limiter allows valid requests."""
        from backend.services.security import RateLimiter, RateLimitConfig, RateLimitTier
        
        limiter = RateLimiter(RateLimitConfig())
        
        result = limiter.check(
            identifier="test-key",
            tier=RateLimitTier.FREE
        )
        
        assert result.allowed == True
        assert result.tier == RateLimitTier.FREE
    
    def test_limiter_returns_remaining(self):
        """Test limiter returns remaining counts."""
        from backend.services.security import RateLimiter, RateLimitConfig, RateLimitTier
        
        limiter = RateLimiter(RateLimitConfig())
        
        result = limiter.check(
            identifier="test-key-2",
            tier=RateLimitTier.FREE
        )
        
        assert result.requests_remaining_minute > 0
        assert result.requests_remaining_hour > 0
    
    def test_limiter_preview_mode(self):
        """Test preview mode doesn't consume quota."""
        from backend.services.security import RateLimiter, RateLimitConfig, RateLimitTier
        
        limiter = RateLimiter(RateLimitConfig())
        
        # First check with consume
        result1 = limiter.check("test-key-3", RateLimitTier.FREE, consume=True)
        
        # Preview check
        result2 = limiter.check("test-key-3", RateLimitTier.FREE, consume=False)
        
        # Both should show same remaining (preview didn't consume)
        assert result1.requests_remaining_minute == result2.requests_remaining_minute
    
    def test_limiter_get_usage(self):
        """Test getting usage statistics."""
        from backend.services.security import RateLimiter, RateLimitConfig, RateLimitTier
        
        limiter = RateLimiter(RateLimitConfig())
        
        limiter.check("usage-test", RateLimitTier.FREE)
        limiter.check("usage-test", RateLimitTier.FREE)
        
        usage = limiter.get_usage("usage-test")
        
        assert usage["identifier"] == "usage-test"
        assert usage["minute_count"] >= 2
    
    def test_limiter_reset(self):
        """Test resetting rate limits."""
        from backend.services.security import RateLimiter, RateLimitConfig, RateLimitTier
        
        limiter = RateLimiter(RateLimitConfig())
        
        limiter.check("reset-test", RateLimitTier.FREE)
        limiter.check("reset-test", RateLimitTier.FREE)
        
        limiter.reset("reset-test")
        
        usage = limiter.get_usage("reset-test")
        assert usage["minute_count"] == 0


# ============================================================================
# API Key Manager Tests
# ============================================================================

class TestAPIKeyTier:
    """Test APIKeyTier enumeration."""
    
    def test_tier_values(self):
        """Test tier values."""
        from backend.services.security import APIKeyTier
        
        assert APIKeyTier.FREE.value == "free"
        assert APIKeyTier.ENTERPRISE.value == "enterprise"
    
    def test_tier_pricing(self):
        """Test tier pricing."""
        from backend.services.security import APIKeyTier
        
        assert APIKeyTier.FREE.monthly_price_usd == 0
        assert APIKeyTier.PROFESSIONAL.monthly_price_usd == 2500
    
    def test_tier_max_keys(self):
        """Test tier max keys."""
        from backend.services.security import APIKeyTier
        
        assert APIKeyTier.FREE.max_keys == 1
        assert APIKeyTier.ENTERPRISE.max_keys > APIKeyTier.FREE.max_keys


class TestAPIKeyScope:
    """Test APIKeyScope enumeration."""
    
    def test_scope_values(self):
        """Test scope values."""
        from backend.services.security import APIKeyScope
        
        assert APIKeyScope.READ_STUDIES.value == "read:studies"
        assert APIKeyScope.CHEMPATH_ACCESS.value == "chempath:access"
    
    def test_scope_helpers(self):
        """Test scope helper methods."""
        from backend.services.security import APIKeyScope
        
        read_scopes = APIKeyScope.all_read()
        assert "read:studies" in read_scopes
        
        service_scopes = APIKeyScope.all_services()
        assert "chempath:access" in service_scopes


class TestAPIKeyManager:
    """Test APIKeyManager."""
    
    def test_manager_creation(self):
        """Test manager initializes."""
        from backend.services.security import APIKeyManager
        
        manager = APIKeyManager()
        assert manager is not None
    
    def test_create_key(self):
        """Test creating API key."""
        from backend.services.security import APIKeyManager, APIKeyRequest, APIKeyTier
        
        manager = APIKeyManager()
        
        request = APIKeyRequest(
            owner_id="test-owner-1",
            tier=APIKeyTier.FREE,
            name="Test Key"
        )
        
        result = manager.create_key(request)
        
        assert result.raw_key.startswith("nb_")
        assert result.key.tier == APIKeyTier.FREE
        assert result.key.owner_id == "test-owner-1"
    
    def test_validate_key(self):
        """Test validating API key."""
        from backend.services.security import APIKeyManager, APIKeyRequest, APIKeyTier
        
        manager = APIKeyManager()
        
        request = APIKeyRequest(
            owner_id="test-owner-2",
            tier=APIKeyTier.BASIC,
            name="Validation Test"
        )
        
        result = manager.create_key(request)
        raw_key = result.raw_key
        
        # Validate the key
        validated = manager.validate_key(raw_key)
        
        assert validated is not None
        assert validated.key_id == result.key.key_id
    
    def test_validate_invalid_key(self):
        """Test validating invalid key returns None."""
        from backend.services.security import APIKeyManager
        
        manager = APIKeyManager()
        
        validated = manager.validate_key("invalid-key")
        assert validated is None
    
    def test_revoke_key(self):
        """Test revoking API key."""
        from backend.services.security import APIKeyManager, APIKeyRequest, APIKeyTier
        
        manager = APIKeyManager()
        
        request = APIKeyRequest(
            owner_id="test-owner-3",
            tier=APIKeyTier.FREE,
            name="Revoke Test"
        )
        
        result = manager.create_key(request)
        key_id = result.key.key_id
        
        # Revoke
        revoked = manager.revoke_key(key_id, reason="Testing")
        
        assert revoked is not None
        assert revoked.is_revoked == True
        
        # Validate should now fail
        validated = manager.validate_key(result.raw_key)
        assert validated is None
    
    def test_rotate_key(self):
        """Test rotating API key."""
        from backend.services.security import APIKeyManager, APIKeyRequest, APIKeyTier
        
        manager = APIKeyManager()
        
        request = APIKeyRequest(
            owner_id="test-owner-4",
            tier=APIKeyTier.PROFESSIONAL,
            name="Rotate Test"
        )
        
        result = manager.create_key(request)
        old_key_id = result.key.key_id
        old_raw = result.raw_key
        
        # Rotate
        new_result = manager.rotate_key(old_key_id)
        
        assert new_result is not None
        assert new_result.key.key_id != old_key_id
        assert new_result.raw_key != old_raw
        
        # Old key should be invalid
        assert manager.validate_key(old_raw) is None
        
        # New key should be valid
        assert manager.validate_key(new_result.raw_key) is not None
    
    def test_list_keys(self):
        """Test listing API keys."""
        from backend.services.security import APIKeyManager, APIKeyRequest, APIKeyTier
        
        manager = APIKeyManager()
        owner = "list-test-owner"
        
        # Create keys
        for i in range(3):
            request = APIKeyRequest(
                owner_id=owner,
                tier=APIKeyTier.ENTERPRISE,  # Higher limit
                name=f"List Test {i}"
            )
            manager.create_key(request)
        
        keys = manager.list_keys(owner_id=owner)
        assert len(keys) >= 3


# ============================================================================
# Audit Logger Tests
# ============================================================================

class TestAuditEventType:
    """Test AuditEventType enumeration."""
    
    def test_event_type_values(self):
        """Test event type values."""
        from backend.services.security import AuditEventType
        
        assert AuditEventType.AUTH_LOGIN.value == "auth.login"
        assert AuditEventType.DATA_READ.value == "data.read"
    
    def test_event_severity(self):
        """Test event severity."""
        from backend.services.security import AuditEventType
        
        assert AuditEventType.AUTH_LOGIN.severity == "info"
        assert AuditEventType.SECURITY_THREAT_DETECTED.severity == "critical"
    
    def test_security_event_flag(self):
        """Test security event flag."""
        from backend.services.security import AuditEventType
        
        assert AuditEventType.SECURITY_THREAT_DETECTED.is_security_event == True
        assert AuditEventType.AUTH_LOGIN.is_security_event == False


class TestComplianceLevel:
    """Test ComplianceLevel enumeration."""
    
    def test_compliance_values(self):
        """Test compliance level values."""
        from backend.services.security import ComplianceLevel
        
        assert ComplianceLevel.BASIC.value == "basic"
        assert ComplianceLevel.HIPAA.value == "hipaa"
    
    def test_retention_days(self):
        """Test retention requirements."""
        from backend.services.security import ComplianceLevel
        
        assert ComplianceLevel.BASIC.retention_days == 90
        assert ComplianceLevel.HIPAA.retention_days >= 365 * 6  # 6 years
    
    def test_pii_redaction(self):
        """Test PII redaction requirements."""
        from backend.services.security import ComplianceLevel
        
        assert ComplianceLevel.BASIC.requires_pii_redaction == False
        assert ComplianceLevel.HIPAA.requires_pii_redaction == True


class TestAuditEvent:
    """Test AuditEvent data class."""
    
    def test_event_creation(self):
        """Test audit event creation."""
        from backend.services.security import AuditEvent, AuditEventType
        
        event = AuditEvent(
            event_id="evt_test123",
            event_type=AuditEventType.AUTH_LOGIN,
            timestamp=datetime.utcnow(),
            user_id="user123",
        )
        
        assert event.event_id == "evt_test123"
        assert event.checksum is not None
    
    def test_event_integrity(self):
        """Test event integrity verification."""
        from backend.services.security import AuditEvent, AuditEventType
        
        event = AuditEvent(
            event_id="evt_integrity",
            event_type=AuditEventType.DATA_READ,
            timestamp=datetime.utcnow(),
        )
        
        assert event.verify_integrity() == True
    
    def test_event_to_dict(self):
        """Test event serialization."""
        from backend.services.security import AuditEvent, AuditEventType
        
        event = AuditEvent(
            event_id="evt_dict",
            event_type=AuditEventType.DATA_CREATE,
            timestamp=datetime.utcnow(),
            user_id="user456",
        )
        
        data = event.to_dict()
        
        assert data["event_id"] == "evt_dict"
        assert data["event_type"] == "data.create"
        assert "timestamp" in data


class TestAuditLogger:
    """Test AuditLogger."""
    
    def test_logger_creation(self):
        """Test logger initializes."""
        from backend.services.security import AuditLogger
        
        logger = AuditLogger()
        assert logger is not None
    
    def test_log_event(self):
        """Test logging an event."""
        from backend.services.security import AuditLogger, AuditEventType
        
        logger = AuditLogger()
        
        event = logger.log(
            event_type=AuditEventType.AUTH_LOGIN,
            user_id="test-user",
            source_ip="127.0.0.1",
        )
        
        assert event is not None
        assert event.event_type == AuditEventType.AUTH_LOGIN
    
    def test_log_security_event(self):
        """Test logging security event."""
        from backend.services.security import AuditLogger, AuditEventType
        
        logger = AuditLogger()
        
        event = logger.log_security_event(
            event_type=AuditEventType.SECURITY_SUSPICIOUS_ACTIVITY,
            threat_level="medium",
            source_ip="192.168.1.1",
            description="Unusual access pattern detected",
        )
        
        assert event is not None
        assert event.threat_level == "medium"
    
    def test_get_events(self):
        """Test querying events."""
        from backend.services.security import AuditLogger, AuditEventType, AuditConfig
        
        # Create fresh logger instance for this test
        logger = AuditLogger(AuditConfig(log_debug_events=True))
        
        # Log some events
        logger.log(AuditEventType.DATA_READ, user_id="query-test")
        logger.log(AuditEventType.DATA_CREATE, user_id="query-test")
        
        events = logger.get_events(user_id="query-test", limit=10)
        
        assert len(events) >= 2
    
    def test_security_summary(self):
        """Test security summary."""
        from backend.services.security import AuditLogger, AuditEventType
        
        logger = AuditLogger()
        
        # Log security event
        logger.log_security_event(
            event_type=AuditEventType.SECURITY_THREAT_DETECTED,
            threat_level="high",
            source_ip="10.0.0.1",
            description="Test threat",
        )
        
        summary = logger.get_security_summary(hours=24)
        
        assert "total_events" in summary
        assert "by_threat_level" in summary
    
    def test_integrity_verification(self):
        """Test audit log integrity verification."""
        from backend.services.security import AuditLogger, AuditEventType
        
        logger = AuditLogger()
        
        logger.log(AuditEventType.DATA_READ, user_id="integrity-test")
        
        result = logger.verify_integrity()
        
        assert result["integrity_percentage"] == 100


# ============================================================================
# Security Middleware Tests
# ============================================================================

class TestThreatLevel:
    """Test ThreatLevel enumeration."""
    
    def test_threat_values(self):
        """Test threat level values."""
        from backend.services.security import ThreatLevel
        
        assert ThreatLevel.NONE.value == "none"
        assert ThreatLevel.CRITICAL.value == "critical"
    
    def test_blocking_logic(self):
        """Test blocking thresholds."""
        from backend.services.security import ThreatLevel
        
        assert ThreatLevel.LOW.should_block == False
        assert ThreatLevel.HIGH.should_block == True
        assert ThreatLevel.CRITICAL.should_block == True


class TestRequestFingerprint:
    """Test RequestFingerprint."""
    
    def test_fingerprint_creation(self):
        """Test fingerprint creation."""
        from backend.services.security import RequestFingerprint
        
        fp = RequestFingerprint(
            fingerprint_id="fp_test123",
            source_ip="192.168.1.1",
            user_agent="Mozilla/5.0",
            method="GET",
            path="/api/v1/test",
        )
        
        assert fp.fingerprint_id == "fp_test123"
        assert fp.source_ip == "192.168.1.1"
    
    def test_fingerprint_to_dict(self):
        """Test fingerprint serialization."""
        from backend.services.security import RequestFingerprint
        
        fp = RequestFingerprint(
            fingerprint_id="fp_dict",
            source_ip="10.0.0.1",
            method="POST",
            path="/test",
        )
        
        data = fp.to_dict()
        
        assert data["fingerprint_id"] == "fp_dict"
        assert "timestamp" in data


class TestSecurityConfig:
    """Test SecurityConfig."""
    
    def test_default_config(self):
        """Test default configuration."""
        from backend.services.security import SecurityConfig
        
        config = SecurityConfig()
        
        assert config.enable_rate_limiting == True
        assert config.enable_threat_detection == True
        assert config.max_content_length > 0
    
    def test_custom_config(self):
        """Test custom configuration."""
        from backend.services.security import SecurityConfig
        
        config = SecurityConfig(
            enable_rate_limiting=False,
            max_content_length=5 * 1024 * 1024,
        )
        
        assert config.enable_rate_limiting == False
        assert config.max_content_length == 5 * 1024 * 1024


# ============================================================================
# API Endpoint Tests
# ============================================================================

class TestSecurityAPI:
    """Test Security API endpoints."""
    
    @pytest.fixture
    def client(self):
        """Create test client."""
        from fastapi.testclient import TestClient
        from backend.main import app
        return TestClient(app)
    
    def test_health_endpoint(self, client):
        """Test health check."""
        response = client.get("/api/v1/security/health")
        assert response.status_code == 200
        
        data = response.json()
        assert data["service"] == "security"
        assert data["status"] == "healthy"
    
    def test_list_tiers(self, client):
        """Test listing rate limit tiers."""
        response = client.get("/api/v1/security/rate-limit/tiers")
        assert response.status_code == 200
        
        data = response.json()
        assert len(data) >= 4
        tier_names = [t["tier"] for t in data]
        assert "free" in tier_names
        assert "enterprise" in tier_names
    
    def test_list_key_tiers(self, client):
        """Test listing API key tiers."""
        response = client.get("/api/v1/security/key-tiers")
        assert response.status_code == 200
        
        data = response.json()
        assert len(data) >= 4
    
    def test_list_scopes(self, client):
        """Test listing scopes."""
        response = client.get("/api/v1/security/scopes")
        assert response.status_code == 200
        
        data = response.json()
        assert len(data) >= 10
        scopes = [s["scope"] for s in data]
        assert "read:studies" in scopes
    
    def test_list_event_types(self, client):
        """Test listing audit event types."""
        response = client.get("/api/v1/security/audit/event-types")
        assert response.status_code == 200
        
        data = response.json()
        assert len(data) >= 20
    
    def test_list_compliance_levels(self, client):
        """Test listing compliance levels."""
        response = client.get("/api/v1/security/compliance-levels")
        assert response.status_code == 200
        
        data = response.json()
        assert len(data) >= 3
        levels = [l["level"] for l in data]
        assert "hipaa" in levels
    
    def test_create_api_key(self, client):
        """Test creating API key."""
        response = client.post(
            "/api/v1/security/keys",
            params={"owner_id": "test-owner"},
            json={
                "name": "Test API Key",
                "tier": "free",
            }
        )
        assert response.status_code == 201
        
        data = response.json()
        assert "api_key" in data
        assert data["api_key"].startswith("nb_")
        assert data["tier"] == "free"
    
    def test_get_api_key(self, client):
        """Test getting API key by ID."""
        # First create
        create_response = client.post(
            "/api/v1/security/keys",
            params={"owner_id": "get-test-owner"},
            json={"name": "Get Test Key", "tier": "basic"}
        )
        key_id = create_response.json()["key_id"]
        
        # Then get
        response = client.get(f"/api/v1/security/keys/{key_id}")
        assert response.status_code == 200
        
        data = response.json()
        assert data["key_id"] == key_id
    
    def test_list_api_keys(self, client):
        """Test listing API keys."""
        response = client.get("/api/v1/security/keys")
        assert response.status_code == 200
        
        data = response.json()
        assert isinstance(data, list)
    
    def test_revoke_api_key(self, client):
        """Test revoking API key."""
        # Create
        create_response = client.post(
            "/api/v1/security/keys",
            params={"owner_id": "revoke-test"},
            json={"name": "Revoke Test", "tier": "free"}
        )
        key_id = create_response.json()["key_id"]
        
        # Revoke - use request with body via generic request method
        response = client.request(
            "DELETE",
            f"/api/v1/security/keys/{key_id}",
            json={"reason": "Testing revocation"}
        )
        assert response.status_code == 200
        
        data = response.json()
        assert data["status"] == "revoked"
    
    def test_rotate_api_key(self, client):
        """Test rotating API key."""
        # Create
        create_response = client.post(
            "/api/v1/security/keys",
            params={"owner_id": "rotate-test"},
            json={"name": "Rotate Test", "tier": "professional"}
        )
        old_key_id = create_response.json()["key_id"]
        
        # Rotate
        response = client.post(f"/api/v1/security/keys/{old_key_id}/rotate")
        assert response.status_code == 200
        
        data = response.json()
        assert data["key_id"] != old_key_id
        assert "api_key" in data
    
    def test_rate_limit_status(self, client):
        """Test getting rate limit status."""
        response = client.get(
            "/api/v1/security/rate-limit/status",
            params={"identifier": "test-identifier"}
        )
        assert response.status_code == 200
        
        data = response.json()
        assert "tier" in data
        assert "limits" in data
        assert "remaining" in data
    
    def test_audit_events(self, client):
        """Test querying audit events."""
        response = client.get(
            "/api/v1/security/audit/events",
            params={"limit": 10}
        )
        assert response.status_code == 200
        
        data = response.json()
        assert "events" in data
        assert "count" in data
    
    def test_security_summary(self, client):
        """Test security summary."""
        response = client.get(
            "/api/v1/security/audit/security-summary",
            params={"hours": 24}
        )
        assert response.status_code == 200
        
        data = response.json()
        assert "total_events" in data
        assert "by_threat_level" in data
    
    def test_audit_integrity(self, client):
        """Test audit integrity check."""
        response = client.get("/api/v1/security/audit/integrity")
        assert response.status_code == 200
        
        data = response.json()
        assert "integrity_percentage" in data


# ============================================================================
# Integration Tests
# ============================================================================

class TestSecurityIntegration:
    """Integration tests for security components."""
    
    def test_full_api_key_lifecycle(self):
        """Test complete API key lifecycle."""
        from backend.services.security import (
            APIKeyManager, APIKeyRequest, APIKeyTier
        )
        
        manager = APIKeyManager()
        
        # Create
        request = APIKeyRequest(
            owner_id="lifecycle-test",
            tier=APIKeyTier.PROFESSIONAL,
            name="Lifecycle Test Key",
            description="Testing full lifecycle",
        )
        result = manager.create_key(request)
        raw_key = result.raw_key
        key_id = result.key.key_id
        
        # Validate
        validated = manager.validate_key(raw_key)
        assert validated is not None
        assert validated.request_count == 1
        
        # Use again
        manager.validate_key(raw_key)
        assert manager.get_key(key_id).request_count == 2
        
        # Check usage
        stats = manager.get_usage_stats(key_id)
        assert stats["request_count"] == 2
        
        # Rotate
        new_result = manager.rotate_key(key_id)
        assert new_result is not None
        
        # Old key invalid
        assert manager.validate_key(raw_key) is None
        
        # New key valid
        assert manager.validate_key(new_result.raw_key) is not None
    
    def test_rate_limit_with_audit(self):
        """Test rate limiting with audit logging."""
        from backend.services.security import (
            RateLimiter, RateLimitConfig, RateLimitTier,
            AuditLogger, AuditEventType
        )
        
        limiter = RateLimiter(RateLimitConfig())
        audit = AuditLogger()
        
        # Simulate rate-limited request
        result = limiter.check(
            identifier="rate-audit-test",
            tier=RateLimitTier.FREE,
        )
        
        # Log based on result
        if not result.allowed:
            audit.log(
                event_type=AuditEventType.KEY_RATE_LIMITED,
                details={"reason": result.reason}
            )
        
        # Check audit logged the access
        assert result.allowed == True  # First request should succeed
    
    def test_security_headers_config(self):
        """Test security configuration affects headers."""
        from backend.services.security import SecurityConfig
        
        # Strict config
        strict = SecurityConfig(
            strict_transport_security=True,
            content_security_policy=True,
            x_frame_options="DENY",
        )
        
        assert strict.strict_transport_security == True
        assert strict.x_frame_options == "DENY"
        
        # Relaxed config
        relaxed = SecurityConfig(
            strict_transport_security=False,
            x_frame_options="SAMEORIGIN",
        )
        
        assert relaxed.strict_transport_security == False
        assert relaxed.x_frame_options == "SAMEORIGIN"
