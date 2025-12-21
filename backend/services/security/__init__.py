"""
NeuroBotanica Security Services - Week 12

Security hardening and compliance module providing:
- Rate limiting with token bucket algorithm
- API key management with tiered access
- Comprehensive audit logging
- HIPAA/SOC2 compliance helpers
- Request fingerprinting and abuse detection
- Security event monitoring

Trade Secret Note:
- Rate limit thresholds and algorithms are protected
- Abuse detection heuristics are proprietary
- Fingerprinting methods are internal-only

Reference: NeuroBotanica MVP Development Plan - Week 12: Security hardening + compliance
"""

from .rate_limiter import (
    RateLimiter,
    RateLimitConfig,
    RateLimitTier,
    RateLimitResult,
    RateLimitExceeded,
    TokenBucket,
    SlidingWindowCounter,
    get_rate_limiter,
)

from .api_key_manager import (
    APIKeyManager,
    APIKey,
    APIKeyTier,
    APIKeyScope,
    APIKeyRequest,
    APIKeyResponse,
    get_api_key_manager,
)

from .audit_logger import (
    AuditLogger,
    AuditEvent,
    AuditEventType,
    AuditConfig,
    ComplianceLevel,
    SecurityEvent,
    get_audit_logger,
)

from .security_middleware import (
    SecurityMiddleware,
    SecurityConfig,
    RequestFingerprint,
    ThreatLevel,
)

__all__ = [
    # Rate Limiter
    "RateLimiter",
    "RateLimitConfig",
    "RateLimitTier",
    "RateLimitResult",
    "RateLimitExceeded",
    "TokenBucket",
    "SlidingWindowCounter",
    "get_rate_limiter",
    # API Key Manager
    "APIKeyManager",
    "APIKey",
    "APIKeyTier",
    "APIKeyScope",
    "APIKeyRequest",
    "APIKeyResponse",
    "get_api_key_manager",
    # Audit Logger
    "AuditLogger",
    "AuditEvent",
    "AuditEventType",
    "AuditConfig",
    "ComplianceLevel",
    "SecurityEvent",
    "get_audit_logger",
    # Security Middleware
    "SecurityMiddleware",
    "SecurityConfig",
    "RequestFingerprint",
    "ThreatLevel",
]
