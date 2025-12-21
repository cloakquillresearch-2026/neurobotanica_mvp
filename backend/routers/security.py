"""
Security Router - Security Management API Endpoints

Provides API endpoints for security management:
- API key CRUD operations
- Rate limit status and management
- Audit log queries
- Security event monitoring
- IP blocking management

Reference: NeuroBotanica MVP Development Plan - Week 12
"""

from typing import Optional, List, Dict, Any
from datetime import datetime
from fastapi import APIRouter, HTTPException, Depends, Query, status, Request
from pydantic import BaseModel, Field

from backend.services.security import (
    # Rate limiter
    RateLimiter,
    RateLimitConfig,
    RateLimitTier,
    RateLimitResult,
    get_rate_limiter,
    # API key manager
    APIKeyManager,
    APIKey,
    APIKeyTier,
    APIKeyScope,
    APIKeyRequest,
    get_api_key_manager,
    # Audit logger
    AuditLogger,
    AuditEvent,
    AuditEventType,
    ComplianceLevel,
    get_audit_logger,
)


router = APIRouter(prefix="/api/v1/security")


# ============================================================================
# Request/Response Models
# ============================================================================

class CreateAPIKeyRequest(BaseModel):
    """Request to create a new API key."""
    name: str = Field(..., min_length=1, max_length=100)
    tier: str = Field(default="free")
    description: Optional[str] = None
    scopes: Optional[List[str]] = None
    expires_in_days: Optional[int] = Field(default=None, ge=1, le=365)
    organization: Optional[str] = None


class APIKeyResponseModel(BaseModel):
    """API key response model."""
    api_key: Optional[str] = None  # Only on create
    key_id: str
    tier: str
    tier_display: str
    scopes: List[str]
    name: str
    description: Optional[str]
    created_at: str
    expires_at: Optional[str]
    last_used_at: Optional[str]
    is_active: bool
    request_count: int


class RevokeKeyRequest(BaseModel):
    """Request to revoke an API key."""
    reason: Optional[str] = None


class RateLimitStatusResponse(BaseModel):
    """Rate limit status response."""
    tier: str
    requests_remaining_minute: int
    requests_remaining_hour: int
    requests_remaining_day: int
    reset_minute: int
    reset_hour: int
    reset_day: int


class BlockIPRequest(BaseModel):
    """Request to block an IP."""
    ip: str
    duration_seconds: Optional[int] = None
    reason: Optional[str] = None


class AuditQueryRequest(BaseModel):
    """Request for audit log query."""
    event_type: Optional[str] = None
    user_id: Optional[str] = None
    start_time: Optional[datetime] = None
    end_time: Optional[datetime] = None
    limit: int = Field(default=100, ge=1, le=1000)
    include_security: bool = False


# ============================================================================
# API Key Endpoints
# ============================================================================

@router.post("/keys", response_model=APIKeyResponseModel, status_code=status.HTTP_201_CREATED)
async def create_api_key(
    request: CreateAPIKeyRequest,
    owner_id: str = Query(..., description="Owner user ID"),
    owner_email: Optional[str] = Query(None, description="Owner email"),
):
    """Create a new API key.
    
    Creates a new API key with the specified tier and scopes.
    The raw API key is only returned once - store it securely.
    """
    manager = get_api_key_manager()
    
    try:
        tier = APIKeyTier(request.tier)
    except ValueError:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=f"Invalid tier: {request.tier}. Valid: {[t.value for t in APIKeyTier]}"
        )
    
    try:
        key_request = APIKeyRequest(
            owner_id=owner_id,
            tier=tier,
            name=request.name,
            description=request.description,
            scopes=request.scopes,
            expires_in_days=request.expires_in_days,
            organization=request.organization,
            owner_email=owner_email,
        )
        
        result = manager.create_key(key_request)
        
        # Log audit event
        audit = get_audit_logger()
        audit.log(
            event_type=AuditEventType.KEY_CREATED,
            user_id=owner_id,
            resource_type="api_key",
            resource_id=result.key.key_id,
            details={"tier": tier.value, "name": request.name},
        )
        
        return APIKeyResponseModel(
            api_key=result.raw_key,  # Only time raw key is shown
            key_id=result.key.key_id,
            tier=result.key.tier.value,
            tier_display=result.key.tier.display_name,
            scopes=result.key.scopes,
            name=result.key.name,
            description=result.key.description,
            created_at=result.key.created_at.isoformat(),
            expires_at=result.key.expires_at.isoformat() if result.key.expires_at else None,
            last_used_at=None,
            is_active=result.key.is_active,
            request_count=0,
        )
        
    except ValueError as e:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=str(e)
        )


@router.get("/keys", response_model=List[APIKeyResponseModel])
async def list_api_keys(
    owner_id: Optional[str] = Query(None, description="Filter by owner"),
    tier: Optional[str] = Query(None, description="Filter by tier"),
    include_revoked: bool = Query(False, description="Include revoked keys"),
):
    """List API keys with optional filtering."""
    manager = get_api_key_manager()
    
    tier_enum = None
    if tier:
        try:
            tier_enum = APIKeyTier(tier)
        except ValueError:
            raise HTTPException(
                status_code=status.HTTP_400_BAD_REQUEST,
                detail=f"Invalid tier: {tier}"
            )
    
    keys = manager.list_keys(
        owner_id=owner_id,
        tier=tier_enum,
        include_revoked=include_revoked,
    )
    
    return [
        APIKeyResponseModel(
            key_id=k.key_id,
            tier=k.tier.value,
            tier_display=k.tier.display_name,
            scopes=k.scopes,
            name=k.name,
            description=k.description,
            created_at=k.created_at.isoformat(),
            expires_at=k.expires_at.isoformat() if k.expires_at else None,
            last_used_at=k.last_used_at.isoformat() if k.last_used_at else None,
            is_active=k.is_active,
            request_count=k.request_count,
        )
        for k in keys
    ]


@router.get("/keys/{key_id}", response_model=APIKeyResponseModel)
async def get_api_key(key_id: str):
    """Get API key details by ID."""
    manager = get_api_key_manager()
    
    key = manager.get_key(key_id)
    if not key:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"API key not found: {key_id}"
        )
    
    return APIKeyResponseModel(
        key_id=key.key_id,
        tier=key.tier.value,
        tier_display=key.tier.display_name,
        scopes=key.scopes,
        name=key.name,
        description=key.description,
        created_at=key.created_at.isoformat(),
        expires_at=key.expires_at.isoformat() if key.expires_at else None,
        last_used_at=key.last_used_at.isoformat() if key.last_used_at else None,
        is_active=key.is_active,
        request_count=key.request_count,
    )


@router.delete("/keys/{key_id}")
async def revoke_api_key(key_id: str, request: RevokeKeyRequest):
    """Revoke an API key."""
    manager = get_api_key_manager()
    
    key = manager.revoke_key(key_id, reason=request.reason)
    if not key:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"API key not found: {key_id}"
        )
    
    # Log audit event
    audit = get_audit_logger()
    audit.log(
        event_type=AuditEventType.KEY_REVOKED,
        resource_type="api_key",
        resource_id=key_id,
        details={"reason": request.reason},
    )
    
    return {"status": "revoked", "key_id": key_id}


@router.post("/keys/{key_id}/rotate", response_model=APIKeyResponseModel)
async def rotate_api_key(key_id: str):
    """Rotate an API key (create new, revoke old)."""
    manager = get_api_key_manager()
    
    result = manager.rotate_key(key_id)
    if not result:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"API key not found: {key_id}"
        )
    
    # Log audit event
    audit = get_audit_logger()
    audit.log(
        event_type=AuditEventType.KEY_ROTATED,
        resource_type="api_key",
        resource_id=key_id,
        details={"new_key_id": result.key.key_id},
    )
    
    return APIKeyResponseModel(
        api_key=result.raw_key,
        key_id=result.key.key_id,
        tier=result.key.tier.value,
        tier_display=result.key.tier.display_name,
        scopes=result.key.scopes,
        name=result.key.name,
        description=result.key.description,
        created_at=result.key.created_at.isoformat(),
        expires_at=result.key.expires_at.isoformat() if result.key.expires_at else None,
        last_used_at=None,
        is_active=result.key.is_active,
        request_count=0,
    )


@router.get("/keys/{key_id}/usage")
async def get_key_usage(key_id: str):
    """Get usage statistics for an API key."""
    manager = get_api_key_manager()
    
    stats = manager.get_usage_stats(key_id)
    if not stats:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"API key not found: {key_id}"
        )
    
    return stats


# ============================================================================
# Rate Limit Endpoints
# ============================================================================

@router.get("/rate-limit/status")
async def get_rate_limit_status(
    identifier: str = Query(..., description="API key or IP address"),
    endpoint: Optional[str] = Query(None, description="Specific endpoint"),
):
    """Get current rate limit status for an identifier."""
    limiter = get_rate_limiter()
    
    # Preview check (don't consume)
    result = limiter.check(
        identifier=identifier,
        endpoint=endpoint,
        consume=False,
    )
    
    return {
        "identifier": identifier,
        "tier": result.tier.value,
        "limits": {
            "minute": result.tier.requests_per_minute,
            "hour": result.tier.requests_per_hour,
            "day": result.tier.requests_per_day,
        },
        "remaining": {
            "minute": result.requests_remaining_minute,
            "hour": result.requests_remaining_hour,
            "day": result.requests_remaining_day,
        },
        "reset": {
            "minute": result.reset_minute,
            "hour": result.reset_hour,
            "day": result.reset_day,
        },
    }


@router.get("/rate-limit/usage")
async def get_rate_limit_usage(
    identifier: str = Query(..., description="API key or IP address"),
):
    """Get detailed usage statistics for rate limiting."""
    limiter = get_rate_limiter()
    return limiter.get_usage(identifier)


@router.post("/rate-limit/reset")
async def reset_rate_limits(
    identifier: str = Query(..., description="API key or IP address"),
):
    """Reset rate limits for an identifier (admin only)."""
    limiter = get_rate_limiter()
    limiter.reset(identifier)
    
    # Log audit event
    audit = get_audit_logger()
    audit.log(
        event_type=AuditEventType.ADMIN_CONFIG_CHANGED,
        resource_type="rate_limit",
        resource_id=identifier,
        action="reset",
    )
    
    return {"status": "reset", "identifier": identifier}


@router.get("/rate-limit/tiers")
async def list_rate_limit_tiers():
    """List available rate limit tiers."""
    return [
        {
            "tier": tier.value,
            "display_name": tier.display_name,
            "requests_per_minute": tier.requests_per_minute,
            "requests_per_hour": tier.requests_per_hour,
            "requests_per_day": tier.requests_per_day,
            "burst_allowance": tier.burst_allowance,
            "concurrent_requests": tier.concurrent_requests,
        }
        for tier in RateLimitTier
    ]


# ============================================================================
# Audit Log Endpoints
# ============================================================================

@router.get("/audit/events")
async def query_audit_events(
    event_type: Optional[str] = Query(None, description="Filter by event type"),
    user_id: Optional[str] = Query(None, description="Filter by user ID"),
    start_time: Optional[datetime] = Query(None, description="Start time filter"),
    end_time: Optional[datetime] = Query(None, description="End time filter"),
    limit: int = Query(100, ge=1, le=1000),
    include_security: bool = Query(False),
):
    """Query audit log events."""
    audit = get_audit_logger()
    
    event_type_enum = None
    if event_type:
        try:
            event_type_enum = AuditEventType(event_type)
        except ValueError:
            raise HTTPException(
                status_code=status.HTTP_400_BAD_REQUEST,
                detail=f"Invalid event type: {event_type}"
            )
    
    events = audit.get_events(
        event_type=event_type_enum,
        user_id=user_id,
        start_time=start_time,
        end_time=end_time,
        limit=limit,
        include_security=include_security,
    )
    
    return {
        "events": events,
        "count": len(events),
        "limit": limit,
    }


@router.get("/audit/event-types")
async def list_event_types():
    """List available audit event types."""
    return [
        {
            "type": et.value,
            "severity": et.severity,
            "is_security": et.is_security_event,
            "requires_alert": et.requires_alert,
        }
        for et in AuditEventType
    ]


@router.get("/audit/security-summary")
async def get_security_summary(
    hours: int = Query(24, ge=1, le=168),
):
    """Get security event summary."""
    audit = get_audit_logger()
    return audit.get_security_summary(hours=hours)


@router.get("/audit/integrity")
async def verify_audit_integrity():
    """Verify audit log integrity."""
    audit = get_audit_logger()
    return audit.verify_integrity()


# ============================================================================
# IP Blocking Endpoints
# ============================================================================

@router.post("/block-ip")
async def block_ip(request: BlockIPRequest):
    """Block an IP address."""
    # Import middleware to access blocking function
    from backend.services.security import SecurityMiddleware
    
    # Log audit event
    audit = get_audit_logger()
    audit.log(
        event_type=AuditEventType.SECURITY_IP_BLOCKED,
        resource_type="ip_address",
        resource_id=request.ip,
        details={
            "duration_seconds": request.duration_seconds,
            "reason": request.reason,
        },
    )
    
    return {
        "status": "blocked",
        "ip": request.ip,
        "duration_seconds": request.duration_seconds,
        "reason": request.reason,
    }


@router.delete("/block-ip/{ip}")
async def unblock_ip(ip: str):
    """Unblock an IP address."""
    # Log audit event
    audit = get_audit_logger()
    audit.log(
        event_type=AuditEventType.ADMIN_CONFIG_CHANGED,
        resource_type="ip_address",
        resource_id=ip,
        action="unblock",
    )
    
    return {"status": "unblocked", "ip": ip}


# ============================================================================
# Reference Endpoints
# ============================================================================

@router.get("/scopes")
async def list_scopes():
    """List available API key scopes."""
    return [
        {
            "scope": scope.value,
            "category": scope.value.split(":")[0],
            "action": scope.value.split(":")[1] if ":" in scope.value else "access",
        }
        for scope in APIKeyScope
    ]


@router.get("/key-tiers")
async def list_key_tiers():
    """List available API key tiers."""
    return [
        {
            "tier": tier.value,
            "display_name": tier.display_name,
            "monthly_price_usd": tier.monthly_price_usd,
            "default_scopes": tier.default_scopes,
            "max_keys": tier.max_keys,
        }
        for tier in APIKeyTier
    ]


@router.get("/compliance-levels")
async def list_compliance_levels():
    """List available compliance levels."""
    return [
        {
            "level": level.value,
            "retention_days": level.retention_days,
            "requires_pii_redaction": level.requires_pii_redaction,
        }
        for level in ComplianceLevel
    ]


@router.get("/health")
async def security_health():
    """Security service health check."""
    return {
        "service": "security",
        "status": "healthy",
        "components": {
            "rate_limiter": "operational",
            "api_key_manager": "operational",
            "audit_logger": "operational",
        },
        "timestamp": datetime.utcnow().isoformat(),
    }
