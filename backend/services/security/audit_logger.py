"""
Audit Logger - Compliance-Grade Event Logging

Provides comprehensive audit logging for NeuroBotanica:
- HIPAA-aligned event logging
- SOC2 compliance support
- Immutable audit trails
- Security event detection and alerting
- Request/response correlation
- PII handling and redaction

Trade Secret Note:
- Event classification rules are protected
- Alert thresholds are confidential
- Anomaly detection algorithms are proprietary

Reference: NeuroBotanica MVP Development Plan - Week 12
"""

from dataclasses import dataclass, field
from datetime import datetime, timedelta
from enum import Enum
from typing import Optional, Dict, List, Any, Callable
from collections import defaultdict
import threading
import hashlib
import json
import uuid
import logging
import time

logger = logging.getLogger(__name__)


class AuditEventType(Enum):
    """Types of audit events."""
    
    # Authentication events
    AUTH_LOGIN = "auth.login"
    AUTH_LOGOUT = "auth.logout"
    AUTH_LOGIN_FAILED = "auth.login_failed"
    AUTH_TOKEN_ISSUED = "auth.token_issued"
    AUTH_TOKEN_REVOKED = "auth.token_revoked"
    AUTH_PASSWORD_CHANGED = "auth.password_changed"
    AUTH_MFA_ENABLED = "auth.mfa_enabled"
    AUTH_MFA_DISABLED = "auth.mfa_disabled"
    
    # API key events
    KEY_CREATED = "key.created"
    KEY_REVOKED = "key.revoked"
    KEY_ROTATED = "key.rotated"
    KEY_USED = "key.used"
    KEY_RATE_LIMITED = "key.rate_limited"
    
    # Data access events
    DATA_READ = "data.read"
    DATA_CREATE = "data.create"
    DATA_UPDATE = "data.update"
    DATA_DELETE = "data.delete"
    DATA_EXPORT = "data.export"
    DATA_BULK_ACCESS = "data.bulk_access"
    
    # Service access events
    SERVICE_CHEMPATH = "service.chempath"
    SERVICE_TOXPATH = "service.toxpath"
    SERVICE_REGPATH = "service.regpath"
    SERVICE_PATENTPATH = "service.patentpath"
    SERVICE_OMNIPATH = "service.omnipath"
    
    # Security events
    SECURITY_THREAT_DETECTED = "security.threat_detected"
    SECURITY_RATE_LIMIT_EXCEEDED = "security.rate_limit_exceeded"
    SECURITY_INVALID_TOKEN = "security.invalid_token"
    SECURITY_PERMISSION_DENIED = "security.permission_denied"
    SECURITY_SUSPICIOUS_ACTIVITY = "security.suspicious_activity"
    SECURITY_IP_BLOCKED = "security.ip_blocked"
    
    # Admin events
    ADMIN_USER_CREATED = "admin.user_created"
    ADMIN_USER_DELETED = "admin.user_deleted"
    ADMIN_PERMISSION_CHANGED = "admin.permission_changed"
    ADMIN_CONFIG_CHANGED = "admin.config_changed"
    
    # System events
    SYSTEM_STARTUP = "system.startup"
    SYSTEM_SHUTDOWN = "system.shutdown"
    SYSTEM_ERROR = "system.error"
    SYSTEM_HEALTH_CHECK = "system.health_check"
    
    @property
    def severity(self) -> str:
        """Get severity level for event type."""
        return _EVENT_SEVERITY.get(self.value, "info")
    
    @property
    def is_security_event(self) -> bool:
        """Check if this is a security event."""
        return self.value.startswith("security.")
    
    @property
    def requires_alert(self) -> bool:
        """Check if event type requires immediate alerting."""
        return self in _ALERT_EVENTS


class ComplianceLevel(Enum):
    """Compliance framework levels."""
    
    BASIC = "basic"
    HIPAA = "hipaa"
    SOC2 = "soc2"
    HIPAA_SOC2 = "hipaa_soc2"
    
    @property
    def retention_days(self) -> int:
        """Minimum retention period in days."""
        return {
            "basic": 90,
            "hipaa": 2190,  # 6 years
            "soc2": 365,
            "hipaa_soc2": 2190,
        }.get(self.value, 90)
    
    @property
    def requires_pii_redaction(self) -> bool:
        """Whether PII should be redacted in logs."""
        return self in [ComplianceLevel.HIPAA, ComplianceLevel.HIPAA_SOC2]


# ============================================================================
# TRADE SECRET: Event Severity Classification
# ============================================================================

_EVENT_SEVERITY = {
    "auth.login": "info",
    "auth.logout": "info",
    "auth.login_failed": "warning",
    "auth.token_issued": "info",
    "auth.token_revoked": "warning",
    "auth.password_changed": "info",
    "auth.mfa_enabled": "info",
    "auth.mfa_disabled": "warning",
    "key.created": "info",
    "key.revoked": "warning",
    "key.rotated": "info",
    "key.used": "debug",
    "key.rate_limited": "warning",
    "data.read": "debug",
    "data.create": "info",
    "data.update": "info",
    "data.delete": "warning",
    "data.export": "info",
    "data.bulk_access": "info",
    "service.chempath": "info",
    "service.toxpath": "info",
    "service.regpath": "info",
    "service.patentpath": "info",
    "service.omnipath": "info",
    "security.threat_detected": "critical",
    "security.rate_limit_exceeded": "warning",
    "security.invalid_token": "warning",
    "security.permission_denied": "warning",
    "security.suspicious_activity": "warning",
    "security.ip_blocked": "critical",
    "admin.user_created": "info",
    "admin.user_deleted": "warning",
    "admin.permission_changed": "warning",
    "admin.config_changed": "warning",
    "system.startup": "info",
    "system.shutdown": "info",
    "system.error": "error",
    "system.health_check": "debug",
}

_ALERT_EVENTS = {
    AuditEventType.SECURITY_THREAT_DETECTED,
    AuditEventType.SECURITY_IP_BLOCKED,
    AuditEventType.AUTH_LOGIN_FAILED,  # After threshold
    AuditEventType.ADMIN_PERMISSION_CHANGED,
    AuditEventType.SYSTEM_ERROR,
}


@dataclass
class SecurityEvent:
    """Security-specific event with threat context."""
    
    event_id: str
    event_type: AuditEventType
    threat_level: str  # low, medium, high, critical
    
    source_ip: str
    user_id: Optional[str] = None
    api_key_id: Optional[str] = None
    
    description: str = ""
    indicators: List[str] = field(default_factory=list)
    recommended_action: Optional[str] = None
    
    timestamp: datetime = field(default_factory=datetime.utcnow)
    is_blocked: bool = False
    
    def to_dict(self) -> Dict[str, Any]:
        return {
            "event_id": self.event_id,
            "event_type": self.event_type.value,
            "threat_level": self.threat_level,
            "source_ip": self.source_ip,
            "user_id": self.user_id,
            "api_key_id": self.api_key_id,
            "description": self.description,
            "indicators": self.indicators,
            "recommended_action": self.recommended_action,
            "timestamp": self.timestamp.isoformat(),
            "is_blocked": self.is_blocked,
        }


@dataclass
class AuditEvent:
    """Immutable audit event record."""
    
    event_id: str
    event_type: AuditEventType
    timestamp: datetime
    
    # Actor information
    user_id: Optional[str] = None
    api_key_id: Optional[str] = None
    session_id: Optional[str] = None
    
    # Request context
    request_id: Optional[str] = None
    source_ip: Optional[str] = None
    user_agent: Optional[str] = None
    endpoint: Optional[str] = None
    method: Optional[str] = None
    
    # Event details
    resource_type: Optional[str] = None
    resource_id: Optional[str] = None
    action: Optional[str] = None
    details: Dict[str, Any] = field(default_factory=dict)
    
    # Response
    status_code: Optional[int] = None
    response_time_ms: Optional[float] = None
    success: bool = True
    error_message: Optional[str] = None
    
    # Metadata
    severity: str = "info"
    tags: List[str] = field(default_factory=list)
    
    # Integrity
    checksum: Optional[str] = None
    
    def __post_init__(self):
        """Calculate checksum for integrity verification."""
        if self.checksum is None:
            self.checksum = self._calculate_checksum()
    
    def _calculate_checksum(self) -> str:
        """Calculate SHA-256 checksum of event data."""
        data = {
            "event_id": self.event_id,
            "event_type": self.event_type.value,
            "timestamp": self.timestamp.isoformat(),
            "user_id": self.user_id,
            "request_id": self.request_id,
            "resource_type": self.resource_type,
            "resource_id": self.resource_id,
            "action": self.action,
            "success": self.success,
        }
        data_str = json.dumps(data, sort_keys=True)
        return hashlib.sha256(data_str.encode()).hexdigest()[:16]
    
    def verify_integrity(self) -> bool:
        """Verify event integrity via checksum."""
        return self.checksum == self._calculate_checksum()
    
    def to_dict(self, redact_pii: bool = False) -> Dict[str, Any]:
        """Convert to dictionary for storage/API response."""
        data = {
            "event_id": self.event_id,
            "event_type": self.event_type.value,
            "timestamp": self.timestamp.isoformat(),
            "user_id": self._redact(self.user_id) if redact_pii else self.user_id,
            "api_key_id": self._redact(self.api_key_id) if redact_pii else self.api_key_id,
            "session_id": self.session_id,
            "request_id": self.request_id,
            "source_ip": self._redact_ip(self.source_ip) if redact_pii else self.source_ip,
            "endpoint": self.endpoint,
            "method": self.method,
            "resource_type": self.resource_type,
            "resource_id": self.resource_id,
            "action": self.action,
            "details": self._redact_details(self.details) if redact_pii else self.details,
            "status_code": self.status_code,
            "response_time_ms": self.response_time_ms,
            "success": self.success,
            "error_message": self.error_message,
            "severity": self.severity,
            "tags": self.tags,
            "checksum": self.checksum,
        }
        return data
    
    def _redact(self, value: Optional[str]) -> Optional[str]:
        """Redact sensitive value."""
        if value is None:
            return None
        if len(value) <= 4:
            return "****"
        return value[:2] + "****" + value[-2:]
    
    def _redact_ip(self, ip: Optional[str]) -> Optional[str]:
        """Redact IP address."""
        if ip is None:
            return None
        parts = ip.split(".")
        if len(parts) == 4:
            return f"{parts[0]}.{parts[1]}.*.*"
        return "****"
    
    def _redact_details(self, details: Dict) -> Dict:
        """Redact sensitive fields from details."""
        sensitive_fields = {"email", "password", "token", "key", "secret", "ssn", "dob"}
        redacted = {}
        for key, value in details.items():
            if any(s in key.lower() for s in sensitive_fields):
                redacted[key] = "[REDACTED]"
            elif isinstance(value, dict):
                redacted[key] = self._redact_details(value)
            else:
                redacted[key] = value
        return redacted


@dataclass
class AuditConfig:
    """Configuration for audit logger."""
    
    compliance_level: ComplianceLevel = ComplianceLevel.BASIC
    
    # Storage
    max_events_in_memory: int = 10000
    flush_interval_seconds: int = 60
    
    # Alerting
    enable_alerts: bool = True
    alert_callback: Optional[Callable[[SecurityEvent], None]] = None
    
    # Failed login tracking
    failed_login_threshold: int = 5
    failed_login_window_seconds: int = 300
    
    # Rate limit violation tracking
    rate_limit_alert_threshold: int = 10
    rate_limit_alert_window_seconds: int = 60
    
    # PII handling
    redact_pii_in_logs: bool = False
    
    # Event filtering
    log_debug_events: bool = False
    log_health_checks: bool = False


class AuditLogger:
    """Main audit logging service with compliance support."""
    
    def __init__(self, config: Optional[AuditConfig] = None):
        self.config = config or AuditConfig()
        
        # Event storage
        self._events: List[AuditEvent] = []
        self._security_events: List[SecurityEvent] = []
        
        # Tracking for anomaly detection
        self._failed_logins: Dict[str, List[float]] = defaultdict(list)
        self._rate_limit_violations: Dict[str, List[float]] = defaultdict(list)
        
        # Thread safety
        self._lock = threading.Lock()
        self._last_flush = time.time()
        
        logger.info(
            "AuditLogger initialized with compliance level: %s",
            self.config.compliance_level.value
        )
    
    def log(
        self,
        event_type: AuditEventType,
        user_id: Optional[str] = None,
        api_key_id: Optional[str] = None,
        request_id: Optional[str] = None,
        source_ip: Optional[str] = None,
        endpoint: Optional[str] = None,
        method: Optional[str] = None,
        resource_type: Optional[str] = None,
        resource_id: Optional[str] = None,
        action: Optional[str] = None,
        details: Optional[Dict] = None,
        status_code: Optional[int] = None,
        response_time_ms: Optional[float] = None,
        success: bool = True,
        error_message: Optional[str] = None,
        tags: Optional[List[str]] = None,
    ) -> AuditEvent:
        """Log an audit event.
        
        Args:
            event_type: Type of event
            user_id: User ID if authenticated
            api_key_id: API key ID if used
            request_id: Correlation ID for request
            source_ip: Client IP address
            endpoint: API endpoint path
            method: HTTP method
            resource_type: Type of resource accessed
            resource_id: ID of resource accessed
            action: Specific action performed
            details: Additional event details
            status_code: HTTP response status
            response_time_ms: Request duration
            success: Whether operation succeeded
            error_message: Error message if failed
            tags: Optional tags for filtering
            
        Returns:
            Created AuditEvent
        """
        # Filter debug/health events if configured
        if not self.config.log_debug_events and event_type.severity == "debug":
            return None
        if not self.config.log_health_checks and event_type == AuditEventType.SYSTEM_HEALTH_CHECK:
            return None
        
        event = AuditEvent(
            event_id=f"evt_{uuid.uuid4().hex[:16]}",
            event_type=event_type,
            timestamp=datetime.utcnow(),
            user_id=user_id,
            api_key_id=api_key_id,
            request_id=request_id,
            source_ip=source_ip,
            endpoint=endpoint,
            method=method,
            resource_type=resource_type,
            resource_id=resource_id,
            action=action,
            details=details or {},
            status_code=status_code,
            response_time_ms=response_time_ms,
            success=success,
            error_message=error_message,
            severity=event_type.severity,
            tags=tags or [],
        )
        
        with self._lock:
            self._events.append(event)
            
            # Check for memory limit
            if len(self._events) > self.config.max_events_in_memory:
                self._flush_old_events()
        
        # Check for security conditions
        self._check_security_conditions(event)
        
        # Log to standard logger
        log_level = {
            "debug": logging.DEBUG,
            "info": logging.INFO,
            "warning": logging.WARNING,
            "error": logging.ERROR,
            "critical": logging.CRITICAL,
        }.get(event.severity, logging.INFO)
        
        logger.log(
            log_level,
            "Audit: %s | user=%s | ip=%s | endpoint=%s | success=%s",
            event_type.value,
            user_id or "anonymous",
            source_ip or "unknown",
            endpoint or "N/A",
            success
        )
        
        return event
    
    def log_security_event(
        self,
        event_type: AuditEventType,
        threat_level: str,
        source_ip: str,
        description: str,
        user_id: Optional[str] = None,
        api_key_id: Optional[str] = None,
        indicators: Optional[List[str]] = None,
        recommended_action: Optional[str] = None,
        is_blocked: bool = False,
    ) -> SecurityEvent:
        """Log a security-specific event.
        
        Args:
            event_type: Security event type
            threat_level: Threat severity (low/medium/high/critical)
            source_ip: Source IP address
            description: Human-readable description
            user_id: User ID if known
            api_key_id: API key if used
            indicators: List of threat indicators
            recommended_action: Suggested response
            is_blocked: Whether request was blocked
            
        Returns:
            Created SecurityEvent
        """
        event = SecurityEvent(
            event_id=f"sec_{uuid.uuid4().hex[:16]}",
            event_type=event_type,
            threat_level=threat_level,
            source_ip=source_ip,
            user_id=user_id,
            api_key_id=api_key_id,
            description=description,
            indicators=indicators or [],
            recommended_action=recommended_action,
            is_blocked=is_blocked,
        )
        
        with self._lock:
            self._security_events.append(event)
        
        # Trigger alert if configured
        if self.config.enable_alerts and event_type.requires_alert:
            self._trigger_alert(event)
        
        logger.warning(
            "SECURITY: %s | level=%s | ip=%s | blocked=%s | %s",
            event_type.value,
            threat_level,
            source_ip,
            is_blocked,
            description
        )
        
        return event
    
    def get_events(
        self,
        event_type: Optional[AuditEventType] = None,
        user_id: Optional[str] = None,
        start_time: Optional[datetime] = None,
        end_time: Optional[datetime] = None,
        limit: int = 100,
        include_security: bool = False,
    ) -> List[Dict]:
        """Query audit events.
        
        Args:
            event_type: Filter by event type
            user_id: Filter by user ID
            start_time: Filter by start time
            end_time: Filter by end time
            limit: Maximum events to return
            include_security: Include security events
            
        Returns:
            List of event dictionaries
        """
        events = self._events.copy()
        
        # Apply filters
        if event_type:
            events = [e for e in events if e.event_type == event_type]
        if user_id:
            events = [e for e in events if e.user_id == user_id]
        if start_time:
            events = [e for e in events if e.timestamp >= start_time]
        if end_time:
            events = [e for e in events if e.timestamp <= end_time]
        
        # Sort by timestamp descending
        events = sorted(events, key=lambda e: e.timestamp, reverse=True)
        
        # Apply limit
        events = events[:limit]
        
        # Convert to dict
        redact = self.config.compliance_level.requires_pii_redaction
        result = [e.to_dict(redact_pii=redact) for e in events]
        
        # Include security events if requested
        if include_security:
            sec_events = [e.to_dict() for e in self._security_events[:limit]]
            result.extend(sec_events)
        
        return result
    
    def get_security_summary(
        self,
        hours: int = 24
    ) -> Dict[str, Any]:
        """Get security event summary.
        
        Args:
            hours: Look back period in hours
            
        Returns:
            Summary statistics
        """
        cutoff = datetime.utcnow() - timedelta(hours=hours)
        
        recent_security = [
            e for e in self._security_events
            if e.timestamp >= cutoff
        ]
        
        # Count by threat level
        by_level = defaultdict(int)
        for e in recent_security:
            by_level[e.threat_level] += 1
        
        # Count by type
        by_type = defaultdict(int)
        for e in recent_security:
            by_type[e.event_type.value] += 1
        
        # Unique IPs
        unique_ips = set(e.source_ip for e in recent_security)
        
        return {
            "period_hours": hours,
            "total_events": len(recent_security),
            "by_threat_level": dict(by_level),
            "by_event_type": dict(by_type),
            "unique_source_ips": len(unique_ips),
            "blocked_requests": sum(1 for e in recent_security if e.is_blocked),
        }
    
    def verify_integrity(self) -> Dict[str, Any]:
        """Verify integrity of audit log.
        
        Returns:
            Integrity check results
        """
        total = len(self._events)
        valid = sum(1 for e in self._events if e.verify_integrity())
        
        return {
            "total_events": total,
            "valid_events": valid,
            "invalid_events": total - valid,
            "integrity_percentage": (valid / total * 100) if total > 0 else 100,
            "compliance_level": self.config.compliance_level.value,
        }
    
    def _check_security_conditions(self, event: AuditEvent):
        """Check for security conditions that require alerting."""
        now = time.time()
        
        # Failed login tracking
        if event.event_type == AuditEventType.AUTH_LOGIN_FAILED:
            identifier = event.source_ip or "unknown"
            self._failed_logins[identifier].append(now)
            
            # Clean old entries
            cutoff = now - self.config.failed_login_window_seconds
            self._failed_logins[identifier] = [
                t for t in self._failed_logins[identifier] if t > cutoff
            ]
            
            # Check threshold
            if len(self._failed_logins[identifier]) >= self.config.failed_login_threshold:
                self.log_security_event(
                    event_type=AuditEventType.SECURITY_SUSPICIOUS_ACTIVITY,
                    threat_level="high",
                    source_ip=event.source_ip or "unknown",
                    description=f"Multiple failed login attempts: {len(self._failed_logins[identifier])}",
                    user_id=event.user_id,
                    indicators=["brute_force", "credential_stuffing"],
                    recommended_action="Consider IP blocking or account lockout",
                )
        
        # Rate limit violation tracking
        if event.event_type == AuditEventType.KEY_RATE_LIMITED:
            identifier = event.api_key_id or event.source_ip or "unknown"
            self._rate_limit_violations[identifier].append(now)
            
            # Clean old entries
            cutoff = now - self.config.rate_limit_alert_window_seconds
            self._rate_limit_violations[identifier] = [
                t for t in self._rate_limit_violations[identifier] if t > cutoff
            ]
            
            # Check threshold
            if len(self._rate_limit_violations[identifier]) >= self.config.rate_limit_alert_threshold:
                self.log_security_event(
                    event_type=AuditEventType.SECURITY_RATE_LIMIT_EXCEEDED,
                    threat_level="medium",
                    source_ip=event.source_ip or "unknown",
                    description=f"Excessive rate limit violations: {len(self._rate_limit_violations[identifier])}",
                    api_key_id=event.api_key_id,
                    indicators=["abuse", "scraping", "ddos"],
                    recommended_action="Review API key usage patterns",
                )
    
    def _trigger_alert(self, event: SecurityEvent):
        """Trigger alert for security event."""
        if self.config.alert_callback:
            try:
                self.config.alert_callback(event)
            except Exception as e:
                logger.error("Alert callback failed: %s", e)
        
        # Always log critical alerts
        if event.threat_level == "critical":
            logger.critical(
                "CRITICAL SECURITY ALERT: %s - %s",
                event.event_type.value,
                event.description
            )
    
    def _flush_old_events(self):
        """Flush old events to make room for new ones."""
        # Keep most recent events
        keep_count = int(self.config.max_events_in_memory * 0.8)
        self._events = self._events[-keep_count:]
        
        logger.debug("Flushed old audit events, keeping %d", keep_count)


# Singleton instance
_audit_logger: Optional[AuditLogger] = None


def get_audit_logger(config: Optional[AuditConfig] = None) -> AuditLogger:
    """Get or create the global audit logger instance."""
    global _audit_logger
    
    if _audit_logger is None:
        _audit_logger = AuditLogger(config)
    
    return _audit_logger
