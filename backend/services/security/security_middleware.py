"""
Security Middleware - Request Protection Layer

Provides comprehensive request-level security for NeuroBotanica:
- Request fingerprinting
- Threat detection
- IP-based blocking
- Request validation
- Security headers

Trade Secret Note:
- Fingerprinting algorithms are protected
- Threat detection rules are confidential
- Blocking thresholds are proprietary

Reference: NeuroBotanica MVP Development Plan - Week 12
"""

from dataclasses import dataclass, field
from datetime import datetime, timedelta
from enum import Enum
from typing import Optional, Dict, List, Set, Any, Tuple
from collections import defaultdict
import time
import hashlib
import threading
import logging
import re

from starlette.middleware.base import BaseHTTPMiddleware
from starlette.requests import Request
from starlette.responses import Response, JSONResponse
from fastapi import status

logger = logging.getLogger(__name__)


class ThreatLevel(Enum):
    """Threat classification levels."""
    
    NONE = "none"
    LOW = "low"
    MEDIUM = "medium"
    HIGH = "high"
    CRITICAL = "critical"
    
    @property
    def should_block(self) -> bool:
        """Whether this threat level should block the request."""
        return self in [ThreatLevel.HIGH, ThreatLevel.CRITICAL]
    
    @property
    def requires_logging(self) -> bool:
        """Whether this threat level requires detailed logging."""
        return self != ThreatLevel.NONE


@dataclass
class RequestFingerprint:
    """Fingerprint of a request for tracking and analysis."""
    
    fingerprint_id: str
    
    # Network info
    source_ip: str
    forwarded_for: Optional[str] = None
    
    # Client info
    user_agent: str = ""
    accept_language: Optional[str] = None
    accept_encoding: Optional[str] = None
    
    # Request info
    method: str = ""
    path: str = ""
    query_string: Optional[str] = None
    content_type: Optional[str] = None
    content_length: Optional[int] = None
    
    # Derived
    ua_hash: str = ""
    headers_hash: str = ""
    
    # Timestamp
    timestamp: datetime = field(default_factory=datetime.utcnow)
    
    def to_dict(self) -> Dict[str, Any]:
        return {
            "fingerprint_id": self.fingerprint_id,
            "source_ip": self.source_ip,
            "user_agent": self.user_agent[:100] if self.user_agent else None,
            "method": self.method,
            "path": self.path,
            "ua_hash": self.ua_hash,
            "timestamp": self.timestamp.isoformat(),
        }


@dataclass
class SecurityConfig:
    """Configuration for security middleware."""
    
    # Enable/disable features
    enable_fingerprinting: bool = True
    enable_rate_limiting: bool = True
    enable_threat_detection: bool = True
    enable_ip_blocking: bool = True
    enable_security_headers: bool = True
    
    # IP blocking
    blocked_ips: Set[str] = field(default_factory=set)
    blocked_ranges: List[str] = field(default_factory=list)
    auto_block_threshold: int = 100  # Requests before auto-block
    auto_block_window_seconds: int = 60
    auto_block_duration_seconds: int = 3600  # 1 hour
    
    # Path protection
    protected_paths: List[str] = field(default_factory=lambda: [
        "/api/v1/admin",
        "/api/v1/internal",
    ])
    
    # Request limits
    max_content_length: int = 10 * 1024 * 1024  # 10MB
    max_header_size: int = 8192
    max_query_string_length: int = 2048
    
    # Security headers
    strict_transport_security: bool = True
    content_security_policy: bool = True
    x_frame_options: str = "DENY"
    x_content_type_options: str = "nosniff"
    
    # Suspicious patterns (regex)
    suspicious_patterns: List[str] = field(default_factory=lambda: [
        r"(?i)(union\s+select|insert\s+into|drop\s+table)",  # SQL injection
        r"(?i)(<script|javascript:|on\w+\s*=)",  # XSS
        r"\.\./\.\./",  # Path traversal
        r"(?i)(exec\s*\(|system\s*\(|eval\s*\()",  # Command injection
    ])


# ============================================================================
# TRADE SECRET: Threat Detection Rules
# ============================================================================

_THREAT_PATTERNS = {
    "sql_injection": {
        "patterns": [
            r"(?i)union\s+select",
            r"(?i)insert\s+into",
            r"(?i)drop\s+(table|database)",
            r"(?i)select\s+.*\s+from",
            r"'\s*or\s*'",
            r";\s*--",
        ],
        "threat_level": ThreatLevel.HIGH,
    },
    "xss": {
        "patterns": [
            r"(?i)<script",
            r"(?i)javascript:",
            r"(?i)on\w+\s*=",
            r"(?i)data:text/html",
        ],
        "threat_level": ThreatLevel.HIGH,
    },
    "path_traversal": {
        "patterns": [
            r"\.\./",
            r"\.\.\\",
            r"%2e%2e/",
            r"%2e%2e\\",
        ],
        "threat_level": ThreatLevel.HIGH,
    },
    "command_injection": {
        "patterns": [
            r"(?i)\|\s*\w+",
            r"(?i);\s*\w+",
            r"(?i)\$\(",
            r"(?i)`.*`",
        ],
        "threat_level": ThreatLevel.CRITICAL,
    },
    "suspicious_user_agent": {
        "patterns": [
            r"(?i)(sqlmap|nikto|nmap|masscan|burp|zap|dirbuster)",
            r"(?i)(python-requests|curl|wget).*scan",
        ],
        "threat_level": ThreatLevel.MEDIUM,
    },
}

# Threshold multipliers for different behaviors
_BEHAVIOR_SCORES = {
    "rapid_requests": 2.0,
    "many_404s": 1.5,
    "auth_failures": 3.0,
    "varied_endpoints": 1.2,
    "suspicious_timing": 1.5,
}


class SecurityMiddleware(BaseHTTPMiddleware):
    """Comprehensive security middleware for request protection."""
    
    def __init__(self, app, config: Optional[SecurityConfig] = None):
        super().__init__(app)
        self.config = config or SecurityConfig()
        
        # Tracking state
        self._request_counts: Dict[str, List[float]] = defaultdict(list)
        self._blocked_until: Dict[str, float] = {}
        self._fingerprints: Dict[str, List[RequestFingerprint]] = defaultdict(list)
        
        # Compiled patterns
        self._compiled_patterns = self._compile_patterns()
        
        # Thread safety
        self._lock = threading.Lock()
        
        logger.info("SecurityMiddleware initialized")
    
    async def dispatch(self, request: Request, call_next) -> Response:
        """Process each request through security checks."""
        start_time = time.time()
        
        # Get client IP
        client_ip = self._get_client_ip(request)
        
        # Check IP blocking
        if self.config.enable_ip_blocking:
            if self._is_ip_blocked(client_ip):
                return self._blocked_response(client_ip, "IP address blocked")
        
        # Create fingerprint
        fingerprint = None
        if self.config.enable_fingerprinting:
            fingerprint = self._create_fingerprint(request, client_ip)
        
        # Threat detection
        threat_level = ThreatLevel.NONE
        threat_details = []
        
        if self.config.enable_threat_detection:
            threat_level, threat_details = await self._detect_threats(request, fingerprint)
            
            if threat_level.should_block:
                self._maybe_auto_block(client_ip)
                return self._blocked_response(
                    client_ip,
                    f"Request blocked: {', '.join(threat_details)}"
                )
        
        # Track request for rate limiting
        if self.config.enable_rate_limiting:
            self._track_request(client_ip)
            
            if self._should_auto_block(client_ip):
                self._auto_block_ip(client_ip)
                return self._blocked_response(client_ip, "Rate limit exceeded - IP temporarily blocked")
        
        # Process request
        response = await call_next(request)
        
        # Add security headers
        if self.config.enable_security_headers:
            self._add_security_headers(response)
        
        # Add processing time header
        process_time = (time.time() - start_time) * 1000
        response.headers["X-Process-Time-Ms"] = f"{process_time:.2f}"
        
        # Log if threat detected
        if threat_level.requires_logging:
            logger.warning(
                "Threat detected: level=%s ip=%s path=%s details=%s",
                threat_level.value,
                client_ip,
                request.url.path,
                threat_details
            )
        
        return response
    
    def _get_client_ip(self, request: Request) -> str:
        """Extract client IP from request."""
        # Check X-Forwarded-For header
        forwarded_for = request.headers.get("X-Forwarded-For")
        if forwarded_for:
            # Take first IP (original client)
            return forwarded_for.split(",")[0].strip()
        
        # Check X-Real-IP header
        real_ip = request.headers.get("X-Real-IP")
        if real_ip:
            return real_ip.strip()
        
        # Fall back to direct client
        if request.client:
            return request.client.host
        
        return "unknown"
    
    def _is_ip_blocked(self, ip: str) -> bool:
        """Check if IP is blocked."""
        # Check explicit block list
        if ip in self.config.blocked_ips:
            return True
        
        # Check auto-blocked IPs
        blocked_until = self._blocked_until.get(ip, 0)
        if blocked_until > time.time():
            return True
        
        # Check blocked ranges (CIDR)
        for ip_range in self.config.blocked_ranges:
            if self._ip_in_range(ip, ip_range):
                return True
        
        return False
    
    def _ip_in_range(self, ip: str, cidr: str) -> bool:
        """Check if IP is in CIDR range."""
        try:
            import ipaddress
            network = ipaddress.ip_network(cidr, strict=False)
            return ipaddress.ip_address(ip) in network
        except (ValueError, TypeError):
            return False
    
    def _create_fingerprint(self, request: Request, client_ip: str) -> RequestFingerprint:
        """Create request fingerprint."""
        user_agent = request.headers.get("User-Agent", "")
        
        # Hash user agent for grouping
        ua_hash = hashlib.md5(user_agent.encode()).hexdigest()[:8]
        
        # Hash relevant headers for fingerprinting
        header_data = f"{request.headers.get('Accept', '')}{request.headers.get('Accept-Language', '')}"
        headers_hash = hashlib.md5(header_data.encode()).hexdigest()[:8]
        
        fingerprint = RequestFingerprint(
            fingerprint_id=f"fp_{ua_hash}_{headers_hash}",
            source_ip=client_ip,
            forwarded_for=request.headers.get("X-Forwarded-For"),
            user_agent=user_agent,
            accept_language=request.headers.get("Accept-Language"),
            accept_encoding=request.headers.get("Accept-Encoding"),
            method=request.method,
            path=request.url.path,
            query_string=str(request.query_params),
            content_type=request.headers.get("Content-Type"),
            content_length=int(request.headers.get("Content-Length", 0)),
            ua_hash=ua_hash,
            headers_hash=headers_hash,
        )
        
        # Store fingerprint
        with self._lock:
            self._fingerprints[client_ip].append(fingerprint)
            # Keep only recent fingerprints
            if len(self._fingerprints[client_ip]) > 100:
                self._fingerprints[client_ip] = self._fingerprints[client_ip][-100:]
        
        return fingerprint
    
    async def _detect_threats(
        self,
        request: Request,
        fingerprint: Optional[RequestFingerprint]
    ) -> Tuple[ThreatLevel, List[str]]:
        """Detect threats in request."""
        threats = []
        max_level = ThreatLevel.NONE
        
        # Check URL path
        path_threats = self._check_patterns(request.url.path)
        threats.extend(path_threats)
        
        # Check query string
        if request.query_params:
            query_str = str(request.query_params)
            query_threats = self._check_patterns(query_str)
            threats.extend(query_threats)
        
        # Check headers
        for header_name, header_value in request.headers.items():
            if header_name.lower() in ["cookie", "authorization"]:
                continue  # Skip sensitive headers
            header_threats = self._check_patterns(header_value)
            threats.extend(header_threats)
        
        # Check user agent
        if fingerprint:
            ua_threats = self._check_user_agent(fingerprint.user_agent)
            threats.extend(ua_threats)
        
        # Check body for POST/PUT/PATCH
        if request.method in ["POST", "PUT", "PATCH"]:
            try:
                body = await request.body()
                if body:
                    body_str = body.decode("utf-8", errors="ignore")
                    body_threats = self._check_patterns(body_str[:10000])  # Limit check size
                    threats.extend(body_threats)
            except Exception:
                pass  # Ignore body read errors
        
        # Check request size limits
        content_length = int(request.headers.get("Content-Length", 0))
        if content_length > self.config.max_content_length:
            threats.append(("size_limit", ThreatLevel.MEDIUM))
        
        # Determine max threat level
        for threat_name, level in threats:
            if level.value > max_level.value:
                max_level = level
        
        threat_names = [t[0] for t in threats]
        return max_level, threat_names
    
    def _check_patterns(self, text: str) -> List[Tuple[str, ThreatLevel]]:
        """Check text against threat patterns."""
        threats = []
        
        for threat_name, threat_info in _THREAT_PATTERNS.items():
            for pattern in self._compiled_patterns.get(threat_name, []):
                if pattern.search(text):
                    threats.append((threat_name, threat_info["threat_level"]))
                    break  # One match per category is enough
        
        return threats
    
    def _check_user_agent(self, user_agent: str) -> List[Tuple[str, ThreatLevel]]:
        """Check user agent for suspicious patterns."""
        threats = []
        
        if not user_agent:
            threats.append(("empty_user_agent", ThreatLevel.LOW))
        
        # Check against suspicious UA patterns
        for pattern in self._compiled_patterns.get("suspicious_user_agent", []):
            if pattern.search(user_agent):
                threats.append(("suspicious_user_agent", ThreatLevel.MEDIUM))
                break
        
        return threats
    
    def _compile_patterns(self) -> Dict[str, List[re.Pattern]]:
        """Compile threat detection patterns."""
        compiled = {}
        for threat_name, threat_info in _THREAT_PATTERNS.items():
            compiled[threat_name] = [
                re.compile(p) for p in threat_info["patterns"]
            ]
        return compiled
    
    def _track_request(self, client_ip: str):
        """Track request for rate limiting."""
        now = time.time()
        
        with self._lock:
            self._request_counts[client_ip].append(now)
            
            # Clean old entries
            cutoff = now - self.config.auto_block_window_seconds
            self._request_counts[client_ip] = [
                t for t in self._request_counts[client_ip] if t > cutoff
            ]
    
    def _should_auto_block(self, client_ip: str) -> bool:
        """Check if IP should be auto-blocked."""
        count = len(self._request_counts.get(client_ip, []))
        return count >= self.config.auto_block_threshold
    
    def _auto_block_ip(self, client_ip: str):
        """Auto-block an IP address."""
        with self._lock:
            self._blocked_until[client_ip] = (
                time.time() + self.config.auto_block_duration_seconds
            )
        
        logger.warning(
            "IP auto-blocked: %s for %d seconds",
            client_ip,
            self.config.auto_block_duration_seconds
        )
    
    def _maybe_auto_block(self, client_ip: str):
        """Maybe block IP after threat detection."""
        # Count recent threats
        with self._lock:
            self._request_counts[client_ip].append(time.time())
    
    def _add_security_headers(self, response: Response):
        """Add security headers to response."""
        if self.config.strict_transport_security:
            response.headers["Strict-Transport-Security"] = (
                "max-age=31536000; includeSubDomains"
            )
        
        if self.config.content_security_policy:
            response.headers["Content-Security-Policy"] = (
                "default-src 'self'; script-src 'self'; style-src 'self' 'unsafe-inline'"
            )
        
        response.headers["X-Frame-Options"] = self.config.x_frame_options
        response.headers["X-Content-Type-Options"] = self.config.x_content_type_options
        response.headers["X-XSS-Protection"] = "1; mode=block"
        response.headers["Referrer-Policy"] = "strict-origin-when-cross-origin"
    
    def _blocked_response(self, client_ip: str, reason: str) -> JSONResponse:
        """Return blocked response."""
        return JSONResponse(
            status_code=status.HTTP_403_FORBIDDEN,
            content={
                "error": "Forbidden",
                "message": reason,
                "blocked_ip": client_ip[:8] + "...",  # Partial IP for privacy
            }
        )
    
    # Public methods for management
    
    def block_ip(self, ip: str, duration_seconds: Optional[int] = None):
        """Manually block an IP address."""
        if duration_seconds:
            with self._lock:
                self._blocked_until[ip] = time.time() + duration_seconds
        else:
            self.config.blocked_ips.add(ip)
        
        logger.info("IP blocked: %s", ip)
    
    def unblock_ip(self, ip: str):
        """Unblock an IP address."""
        self.config.blocked_ips.discard(ip)
        with self._lock:
            self._blocked_until.pop(ip, None)
        
        logger.info("IP unblocked: %s", ip)
    
    def get_blocked_ips(self) -> List[Dict[str, Any]]:
        """Get list of blocked IPs."""
        now = time.time()
        blocked = []
        
        # Permanent blocks
        for ip in self.config.blocked_ips:
            blocked.append({
                "ip": ip,
                "type": "permanent",
                "expires": None,
            })
        
        # Temporary blocks
        for ip, until in self._blocked_until.items():
            if until > now:
                blocked.append({
                    "ip": ip,
                    "type": "temporary",
                    "expires": datetime.fromtimestamp(until).isoformat(),
                    "remaining_seconds": int(until - now),
                })
        
        return blocked
    
    def get_request_stats(self, ip: str) -> Dict[str, Any]:
        """Get request statistics for an IP."""
        now = time.time()
        counts = self._request_counts.get(ip, [])
        fingerprints = self._fingerprints.get(ip, [])
        
        # Clean and count
        recent_counts = [t for t in counts if t > now - self.config.auto_block_window_seconds]
        
        return {
            "ip": ip,
            "requests_in_window": len(recent_counts),
            "window_seconds": self.config.auto_block_window_seconds,
            "threshold": self.config.auto_block_threshold,
            "is_blocked": self._is_ip_blocked(ip),
            "fingerprint_count": len(fingerprints),
            "unique_user_agents": len(set(fp.ua_hash for fp in fingerprints)),
        }
