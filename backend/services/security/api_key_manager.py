"""
API Key Manager - Tiered Access Control

Provides comprehensive API key management for NeuroBotanica:
- Tiered access with different quotas and permissions
- Scoped access to specific endpoints/features
- Key generation with cryptographic security
- Key rotation and revocation
- Usage tracking and analytics

Trade Secret Note:
- Key generation algorithms are protected
- Tier pricing thresholds are confidential
- Scope mappings are proprietary

Reference: NeuroBotanica MVP Development Plan - Week 12
"""

from dataclasses import dataclass, field
from datetime import datetime, timedelta
from enum import Enum
from typing import Optional, Dict, List, Set, Any
import secrets
import hashlib
import base64
import threading
import logging
import uuid
import json

logger = logging.getLogger(__name__)


class APIKeyTier(Enum):
    """API key tiers with feature access levels."""
    
    FREE = "free"
    BASIC = "basic"
    PROFESSIONAL = "professional"
    ENTERPRISE = "enterprise"
    INTERNAL = "internal"
    
    @property
    def display_name(self) -> str:
        return {
            "free": "Free Tier",
            "basic": "Basic ($500/mo)",
            "professional": "Professional ($2,500/mo)",
            "enterprise": "Enterprise ($7,500+/mo)",
            "internal": "Internal Service",
        }.get(self.value, self.value.title())
    
    @property
    def monthly_price_usd(self) -> Optional[int]:
        return {
            "free": 0,
            "basic": 500,
            "professional": 2500,
            "enterprise": 7500,
            "internal": None,
        }.get(self.value)
    
    @property
    def default_scopes(self) -> List[str]:
        """Default scopes for this tier - TRADE SECRET."""
        return _TIER_SCOPES.get(self.value, [])
    
    @property
    def max_keys(self) -> int:
        """Maximum API keys per account."""
        return {
            "free": 1,
            "basic": 3,
            "professional": 10,
            "enterprise": 50,
            "internal": 100,
        }.get(self.value, 1)


class APIKeyScope(Enum):
    """API key access scopes."""
    
    # Read operations
    READ_STUDIES = "read:studies"
    READ_COMPOUNDS = "read:compounds"
    READ_FDA = "read:fda"
    READ_CONFORMERS = "read:conformers"
    READ_EVIDENCE = "read:evidence"
    
    # Write operations
    WRITE_STUDIES = "write:studies"
    WRITE_COMPOUNDS = "write:compounds"
    WRITE_ANALYSES = "write:analyses"
    
    # Service-specific scopes
    CHEMPATH_ACCESS = "chempath:access"
    TOXPATH_ACCESS = "toxpath:access"
    REGPATH_ACCESS = "regpath:access"
    PATENTPATH_ACCESS = "patentpath:access"
    OMNIPATH_ACCESS = "omnipath:access"
    
    # Admin scopes
    ADMIN_USERS = "admin:users"
    ADMIN_KEYS = "admin:keys"
    ADMIN_AUDIT = "admin:audit"
    
    @classmethod
    def all_read(cls) -> List[str]:
        return [s.value for s in cls if s.value.startswith("read:")]
    
    @classmethod
    def all_write(cls) -> List[str]:
        return [s.value for s in cls if s.value.startswith("write:")]
    
    @classmethod
    def all_services(cls) -> List[str]:
        return [s.value for s in cls if ":access" in s.value]


# ============================================================================
# TRADE SECRET: Tier Scopes Mapping
# ============================================================================

_TIER_SCOPES = {
    "free": [
        "read:studies",
        "read:compounds",
        "read:evidence",
    ],
    "basic": [
        "read:studies",
        "read:compounds",
        "read:fda",
        "read:conformers",
        "read:evidence",
        "write:analyses",
        "chempath:access",
    ],
    "professional": [
        "read:studies",
        "read:compounds",
        "read:fda",
        "read:conformers",
        "read:evidence",
        "write:studies",
        "write:compounds",
        "write:analyses",
        "chempath:access",
        "toxpath:access",
        "regpath:access",
    ],
    "enterprise": [
        "read:studies",
        "read:compounds",
        "read:fda",
        "read:conformers",
        "read:evidence",
        "write:studies",
        "write:compounds",
        "write:analyses",
        "chempath:access",
        "toxpath:access",
        "regpath:access",
        "patentpath:access",
        "omnipath:access",
    ],
    "internal": [
        # All scopes
        "read:studies",
        "read:compounds",
        "read:fda",
        "read:conformers",
        "read:evidence",
        "write:studies",
        "write:compounds",
        "write:analyses",
        "chempath:access",
        "toxpath:access",
        "regpath:access",
        "patentpath:access",
        "omnipath:access",
        "admin:users",
        "admin:keys",
        "admin:audit",
    ],
}


@dataclass
class APIKey:
    """API key data structure."""
    
    key_id: str
    key_hash: str  # Never store raw key
    tier: APIKeyTier
    scopes: List[str]
    
    # Ownership
    owner_id: str
    owner_email: Optional[str] = None
    organization: Optional[str] = None
    
    # Metadata
    name: str = "Default API Key"
    description: Optional[str] = None
    
    # Lifecycle
    created_at: datetime = field(default_factory=datetime.utcnow)
    expires_at: Optional[datetime] = None
    last_used_at: Optional[datetime] = None
    
    # Status
    is_active: bool = True
    is_revoked: bool = False
    revoked_at: Optional[datetime] = None
    revoked_reason: Optional[str] = None
    
    # Usage
    request_count: int = 0
    last_ip: Optional[str] = None
    
    # Rate limit overrides (optional)
    custom_rpm: Optional[int] = None
    custom_rph: Optional[int] = None
    custom_rpd: Optional[int] = None
    
    def to_dict(self, include_sensitive: bool = False) -> Dict[str, Any]:
        """Convert to dictionary for API response."""
        data = {
            "key_id": self.key_id,
            "tier": self.tier.value,
            "tier_display": self.tier.display_name,
            "scopes": self.scopes,
            "name": self.name,
            "description": self.description,
            "owner_id": self.owner_id,
            "organization": self.organization,
            "created_at": self.created_at.isoformat(),
            "expires_at": self.expires_at.isoformat() if self.expires_at else None,
            "last_used_at": self.last_used_at.isoformat() if self.last_used_at else None,
            "is_active": self.is_active,
            "is_revoked": self.is_revoked,
            "request_count": self.request_count,
        }
        
        if include_sensitive:
            data["key_hash"] = self.key_hash
            data["last_ip"] = self.last_ip
            data["revoked_reason"] = self.revoked_reason
        
        return data
    
    def has_scope(self, scope: str) -> bool:
        """Check if key has a specific scope."""
        return scope in self.scopes
    
    def has_any_scope(self, scopes: List[str]) -> bool:
        """Check if key has any of the specified scopes."""
        return any(s in self.scopes for s in scopes)
    
    def has_all_scopes(self, scopes: List[str]) -> bool:
        """Check if key has all specified scopes."""
        return all(s in self.scopes for s in scopes)
    
    @property
    def is_expired(self) -> bool:
        """Check if key is expired."""
        if self.expires_at is None:
            return False
        return datetime.utcnow() > self.expires_at
    
    @property
    def is_valid(self) -> bool:
        """Check if key is currently valid."""
        return self.is_active and not self.is_revoked and not self.is_expired


@dataclass
class APIKeyRequest:
    """Request to create a new API key."""
    
    owner_id: str
    tier: APIKeyTier = APIKeyTier.FREE
    name: str = "API Key"
    description: Optional[str] = None
    scopes: Optional[List[str]] = None  # Override default tier scopes
    expires_in_days: Optional[int] = None
    organization: Optional[str] = None
    owner_email: Optional[str] = None


@dataclass
class APIKeyResponse:
    """Response after creating a new API key."""
    
    key: APIKey
    raw_key: str  # Only returned once on creation
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for API response."""
        return {
            "api_key": self.raw_key,  # Only time raw key is visible
            "key_id": self.key.key_id,
            "tier": self.key.tier.value,
            "scopes": self.key.scopes,
            "expires_at": self.key.expires_at.isoformat() if self.key.expires_at else None,
            "warning": "Store this API key securely. It will not be shown again.",
        }


class APIKeyManager:
    """Manages API keys with secure generation and validation."""
    
    # Key format: nb_[tier prefix]_[random bytes]
    KEY_PREFIX = "nb"
    KEY_LENGTH = 32  # bytes for random part
    
    def __init__(self):
        # In-memory storage (would be database in production)
        self._keys: Dict[str, APIKey] = {}  # key_hash -> APIKey
        self._keys_by_id: Dict[str, str] = {}  # key_id -> key_hash
        self._keys_by_owner: Dict[str, List[str]] = {}  # owner_id -> [key_hashes]
        self._lock = threading.Lock()
        
        logger.info("APIKeyManager initialized")
    
    def create_key(self, request: APIKeyRequest) -> APIKeyResponse:
        """Create a new API key.
        
        Args:
            request: API key creation request
            
        Returns:
            APIKeyResponse with raw key (only shown once)
        """
        with self._lock:
            # Check key limit for owner
            owner_keys = self._keys_by_owner.get(request.owner_id, [])
            if len(owner_keys) >= request.tier.max_keys:
                raise ValueError(
                    f"Maximum keys ({request.tier.max_keys}) reached for tier {request.tier.value}"
                )
            
            # Generate key
            raw_key = self._generate_key(request.tier)
            key_hash = self._hash_key(raw_key)
            key_id = f"key_{uuid.uuid4().hex[:12]}"
            
            # Determine scopes
            scopes = request.scopes or request.tier.default_scopes
            
            # Create API key object
            api_key = APIKey(
                key_id=key_id,
                key_hash=key_hash,
                tier=request.tier,
                scopes=scopes,
                owner_id=request.owner_id,
                owner_email=request.owner_email,
                organization=request.organization,
                name=request.name,
                description=request.description,
                expires_at=(
                    datetime.utcnow() + timedelta(days=request.expires_in_days)
                    if request.expires_in_days
                    else None
                ),
            )
            
            # Store key
            self._keys[key_hash] = api_key
            self._keys_by_id[key_id] = key_hash
            
            if request.owner_id not in self._keys_by_owner:
                self._keys_by_owner[request.owner_id] = []
            self._keys_by_owner[request.owner_id].append(key_hash)
            
            logger.info(
                "API key created: id=%s, tier=%s, owner=%s",
                key_id, request.tier.value, request.owner_id[:8] + "..."
            )
            
            return APIKeyResponse(key=api_key, raw_key=raw_key)
    
    def validate_key(self, raw_key: str) -> Optional[APIKey]:
        """Validate an API key and return its data.
        
        Args:
            raw_key: The raw API key string
            
        Returns:
            APIKey if valid, None otherwise
        """
        # Quick format validation
        if not raw_key or not raw_key.startswith(f"{self.KEY_PREFIX}_"):
            return None
        
        key_hash = self._hash_key(raw_key)
        
        api_key = self._keys.get(key_hash)
        if not api_key:
            return None
        
        if not api_key.is_valid:
            return None
        
        # Update usage
        with self._lock:
            api_key.last_used_at = datetime.utcnow()
            api_key.request_count += 1
        
        return api_key
    
    def get_key(self, key_id: str) -> Optional[APIKey]:
        """Get API key by ID."""
        key_hash = self._keys_by_id.get(key_id)
        if key_hash:
            return self._keys.get(key_hash)
        return None
    
    def list_keys(
        self,
        owner_id: Optional[str] = None,
        tier: Optional[APIKeyTier] = None,
        include_revoked: bool = False
    ) -> List[APIKey]:
        """List API keys with optional filtering."""
        keys = []
        
        if owner_id:
            key_hashes = self._keys_by_owner.get(owner_id, [])
            keys = [self._keys[h] for h in key_hashes if h in self._keys]
        else:
            keys = list(self._keys.values())
        
        # Filter by tier
        if tier:
            keys = [k for k in keys if k.tier == tier]
        
        # Filter out revoked unless requested
        if not include_revoked:
            keys = [k for k in keys if not k.is_revoked]
        
        return keys
    
    def revoke_key(
        self,
        key_id: str,
        reason: Optional[str] = None
    ) -> Optional[APIKey]:
        """Revoke an API key.
        
        Args:
            key_id: Key ID to revoke
            reason: Optional reason for revocation
            
        Returns:
            Revoked APIKey or None if not found
        """
        api_key = self.get_key(key_id)
        if not api_key:
            return None
        
        with self._lock:
            api_key.is_revoked = True
            api_key.is_active = False
            api_key.revoked_at = datetime.utcnow()
            api_key.revoked_reason = reason
        
        logger.info("API key revoked: id=%s, reason=%s", key_id, reason)
        
        return api_key
    
    def rotate_key(self, key_id: str) -> Optional[APIKeyResponse]:
        """Rotate an API key (create new, revoke old).
        
        Args:
            key_id: Key ID to rotate
            
        Returns:
            New APIKeyResponse or None if not found
        """
        old_key = self.get_key(key_id)
        if not old_key:
            return None
        
        # Create new key with same properties
        request = APIKeyRequest(
            owner_id=old_key.owner_id,
            tier=old_key.tier,
            name=f"{old_key.name} (rotated)",
            description=old_key.description,
            scopes=old_key.scopes,
            organization=old_key.organization,
            owner_email=old_key.owner_email,
        )
        
        # Revoke old key
        self.revoke_key(key_id, reason="Rotated")
        
        # Create new key
        return self.create_key(request)
    
    def update_last_ip(self, key_id: str, ip: str):
        """Update last used IP for a key."""
        api_key = self.get_key(key_id)
        if api_key:
            api_key.last_ip = ip
    
    def get_usage_stats(self, key_id: str) -> Optional[Dict[str, Any]]:
        """Get usage statistics for a key."""
        api_key = self.get_key(key_id)
        if not api_key:
            return None
        
        return {
            "key_id": key_id,
            "tier": api_key.tier.value,
            "request_count": api_key.request_count,
            "last_used_at": api_key.last_used_at.isoformat() if api_key.last_used_at else None,
            "last_ip": api_key.last_ip,
            "created_at": api_key.created_at.isoformat(),
            "age_days": (datetime.utcnow() - api_key.created_at).days,
        }
    
    def _generate_key(self, tier: APIKeyTier) -> str:
        """Generate a secure API key.
        
        Format: nb_[tier]_[random base64]
        Example: nb_pro_xK7dR9sT2m...
        """
        tier_prefix = {
            APIKeyTier.FREE: "free",
            APIKeyTier.BASIC: "bas",
            APIKeyTier.PROFESSIONAL: "pro",
            APIKeyTier.ENTERPRISE: "ent",
            APIKeyTier.INTERNAL: "int",
        }.get(tier, "key")
        
        # Generate random bytes
        random_bytes = secrets.token_bytes(self.KEY_LENGTH)
        random_part = base64.urlsafe_b64encode(random_bytes).decode().rstrip("=")
        
        return f"{self.KEY_PREFIX}_{tier_prefix}_{random_part}"
    
    def _hash_key(self, raw_key: str) -> str:
        """Hash API key for secure storage."""
        return hashlib.sha256(raw_key.encode()).hexdigest()


# Singleton instance
_api_key_manager: Optional[APIKeyManager] = None


def get_api_key_manager() -> APIKeyManager:
    """Get or create the global API key manager instance."""
    global _api_key_manager
    
    if _api_key_manager is None:
        _api_key_manager = APIKeyManager()
    
    return _api_key_manager
