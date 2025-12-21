"""
OmniPath Client Service
Integration with OmniPath Manifest Infrastructure

Supports NeuroBotanica patent claims:
- Data provenance tracking for all database writes
- Traditional knowledge attribution
- Benefit sharing compliance
- Consent management

Reference: NeuroBotanica MVP Development Plan - Week 2
"""
from typing import Dict, Any, Optional, List
from dataclasses import dataclass, field
from datetime import datetime
from enum import Enum
import hashlib
import json
import logging
import os
from abc import ABC, abstractmethod

logger = logging.getLogger(__name__)


class ConsentLevel(Enum):
    """OmniPath consent levels."""
    FULL = "full"  # All data can be used and shared
    RESEARCH = "research_only"  # Research use only
    LIMITED = "limited"  # Anonymized aggregates only
    REVOKED = "revoked"  # No data use permitted


class OperationType(Enum):
    """Database operation types for provenance tracking."""
    CREATE = "CREATE"
    READ = "READ"
    UPDATE = "UPDATE"
    DELETE = "DELETE"
    AGGREGATE = "AGGREGATE"


@dataclass
class OmniPathToken:
    """OmniPath authentication token."""
    token_id: str
    user_id: str
    consent_level: ConsentLevel
    permissions: List[str]
    issued_at: datetime
    expires_at: datetime
    metadata: Dict[str, Any] = field(default_factory=dict)
    
    @property
    def is_expired(self) -> bool:
        return datetime.now() > self.expires_at
    
    @property
    def is_valid(self) -> bool:
        return not self.is_expired and self.consent_level != ConsentLevel.REVOKED
    
    def has_permission(self, permission: str) -> bool:
        return permission in self.permissions or "admin" in self.permissions


@dataclass
class ProvenanceManifest:
    """Data provenance manifest for OmniPath tracking."""
    manifest_id: str
    timestamp: datetime
    operation: OperationType
    data_hash: str
    data_type: str
    source: str
    user_id: Optional[str]
    traditional_knowledge_source: Optional[str]
    benefit_sharing_percentage: Optional[float]
    metadata: Dict[str, Any] = field(default_factory=dict)


class OmniPathClientBase(ABC):
    """Abstract base class for OmniPath client implementations."""
    
    @abstractmethod
    def validate_token(self, token: str) -> Optional[OmniPathToken]:
        """Validate an OmniPath token."""
        pass
    
    @abstractmethod
    def create_manifest(self, data: Dict[str, Any], operation: OperationType) -> ProvenanceManifest:
        """Create a provenance manifest for a data operation."""
        pass
    
    @abstractmethod
    def check_consent(self, user_id: str, data_type: str) -> ConsentLevel:
        """Check consent level for a user and data type."""
        pass
    
    @abstractmethod
    def record_traditional_knowledge_use(
        self,
        knowledge_source: str,
        usage_type: str,
        benefit_percentage: float
    ) -> str:
        """Record use of traditional knowledge for attribution."""
        pass


class OmniPathClient(OmniPathClientBase):
    """OmniPath integration client for NeuroBotanica.
    
    Handles:
    - Token validation and authentication
    - Consent management
    - Provenance tracking for all data operations
    - Traditional knowledge attribution
    - Benefit sharing calculations
    
    In MVP mode, this uses local storage. In production, it connects
    to the OmniPath distributed infrastructure.
    """
    
    def __init__(
        self,
        api_base_url: Optional[str] = None,
        api_key: Optional[str] = None,
        mode: str = "development"  # "development", "staging", "production"
    ):
        """Initialize OmniPath client.
        
        Args:
            api_base_url: OmniPath API base URL (uses env if not provided)
            api_key: OmniPath API key (uses env if not provided)
            mode: Operating mode (development uses local storage)
        """
        self.api_base_url = api_base_url or os.getenv("OMNIPATH_API_URL", "https://api.omnipath.io/v1")
        self.api_key = api_key or os.getenv("OMNIPATH_API_KEY", "")
        self.mode = mode
        
        # Local manifest storage for development
        self._local_manifests: Dict[str, ProvenanceManifest] = {}
        self._local_tokens: Dict[str, OmniPathToken] = {}
        self._consent_registry: Dict[str, ConsentLevel] = {}
        self._tk_usage_log: List[Dict] = []
        
        # Counter for generating manifest IDs
        self._manifest_counter = 0
        
        logger.info(f"OmniPath client initialized in {mode} mode")
        
        if mode == "development":
            self._setup_development_tokens()
    
    def _setup_development_tokens(self):
        """Set up development tokens for testing."""
        # Create a default development token
        dev_token = OmniPathToken(
            token_id="dev_token_001",
            user_id="dev_user",
            consent_level=ConsentLevel.FULL,
            permissions=["read", "write", "aggregate", "admin"],
            issued_at=datetime.now(),
            expires_at=datetime(2026, 12, 31),
            metadata={"environment": "development"}
        )
        self._local_tokens["dev_token_001"] = dev_token
        self._local_tokens["Bearer dev_token_001"] = dev_token
        
        # Research-only token
        research_token = OmniPathToken(
            token_id="research_token_001",
            user_id="researcher_user",
            consent_level=ConsentLevel.RESEARCH,
            permissions=["read", "aggregate"],
            issued_at=datetime.now(),
            expires_at=datetime(2026, 12, 31),
            metadata={"environment": "development", "institution": "University of Nevada"}
        )
        self._local_tokens["research_token_001"] = research_token
    
    def validate_token(self, token: str) -> Optional[OmniPathToken]:
        """Validate an OmniPath authentication token.
        
        Args:
            token: The token string to validate
            
        Returns:
            OmniPathToken if valid, None if invalid
        """
        if self.mode == "development":
            # Clean up token (remove "Bearer " prefix if present)
            clean_token = token.replace("Bearer ", "").strip()
            
            omnipath_token = self._local_tokens.get(clean_token)
            if omnipath_token and omnipath_token.is_valid:
                logger.debug(f"Token validated: {omnipath_token.token_id}")
                return omnipath_token
            return None
        
        # Production: Call OmniPath API
        # TODO: Implement actual API call
        logger.warning("Production token validation not yet implemented")
        return None
    
    def create_manifest(
        self,
        data: Dict[str, Any],
        operation: OperationType,
        user_id: Optional[str] = None,
        tk_source: Optional[str] = None,
        benefit_percentage: Optional[float] = None
    ) -> ProvenanceManifest:
        """Create a provenance manifest for a data operation.
        
        Args:
            data: The data being operated on
            operation: Type of operation (CREATE, UPDATE, etc.)
            user_id: ID of user performing operation
            tk_source: Traditional knowledge source if applicable
            benefit_percentage: Benefit sharing percentage if TK involved
            
        Returns:
            ProvenanceManifest with unique ID
        """
        self._manifest_counter += 1
        manifest_id = f"NB-MNF-{datetime.now().strftime('%Y%m%d')}-{self._manifest_counter:06d}"
        
        data_hash = self._hash_data(data)
        data_type = data.get("__type__", type(data).__name__)
        
        manifest = ProvenanceManifest(
            manifest_id=manifest_id,
            timestamp=datetime.now(),
            operation=operation,
            data_hash=data_hash,
            data_type=data_type,
            source="neurobotanica_mvp",
            user_id=user_id,
            traditional_knowledge_source=tk_source,
            benefit_sharing_percentage=benefit_percentage,
            metadata={
                "version": "1.0",
                "mode": self.mode,
                "data_size_bytes": len(json.dumps(data, default=str))
            }
        )
        
        if self.mode == "development":
            self._local_manifests[manifest_id] = manifest
            logger.info(f"Created manifest: {manifest_id} for {operation.value} on {data_type}")
        else:
            # TODO: Send to OmniPath API
            pass
        
        return manifest
    
    def get_manifest(self, manifest_id: str) -> Optional[ProvenanceManifest]:
        """Retrieve a provenance manifest by ID."""
        if self.mode == "development":
            return self._local_manifests.get(manifest_id)
        
        # TODO: Fetch from OmniPath API
        return None
    
    def check_consent(self, user_id: str, data_type: str) -> ConsentLevel:
        """Check consent level for a user and data type.
        
        Args:
            user_id: User identifier
            data_type: Type of data being accessed
            
        Returns:
            ConsentLevel for this user/data combination
        """
        if self.mode == "development":
            # In development, check local registry or return FULL
            key = f"{user_id}:{data_type}"
            return self._consent_registry.get(key, ConsentLevel.FULL)
        
        # TODO: Query OmniPath consent registry
        return ConsentLevel.LIMITED
    
    def set_consent(self, user_id: str, data_type: str, level: ConsentLevel):
        """Set consent level for a user and data type."""
        key = f"{user_id}:{data_type}"
        self._consent_registry[key] = level
        logger.info(f"Consent set: {user_id} -> {data_type} = {level.value}")
    
    def record_traditional_knowledge_use(
        self,
        knowledge_source: str,
        usage_type: str,
        benefit_percentage: float,
        compound_name: Optional[str] = None,
        condition: Optional[str] = None
    ) -> str:
        """Record use of traditional knowledge for attribution.
        
        Args:
            knowledge_source: Source of traditional knowledge (e.g., "Ayurvedic", "Indigenous North American")
            usage_type: How the knowledge is being used
            benefit_percentage: Percentage of benefits to be shared
            compound_name: Associated compound if applicable
            condition: Medical condition if applicable
            
        Returns:
            Attribution record ID
        """
        record_id = f"NB-TKA-{datetime.now().strftime('%Y%m%d%H%M%S')}"
        
        record = {
            "id": record_id,
            "timestamp": datetime.now().isoformat(),
            "knowledge_source": knowledge_source,
            "usage_type": usage_type,
            "benefit_percentage": benefit_percentage,
            "compound_name": compound_name,
            "condition": condition,
            "status": "active"
        }
        
        self._tk_usage_log.append(record)
        logger.info(f"TK usage recorded: {record_id} - {knowledge_source} ({benefit_percentage}%)")
        
        return record_id
    
    def get_tk_attribution_summary(self) -> Dict[str, Any]:
        """Get summary of traditional knowledge attributions."""
        if not self._tk_usage_log:
            return {"total_records": 0, "sources": {}, "total_benefit_pool": 0.0}
        
        sources = {}
        total_benefit = 0.0
        
        for record in self._tk_usage_log:
            source = record["knowledge_source"]
            if source not in sources:
                sources[source] = {"count": 0, "benefit_percentage": 0.0}
            sources[source]["count"] += 1
            sources[source]["benefit_percentage"] += record["benefit_percentage"]
            total_benefit += record["benefit_percentage"]
        
        return {
            "total_records": len(self._tk_usage_log),
            "sources": sources,
            "total_benefit_pool": total_benefit,
            "records": self._tk_usage_log
        }
    
    def _hash_data(self, data: Dict[str, Any]) -> str:
        """Generate SHA-256 hash of data for provenance."""
        # Convert to deterministic string representation
        data_str = json.dumps(data, sort_keys=True, default=str)
        return hashlib.sha256(data_str.encode()).hexdigest()
    
    def get_audit_trail(
        self,
        data_type: Optional[str] = None,
        user_id: Optional[str] = None,
        operation: Optional[OperationType] = None,
        start_date: Optional[datetime] = None,
        end_date: Optional[datetime] = None,
        limit: int = 100
    ) -> List[ProvenanceManifest]:
        """Query audit trail with filters.
        
        Args:
            data_type: Filter by data type
            user_id: Filter by user
            operation: Filter by operation type
            start_date: Filter from this date
            end_date: Filter until this date
            limit: Maximum records to return
            
        Returns:
            List of matching ProvenanceManifest records
        """
        results = []
        
        for manifest in self._local_manifests.values():
            # Apply filters
            if data_type and manifest.data_type != data_type:
                continue
            if user_id and manifest.user_id != user_id:
                continue
            if operation and manifest.operation != operation:
                continue
            if start_date and manifest.timestamp < start_date:
                continue
            if end_date and manifest.timestamp > end_date:
                continue
            
            results.append(manifest)
            
            if len(results) >= limit:
                break
        
        # Sort by timestamp descending
        results.sort(key=lambda m: m.timestamp, reverse=True)
        
        return results
    
    def verify_data_integrity(self, manifest_id: str, current_data: Dict[str, Any]) -> Dict[str, Any]:
        """Verify data integrity against stored manifest.
        
        Args:
            manifest_id: ID of the manifest to check against
            current_data: Current state of the data
            
        Returns:
            Verification result with integrity status
        """
        manifest = self.get_manifest(manifest_id)
        
        if not manifest:
            return {
                "verified": False,
                "error": f"Manifest {manifest_id} not found"
            }
        
        current_hash = self._hash_data(current_data)
        
        return {
            "verified": current_hash == manifest.data_hash,
            "manifest_id": manifest_id,
            "original_hash": manifest.data_hash,
            "current_hash": current_hash,
            "timestamp": manifest.timestamp.isoformat(),
            "operation": manifest.operation.value
        }


class OmniPathMockClient(OmniPathClientBase):
    """Mock OmniPath client for testing."""
    
    def validate_token(self, token: str) -> Optional[OmniPathToken]:
        """Always returns a valid token for testing."""
        return OmniPathToken(
            token_id="mock_token",
            user_id="mock_user",
            consent_level=ConsentLevel.FULL,
            permissions=["read", "write", "aggregate"],
            issued_at=datetime.now(),
            expires_at=datetime(2030, 12, 31)
        )
    
    def create_manifest(self, data: Dict[str, Any], operation: OperationType) -> ProvenanceManifest:
        """Create a mock manifest."""
        return ProvenanceManifest(
            manifest_id="MOCK-001",
            timestamp=datetime.now(),
            operation=operation,
            data_hash="mock_hash",
            data_type="mock",
            source="mock",
            user_id="mock_user",
            traditional_knowledge_source=None,
            benefit_sharing_percentage=None
        )
    
    def check_consent(self, user_id: str, data_type: str) -> ConsentLevel:
        """Always return FULL consent for testing."""
        return ConsentLevel.FULL
    
    def record_traditional_knowledge_use(
        self,
        knowledge_source: str,
        usage_type: str,
        benefit_percentage: float
    ) -> str:
        """Return a mock attribution ID."""
        return "MOCK-TKA-001"


# Global client instance
_omnipath_client: Optional[OmniPathClient] = None


def get_omnipath_client() -> OmniPathClient:
    """Get or create the global OmniPath client instance."""
    global _omnipath_client
    
    if _omnipath_client is None:
        mode = os.getenv("OMNIPATH_MODE", "development")
        _omnipath_client = OmniPathClient(mode=mode)
    
    return _omnipath_client


def reset_omnipath_client():
    """Reset the global client (useful for testing)."""
    global _omnipath_client
    _omnipath_client = None
