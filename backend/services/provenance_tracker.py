"""
Provenance Tracker Service
Track data provenance for OmniPath manifest system

Supports NeuroBotanica patent claims:
- Immutable audit trail for all data operations
- Traditional knowledge attribution tracking
- Benefit sharing compliance
- Data integrity verification

Reference: NeuroBotanica MVP Development Plan - Week 2 Task 2.3
"""
from typing import Dict, Any, Optional, List, Callable
from dataclasses import dataclass, asdict
from datetime import datetime
from functools import wraps
import hashlib
import json
import logging
from contextlib import contextmanager

from sqlalchemy.orm import Session
from sqlalchemy import event

from .omnipath_client import (
    OmniPathClient,
    get_omnipath_client,
    OperationType,
    ProvenanceManifest,
    ConsentLevel
)

logger = logging.getLogger(__name__)


@dataclass
class AuditEntry:
    """Single audit trail entry."""
    entry_id: str
    timestamp: datetime
    operation: str
    table_name: str
    record_id: Optional[int]
    user_id: Optional[str]
    data_before: Optional[Dict[str, Any]]
    data_after: Optional[Dict[str, Any]]
    manifest_id: Optional[str]
    ip_address: Optional[str]
    user_agent: Optional[str]


class ProvenanceTracker:
    """Track data provenance for OmniPath manifest system.
    
    Provides:
    - Automatic manifest creation for all database writes
    - Audit trail maintenance
    - Data integrity verification
    - Traditional knowledge attribution
    - Consent checking before data access
    
    Usage:
        tracker = ProvenanceTracker()
        
        # Track a create operation
        manifest = tracker.track_create(
            data={"name": "CBD", "smiles": "..."},
            table_name="cannabinoids",
            user_id="user_123"
        )
        
        # Track an update
        manifest = tracker.track_update(
            record_id=1,
            data_before=old_data,
            data_after=new_data,
            table_name="cannabinoids",
            user_id="user_123"
        )
    """
    
    def __init__(self, omnipath_client: Optional[OmniPathClient] = None):
        """Initialize provenance tracker.
        
        Args:
            omnipath_client: OmniPath client instance (creates default if None)
        """
        self.omnipath_client = omnipath_client or get_omnipath_client()
        self._audit_log: List[AuditEntry] = []
        self._entry_counter = 0
        
        logger.info("ProvenanceTracker initialized")
    
    def track_create(
        self,
        data: Dict[str, Any],
        table_name: str,
        user_id: Optional[str] = None,
        record_id: Optional[int] = None,
        tk_source: Optional[str] = None,
        benefit_percentage: Optional[float] = None,
        request_context: Optional[Dict] = None
    ) -> ProvenanceManifest:
        """Track a CREATE operation.
        
        Args:
            data: The data being created
            table_name: Database table name
            user_id: User performing the operation
            record_id: ID of created record (if known)
            tk_source: Traditional knowledge source if applicable
            benefit_percentage: Benefit sharing percentage
            request_context: HTTP request context (IP, user agent, etc.)
            
        Returns:
            ProvenanceManifest for the operation
        """
        # Add type hint for provenance
        data_with_type = {**data, "__type__": table_name}
        
        # Create manifest
        manifest = self.omnipath_client.create_manifest(
            data=data_with_type,
            operation=OperationType.CREATE,
            user_id=user_id,
            tk_source=tk_source,
            benefit_percentage=benefit_percentage
        )
        
        # Log TK usage if applicable
        if tk_source and benefit_percentage:
            self.omnipath_client.record_traditional_knowledge_use(
                knowledge_source=tk_source,
                usage_type=f"data_creation:{table_name}",
                benefit_percentage=benefit_percentage,
                compound_name=data.get("name")
            )
        
        # Create audit entry
        self._create_audit_entry(
            operation="CREATE",
            table_name=table_name,
            record_id=record_id,
            user_id=user_id,
            data_before=None,
            data_after=data,
            manifest_id=manifest.manifest_id,
            request_context=request_context
        )
        
        return manifest
    
    def track_update(
        self,
        record_id: int,
        data_before: Dict[str, Any],
        data_after: Dict[str, Any],
        table_name: str,
        user_id: Optional[str] = None,
        request_context: Optional[Dict] = None
    ) -> ProvenanceManifest:
        """Track an UPDATE operation.
        
        Args:
            record_id: ID of the record being updated
            data_before: State before update
            data_after: State after update
            table_name: Database table name
            user_id: User performing the operation
            request_context: HTTP request context
            
        Returns:
            ProvenanceManifest for the operation
        """
        # Create manifest with change delta
        change_data = {
            "__type__": table_name,
            "record_id": record_id,
            "changes": self._calculate_diff(data_before, data_after)
        }
        
        manifest = self.omnipath_client.create_manifest(
            data=change_data,
            operation=OperationType.UPDATE,
            user_id=user_id
        )
        
        # Create audit entry
        self._create_audit_entry(
            operation="UPDATE",
            table_name=table_name,
            record_id=record_id,
            user_id=user_id,
            data_before=data_before,
            data_after=data_after,
            manifest_id=manifest.manifest_id,
            request_context=request_context
        )
        
        return manifest
    
    def track_delete(
        self,
        record_id: int,
        data: Dict[str, Any],
        table_name: str,
        user_id: Optional[str] = None,
        request_context: Optional[Dict] = None
    ) -> ProvenanceManifest:
        """Track a DELETE operation.
        
        Args:
            record_id: ID of the record being deleted
            data: The data being deleted (for audit trail)
            table_name: Database table name
            user_id: User performing the operation
            request_context: HTTP request context
            
        Returns:
            ProvenanceManifest for the operation
        """
        data_with_type = {**data, "__type__": table_name, "record_id": record_id}
        
        manifest = self.omnipath_client.create_manifest(
            data=data_with_type,
            operation=OperationType.DELETE,
            user_id=user_id
        )
        
        # Create audit entry
        self._create_audit_entry(
            operation="DELETE",
            table_name=table_name,
            record_id=record_id,
            user_id=user_id,
            data_before=data,
            data_after=None,
            manifest_id=manifest.manifest_id,
            request_context=request_context
        )
        
        return manifest
    
    def track_read(
        self,
        record_ids: List[int],
        table_name: str,
        user_id: Optional[str] = None,
        query_params: Optional[Dict] = None
    ) -> ProvenanceManifest:
        """Track a READ operation (for sensitive data access logging).
        
        Args:
            record_ids: IDs of records being accessed
            table_name: Database table name
            user_id: User performing the operation
            query_params: Query parameters used
            
        Returns:
            ProvenanceManifest for the operation
        """
        read_data = {
            "__type__": f"{table_name}_read",
            "record_ids": record_ids,
            "query_params": query_params or {}
        }
        
        manifest = self.omnipath_client.create_manifest(
            data=read_data,
            operation=OperationType.READ,
            user_id=user_id
        )
        
        return manifest
    
    def check_consent_for_access(
        self,
        user_id: str,
        table_name: str,
        operation: OperationType
    ) -> tuple[bool, ConsentLevel]:
        """Check if user has consent to perform operation.
        
        Args:
            user_id: User requesting access
            table_name: Table being accessed
            operation: Type of operation
            
        Returns:
            Tuple of (allowed, consent_level)
        """
        consent_level = self.omnipath_client.check_consent(user_id, table_name)
        
        # Define access rules based on consent level
        allowed = False
        
        if consent_level == ConsentLevel.FULL:
            allowed = True
        elif consent_level == ConsentLevel.RESEARCH:
            allowed = operation in [OperationType.READ, OperationType.AGGREGATE]
        elif consent_level == ConsentLevel.LIMITED:
            allowed = operation == OperationType.AGGREGATE
        elif consent_level == ConsentLevel.REVOKED:
            allowed = False
        
        if not allowed:
            logger.warning(
                f"Access denied: user={user_id}, table={table_name}, "
                f"operation={operation.value}, consent={consent_level.value}"
            )
        
        return allowed, consent_level
    
    def _create_audit_entry(
        self,
        operation: str,
        table_name: str,
        record_id: Optional[int],
        user_id: Optional[str],
        data_before: Optional[Dict],
        data_after: Optional[Dict],
        manifest_id: str,
        request_context: Optional[Dict] = None
    ):
        """Create an audit trail entry."""
        self._entry_counter += 1
        entry_id = f"AUD-{datetime.now().strftime('%Y%m%d%H%M%S')}-{self._entry_counter:06d}"
        
        entry = AuditEntry(
            entry_id=entry_id,
            timestamp=datetime.now(),
            operation=operation,
            table_name=table_name,
            record_id=record_id,
            user_id=user_id,
            data_before=data_before,
            data_after=data_after,
            manifest_id=manifest_id,
            ip_address=request_context.get("ip") if request_context else None,
            user_agent=request_context.get("user_agent") if request_context else None
        )
        
        self._audit_log.append(entry)
        logger.debug(f"Audit entry created: {entry_id}")
    
    def _calculate_diff(
        self,
        before: Dict[str, Any],
        after: Dict[str, Any]
    ) -> Dict[str, Any]:
        """Calculate the difference between two data states."""
        changes = {}
        
        all_keys = set(before.keys()) | set(after.keys())
        
        for key in all_keys:
            old_val = before.get(key)
            new_val = after.get(key)
            
            if old_val != new_val:
                changes[key] = {
                    "old": old_val,
                    "new": new_val
                }
        
        return changes
    
    def get_audit_log(
        self,
        table_name: Optional[str] = None,
        record_id: Optional[int] = None,
        user_id: Optional[str] = None,
        operation: Optional[str] = None,
        start_date: Optional[datetime] = None,
        end_date: Optional[datetime] = None,
        limit: int = 100
    ) -> List[Dict[str, Any]]:
        """Query the audit log with filters.
        
        Returns:
            List of audit entries as dictionaries
        """
        results = []
        
        for entry in self._audit_log:
            if table_name and entry.table_name != table_name:
                continue
            if record_id and entry.record_id != record_id:
                continue
            if user_id and entry.user_id != user_id:
                continue
            if operation and entry.operation != operation:
                continue
            if start_date and entry.timestamp < start_date:
                continue
            if end_date and entry.timestamp > end_date:
                continue
            
            results.append(asdict(entry))
            
            if len(results) >= limit:
                break
        
        # Sort by timestamp descending
        results.sort(key=lambda x: x["timestamp"], reverse=True)
        
        # Convert datetime objects to ISO strings
        for result in results:
            result["timestamp"] = result["timestamp"].isoformat()
        
        return results
    
    def get_record_history(
        self,
        table_name: str,
        record_id: int
    ) -> List[Dict[str, Any]]:
        """Get complete history for a specific record.
        
        Args:
            table_name: Database table name
            record_id: Record ID
            
        Returns:
            List of all operations on this record in chronological order
        """
        history = self.get_audit_log(
            table_name=table_name,
            record_id=record_id,
            limit=1000
        )
        
        # Return in chronological order
        history.reverse()
        
        return history
    
    def verify_record_integrity(
        self,
        table_name: str,
        record_id: int,
        current_data: Dict[str, Any]
    ) -> Dict[str, Any]:
        """Verify integrity of a record against its provenance trail.
        
        Args:
            table_name: Database table name
            record_id: Record ID
            current_data: Current state of the record
            
        Returns:
            Verification result
        """
        history = self.get_record_history(table_name, record_id)
        
        if not history:
            return {
                "verified": False,
                "error": "No provenance history found",
                "record_id": record_id
            }
        
        # Get the most recent manifest
        latest_entry = history[-1]
        manifest_id = latest_entry.get("manifest_id")
        
        if not manifest_id:
            return {
                "verified": False,
                "error": "No manifest ID in history",
                "record_id": record_id
            }
        
        # Verify against OmniPath manifest
        verification = self.omnipath_client.verify_data_integrity(
            manifest_id=manifest_id,
            current_data=current_data
        )
        
        verification["record_id"] = record_id
        verification["table_name"] = table_name
        verification["history_entries"] = len(history)
        
        return verification


# SQLAlchemy event listeners for automatic provenance tracking
def setup_provenance_listeners(base_class, tracker: ProvenanceTracker):
    """Set up SQLAlchemy event listeners for automatic provenance tracking.
    
    Args:
        base_class: SQLAlchemy declarative base
        tracker: ProvenanceTracker instance
    """
    
    @event.listens_for(base_class, "after_insert", propagate=True)
    def after_insert(mapper, connection, target):
        """Track INSERT operations."""
        table_name = target.__tablename__
        data = {c.name: getattr(target, c.name) for c in target.__table__.columns}
        record_id = getattr(target, "id", None)
        
        tracker.track_create(
            data=data,
            table_name=table_name,
            record_id=record_id
        )
    
    @event.listens_for(base_class, "after_update", propagate=True)
    def after_update(mapper, connection, target):
        """Track UPDATE operations."""
        table_name = target.__tablename__
        
        # Get current (after) state
        data_after = {c.name: getattr(target, c.name) for c in target.__table__.columns}
        record_id = getattr(target, "id", None)
        
        # Note: Getting "before" state requires additional setup
        # For now, we track the current state
        tracker.track_update(
            record_id=record_id,
            data_before={},  # Would need history tracking
            data_after=data_after,
            table_name=table_name
        )
    
    @event.listens_for(base_class, "after_delete", propagate=True)
    def after_delete(mapper, connection, target):
        """Track DELETE operations."""
        table_name = target.__tablename__
        data = {c.name: getattr(target, c.name) for c in target.__table__.columns}
        record_id = getattr(target, "id", None)
        
        tracker.track_delete(
            record_id=record_id,
            data=data,
            table_name=table_name
        )
    
    logger.info("Provenance event listeners registered")


# Decorator for tracking function calls
def track_provenance(table_name: str, operation: OperationType):
    """Decorator to track provenance for a function.
    
    Usage:
        @track_provenance("cannabinoids", OperationType.CREATE)
        def create_compound(data):
            ...
    """
    def decorator(func: Callable):
        @wraps(func)
        def wrapper(*args, **kwargs):
            tracker = get_provenance_tracker()
            
            # Execute function
            result = func(*args, **kwargs)
            
            # Track based on operation type
            if operation == OperationType.CREATE:
                tracker.track_create(
                    data=kwargs.get("data", {}),
                    table_name=table_name
                )
            
            return result
        return wrapper
    return decorator


# Context manager for tracked transactions
@contextmanager
def tracked_transaction(
    tracker: ProvenanceTracker,
    table_name: str,
    operation: OperationType,
    user_id: Optional[str] = None
):
    """Context manager for tracking database transactions.
    
    Usage:
        with tracked_transaction(tracker, "cannabinoids", OperationType.CREATE) as ctx:
            # Perform database operations
            ctx["data"] = {...}
    """
    context = {
        "data": None,
        "record_id": None,
        "manifest": None
    }
    
    try:
        yield context
        
        # Track after successful completion
        if context["data"]:
            if operation == OperationType.CREATE:
                context["manifest"] = tracker.track_create(
                    data=context["data"],
                    table_name=table_name,
                    record_id=context.get("record_id"),
                    user_id=user_id
                )
    except Exception as e:
        logger.error(f"Transaction failed: {e}")
        raise


# Global tracker instance
_provenance_tracker: Optional[ProvenanceTracker] = None


def get_provenance_tracker() -> ProvenanceTracker:
    """Get or create the global provenance tracker instance."""
    global _provenance_tracker
    
    if _provenance_tracker is None:
        _provenance_tracker = ProvenanceTracker()
    
    return _provenance_tracker


def reset_provenance_tracker():
    """Reset the global tracker (useful for testing)."""
    global _provenance_tracker
    _provenance_tracker = None
