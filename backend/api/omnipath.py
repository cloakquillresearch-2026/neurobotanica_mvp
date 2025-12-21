"""
OmniPath API Endpoints
Provides REST API for OmniPath integration features

Supports NeuroBotanica patent claims:
- Consent management endpoints
- Audit trail query endpoints
- Traditional knowledge attribution
- Provenance verification

Reference: NeuroBotanica MVP Development Plan - Week 2 Task 2.2
"""
from typing import Optional, List, Dict, Any
from datetime import datetime, timedelta
from pydantic import BaseModel, Field
from fastapi import APIRouter, HTTPException, status, Depends, Query, Request

from ..services.omnipath_client import (
    get_omnipath_client,
    OmniPathToken,
    ConsentLevel,
    OperationType
)
from ..services.provenance_tracker import get_provenance_tracker
from ..middleware.token_validation import (
    get_current_user,
    get_admin_user,
    require_permissions,
    rate_limit
)

router = APIRouter(prefix="/omnipath", tags=["OmniPath"])


# ==================== Request/Response Models ====================

class ConsentRequest(BaseModel):
    """Request to update consent level."""
    consent_level: str = Field(..., description="Consent level: FULL, RESEARCH, LIMITED, REVOKED")
    data_types: Optional[List[str]] = Field(
        default=None,
        description="Specific data types to apply consent to"
    )
    
    class Config:
        json_schema_extra = {
            "example": {
                "consent_level": "RESEARCH",
                "data_types": ["treatment_data", "terpene_analysis"]
            }
        }


class ConsentResponse(BaseModel):
    """Response for consent operations."""
    user_id: str
    consent_level: str
    data_types: List[str]
    updated_at: str
    message: str


class AuditQueryRequest(BaseModel):
    """Request for querying audit trail."""
    table_name: Optional[str] = None
    record_id: Optional[int] = None
    user_id: Optional[str] = None
    operation: Optional[str] = None
    start_date: Optional[datetime] = None
    end_date: Optional[datetime] = None
    limit: int = Field(default=100, le=1000)


class AuditEntry(BaseModel):
    """Single audit trail entry."""
    entry_id: str
    timestamp: str
    operation: str
    table_name: str
    record_id: Optional[int]
    user_id: Optional[str]
    manifest_id: Optional[str]
    ip_address: Optional[str]


class AuditTrailResponse(BaseModel):
    """Response for audit trail queries."""
    total_entries: int
    entries: List[AuditEntry]
    query_time_ms: float


class TKUsageRecord(BaseModel):
    """Traditional knowledge usage record."""
    knowledge_source: str
    usage_type: str
    benefit_percentage: float
    compound_name: Optional[str]
    timestamp: str


class TKAttributionResponse(BaseModel):
    """Response for TK attribution queries."""
    total_usages: int
    total_benefit_value: float
    records: List[TKUsageRecord]


class ProvenanceVerifyRequest(BaseModel):
    """Request to verify data provenance."""
    table_name: str
    record_id: int
    current_data: Dict[str, Any]


class ProvenanceVerifyResponse(BaseModel):
    """Response for provenance verification."""
    verified: bool
    record_id: int
    table_name: str
    history_entries: int
    error: Optional[str] = None
    manifest_id: Optional[str] = None


class TokenInfoResponse(BaseModel):
    """Response with token information."""
    user_id: str
    permissions: List[str]
    expires_at: str
    consent_level: str


# ==================== Consent Management Endpoints ====================

@router.get("/consent", response_model=ConsentResponse)
async def get_consent_status(
    user: OmniPathToken = Depends(get_current_user)
):
    """Get current user's consent status.
    
    Returns the user's consent level and applicable data types.
    """
    client = get_omnipath_client()
    
    # Get consent level for all data types
    consent_level = client.check_consent(user.user_id, "all")
    
    return ConsentResponse(
        user_id=user.user_id,
        consent_level=consent_level.value,
        data_types=["treatment_data", "terpene_analysis", "research_data"],
        updated_at=datetime.now().isoformat(),
        message="Consent status retrieved successfully"
    )


@router.put("/consent", response_model=ConsentResponse)
async def update_consent(
    request: ConsentRequest,
    user: OmniPathToken = Depends(get_current_user)
):
    """Update user's consent level.
    
    Allows users to modify their consent preferences for data usage.
    Changes are tracked in the audit trail.
    """
    client = get_omnipath_client()
    tracker = get_provenance_tracker()
    
    # Validate consent level
    try:
        consent_level = ConsentLevel[request.consent_level.upper()]
    except KeyError:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=f"Invalid consent level. Must be one of: {[c.value for c in ConsentLevel]}"
        )
    
    # Update consent in OmniPath client
    data_types = request.data_types or ["all"]
    for data_type in data_types:
        client._consent_registry[f"{user.user_id}:{data_type}"] = consent_level
    
    # Track the consent change
    tracker.track_update(
        record_id=0,  # Consent doesn't have a record ID
        data_before={"consent_level": "previous"},  # Would need to track this
        data_after={"consent_level": consent_level.value, "data_types": data_types},
        table_name="consent_registry",
        user_id=user.user_id
    )
    
    return ConsentResponse(
        user_id=user.user_id,
        consent_level=consent_level.value,
        data_types=data_types,
        updated_at=datetime.now().isoformat(),
        message="Consent updated successfully"
    )


@router.delete("/consent")
async def revoke_consent(
    user: OmniPathToken = Depends(get_current_user)
):
    """Revoke all consent.
    
    Sets consent level to REVOKED for all data types.
    This action is logged and cannot be automatically reversed.
    """
    client = get_omnipath_client()
    tracker = get_provenance_tracker()
    
    # Revoke all consent
    for key in list(client._consent_registry.keys()):
        if key.startswith(f"{user.user_id}:"):
            client._consent_registry[key] = ConsentLevel.REVOKED
    
    # Track revocation
    tracker.track_delete(
        record_id=0,
        data={"user_id": user.user_id, "action": "revoke_all"},
        table_name="consent_registry",
        user_id=user.user_id
    )
    
    return {
        "message": "All consent revoked",
        "user_id": user.user_id,
        "revoked_at": datetime.now().isoformat()
    }


# ==================== Audit Trail Endpoints ====================

@router.get("/audit", response_model=AuditTrailResponse)
async def get_audit_trail(
    table_name: Optional[str] = Query(None, description="Filter by table name"),
    record_id: Optional[int] = Query(None, description="Filter by record ID"),
    operation: Optional[str] = Query(None, description="Filter by operation type"),
    start_date: Optional[datetime] = Query(None, description="Start date filter"),
    end_date: Optional[datetime] = Query(None, description="End date filter"),
    limit: int = Query(100, le=1000, description="Maximum results"),
    user: OmniPathToken = Depends(get_current_user)
):
    """Query the audit trail.
    
    Returns audit entries based on provided filters.
    Admin users can see all entries, regular users only see their own.
    """
    import time
    start_time = time.time()
    
    tracker = get_provenance_tracker()
    
    # Regular users can only see their own audit entries
    query_user_id = None
    if "admin:all" not in user.permissions:
        query_user_id = user.user_id
    
    entries = tracker.get_audit_log(
        table_name=table_name,
        record_id=record_id,
        user_id=query_user_id,
        operation=operation,
        start_date=start_date,
        end_date=end_date,
        limit=limit
    )
    
    query_time = (time.time() - start_time) * 1000
    
    return AuditTrailResponse(
        total_entries=len(entries),
        entries=[AuditEntry(**{k: v for k, v in e.items() if k in AuditEntry.model_fields}) for e in entries],
        query_time_ms=round(query_time, 2)
    )


@router.get("/audit/record/{table_name}/{record_id}")
async def get_record_history(
    table_name: str,
    record_id: int,
    user: OmniPathToken = Depends(get_current_user)
):
    """Get complete history for a specific record.
    
    Returns all operations performed on the specified record in chronological order.
    """
    tracker = get_provenance_tracker()
    
    history = tracker.get_record_history(table_name, record_id)
    
    return {
        "table_name": table_name,
        "record_id": record_id,
        "total_operations": len(history),
        "history": history
    }


@router.get("/audit/stats")
async def get_audit_stats(
    days: int = Query(7, description="Number of days to analyze"),
    user: OmniPathToken = Depends(get_admin_user)
):
    """Get audit trail statistics (admin only).
    
    Returns aggregated statistics about data operations.
    """
    tracker = get_provenance_tracker()
    
    end_date = datetime.now()
    start_date = end_date - timedelta(days=days)
    
    entries = tracker.get_audit_log(
        start_date=start_date,
        end_date=end_date,
        limit=10000
    )
    
    # Calculate stats
    operations_by_type = {}
    operations_by_table = {}
    operations_by_user = {}
    
    for entry in entries:
        op = entry.get("operation", "unknown")
        table = entry.get("table_name", "unknown")
        uid = entry.get("user_id", "unknown")
        
        operations_by_type[op] = operations_by_type.get(op, 0) + 1
        operations_by_table[table] = operations_by_table.get(table, 0) + 1
        operations_by_user[uid] = operations_by_user.get(uid, 0) + 1
    
    return {
        "period": {
            "start": start_date.isoformat(),
            "end": end_date.isoformat(),
            "days": days
        },
        "total_operations": len(entries),
        "operations_by_type": operations_by_type,
        "operations_by_table": operations_by_table,
        "unique_users": len(operations_by_user),
        "top_users": sorted(
            operations_by_user.items(),
            key=lambda x: x[1],
            reverse=True
        )[:10]
    }


# ==================== Traditional Knowledge Attribution ====================

@router.get("/tk/attribution", response_model=TKAttributionResponse)
async def get_tk_attribution(
    knowledge_source: Optional[str] = Query(None, description="Filter by knowledge source"),
    compound_name: Optional[str] = Query(None, description="Filter by compound name"),
    limit: int = Query(100, le=500, description="Maximum results"),
    user: OmniPathToken = Depends(get_current_user)
):
    """Get traditional knowledge attribution records.
    
    Returns records of traditional knowledge usage and benefit sharing.
    """
    client = get_omnipath_client()
    
    records = []
    total_benefit = 0.0
    
    for record in client._tk_usage_log:
        # Apply filters
        if knowledge_source and record.get("knowledge_source") != knowledge_source:
            continue
        if compound_name and record.get("compound_name") != compound_name:
            continue
        
        records.append(TKUsageRecord(
            knowledge_source=record.get("knowledge_source", ""),
            usage_type=record.get("usage_type", ""),
            benefit_percentage=record.get("benefit_percentage", 0.0),
            compound_name=record.get("compound_name"),
            timestamp=record.get("timestamp", "")
        ))
        
        total_benefit += record.get("benefit_percentage", 0.0)
        
        if len(records) >= limit:
            break
    
    return TKAttributionResponse(
        total_usages=len(records),
        total_benefit_value=total_benefit,
        records=records
    )


@router.post("/tk/record")
async def record_tk_usage(
    knowledge_source: str,
    usage_type: str,
    benefit_percentage: float = Query(..., ge=0, le=100),
    compound_name: Optional[str] = None,
    user: OmniPathToken = Depends(get_current_user)
):
    """Record a traditional knowledge usage.
    
    Tracks usage of traditional knowledge for benefit sharing compliance.
    """
    if "write:tk_records" not in user.permissions and "admin:all" not in user.permissions:
        raise HTTPException(
            status_code=status.HTTP_403_FORBIDDEN,
            detail="Missing permission: write:tk_records"
        )
    
    client = get_omnipath_client()
    
    record = client.record_traditional_knowledge_use(
        knowledge_source=knowledge_source,
        usage_type=usage_type,
        benefit_percentage=benefit_percentage,
        compound_name=compound_name
    )
    
    return {
        "message": "TK usage recorded successfully",
        "record": record
    }


@router.get("/tk/sources")
async def get_tk_sources(
    user: OmniPathToken = Depends(get_current_user)
):
    """Get list of all traditional knowledge sources.
    
    Returns unique knowledge sources from attribution records.
    """
    client = get_omnipath_client()
    
    sources = set()
    source_stats = {}
    
    for record in client._tk_usage_log:
        source = record.get("knowledge_source", "")
        sources.add(source)
        
        if source not in source_stats:
            source_stats[source] = {"usages": 0, "total_benefit": 0.0}
        
        source_stats[source]["usages"] += 1
        source_stats[source]["total_benefit"] += record.get("benefit_percentage", 0.0)
    
    return {
        "total_sources": len(sources),
        "sources": sorted(list(sources)),
        "source_statistics": source_stats
    }


# ==================== Provenance Verification ====================

@router.post("/verify", response_model=ProvenanceVerifyResponse)
async def verify_provenance(
    request: ProvenanceVerifyRequest,
    user: OmniPathToken = Depends(get_current_user)
):
    """Verify data integrity against provenance records.
    
    Checks if current data matches the provenance trail.
    """
    tracker = get_provenance_tracker()
    
    result = tracker.verify_record_integrity(
        table_name=request.table_name,
        record_id=request.record_id,
        current_data=request.current_data
    )
    
    return ProvenanceVerifyResponse(**result)


@router.get("/manifest/{manifest_id}")
async def get_manifest(
    manifest_id: str,
    user: OmniPathToken = Depends(get_current_user)
):
    """Get a specific provenance manifest by ID.
    
    Returns the full manifest details for verification.
    """
    client = get_omnipath_client()
    
    manifest = client._local_manifests.get(manifest_id)
    
    if not manifest:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"Manifest not found: {manifest_id}"
        )
    
    return {
        "manifest_id": manifest.manifest_id,
        "data_hash": manifest.data_hash,
        "timestamp": manifest.timestamp.isoformat(),
        "operation": manifest.operation.value,
        "user_id": manifest.user_id,
        "signature": manifest.signature
    }


# ==================== Token Information ====================

@router.get("/token/info", response_model=TokenInfoResponse)
async def get_token_info(
    user: OmniPathToken = Depends(get_current_user)
):
    """Get information about the current token.
    
    Returns token details including permissions and expiry.
    """
    client = get_omnipath_client()
    consent_level = client.check_consent(user.user_id, "all")
    
    return TokenInfoResponse(
        user_id=user.user_id,
        permissions=user.permissions,
        expires_at=user.expires_at.isoformat(),
        consent_level=consent_level.value
    )


@router.get("/token/permissions")
async def list_available_permissions(
    user: OmniPathToken = Depends(get_current_user)
):
    """List all available permission types.
    
    Returns the full list of permissions that can be assigned.
    """
    permissions = {
        "read:cannabinoids": "Read cannabinoid data",
        "write:cannabinoids": "Create/update cannabinoid data",
        "delete:cannabinoids": "Delete cannabinoid data",
        "read:patients": "Read patient data",
        "write:patients": "Create/update patient data",
        "read:treatments": "Read treatment data",
        "write:treatments": "Create/update treatment data",
        "read:research_data": "Access research data",
        "read:terpene_analysis": "Access terpene analysis results",
        "write:tk_records": "Record traditional knowledge usage",
        "admin:all": "Full administrative access"
    }
    
    return {
        "available_permissions": permissions,
        "your_permissions": user.permissions
    }


# ==================== Health Check ====================

@router.get("/health")
async def omnipath_health():
    """Check OmniPath integration health.
    
    Verifies connectivity and performance of OmniPath services.
    """
    import time
    
    client = get_omnipath_client()
    tracker = get_provenance_tracker()
    
    # Test token validation speed
    start = time.time()
    _ = client.validate_token("test_token")
    token_time = (time.time() - start) * 1000
    
    # Test manifest creation speed
    start = time.time()
    _ = client.create_manifest(
        data={"test": "data"},
        operation=OperationType.READ
    )
    manifest_time = (time.time() - start) * 1000
    
    return {
        "status": "healthy",
        "mode": client.mode,
        "metrics": {
            "token_validation_ms": round(token_time, 2),
            "manifest_creation_ms": round(manifest_time, 2),
            "token_validation_target_met": token_time < 100,
            "local_manifests_count": len(client._local_manifests),
            "audit_entries_count": len(tracker._audit_log),
            "tk_records_count": len(client._tk_usage_log)
        }
    }
