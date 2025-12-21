"""
MetaPath API Router

REST API endpoints for 31-Pathway Orchestration, DAO Governance,
Emergency Coordination, and Transparency Dashboard.

Trade Secret: API design for coordinated multi-pathway operations.
"""

from fastapi import APIRouter, HTTPException, Query, Path, Body
from typing import Dict, List, Optional, Any
from pydantic import BaseModel, Field
from datetime import datetime

from backend.services.metapath import (
    MetaPathOrchestrator,
    DAOGovernance,
    EmergencyCoordinator,
    TransparencyDashboard,
    PathwayID,
    PathwayLayer,
    PriorityLevel,
    WorkflowStatus,
    MemberRole,
    ProposalType,
    ProposalStatus,
    CrisisLevel,
    CrisisType,
    MetricCategory,
    PrivacyLevel,
    AlertSeverity,
)


router = APIRouter(prefix="/api/v1/metapath", tags=["metapath"])


# Initialize services
_orchestrator = MetaPathOrchestrator()
_governance = DAOGovernance()
_emergency = EmergencyCoordinator()
_dashboard = TransparencyDashboard()


# ==============================================================================
# Request/Response Models
# ==============================================================================

class WorkflowStepRequest(BaseModel):
    """Request model for workflow step."""
    pathway: str
    operation: str
    parameters: Dict[str, Any] = Field(default_factory=dict)
    depends_on: List[str] = Field(default_factory=list)
    involves_tk: bool = False
    tk_communities: List[str] = Field(default_factory=list)


class CreateWorkflowRequest(BaseModel):
    """Request to create workflow."""
    name: str
    description: str = ""
    steps: List[WorkflowStepRequest]
    priority: str = "normal"
    involves_tk: bool = False
    created_by: str = "api_user"
    dao_proposal_id: Optional[str] = None


class WorkflowResponse(BaseModel):
    """Workflow response."""
    workflow_id: str
    name: str
    status: str
    priority: str
    priority_score: float
    steps_count: int
    involves_tk: bool
    created_at: str
    message: str = ""


class RegisterMemberRequest(BaseModel):
    """Request to register DAO member."""
    wallet_address: str
    display_name: str
    governance_tokens: int = 0
    roles: List[str] = Field(default_factory=list)
    tk_communities: List[str] = Field(default_factory=list)
    is_validator: bool = False


class MemberResponse(BaseModel):
    """DAO member response."""
    member_id: str
    wallet_address: str
    display_name: str
    roles: List[str]
    governance_tokens: int
    is_tk_holder: bool
    voting_power: Dict[str, Any]


class CreateProposalRequest(BaseModel):
    """Request to create proposal."""
    creator_id: str
    title: str
    description: str
    proposal_type: str
    involves_tk: bool = False
    affected_communities: List[str] = Field(default_factory=list)
    affected_pathways: List[str] = Field(default_factory=list)
    execution_data: Dict[str, Any] = Field(default_factory=dict)


class ProposalResponse(BaseModel):
    """Proposal response."""
    proposal_id: str
    title: str
    status: str
    proposal_type: str
    impact_level: str
    involves_tk: bool
    requires_elder_endorsement: bool
    total_for_power: float
    total_against_power: float
    quorum_reached: bool


class CastVoteRequest(BaseModel):
    """Request to cast vote."""
    member_id: str
    support: bool
    tokens_to_use: int
    rationale: Optional[str] = None


class VoteResponse(BaseModel):
    """Vote response."""
    vote_id: str
    proposal_id: str
    support: bool
    base_vote_power: float
    multiplier: float
    final_vote_power: float


class CreateEmergencyRequest(BaseModel):
    """Request to create emergency response."""
    crisis_type: str
    crisis_level: str
    title: str
    description: str
    involves_tk: bool = False
    affected_communities: List[str] = Field(default_factory=list)


class EmergencyResponse(BaseModel):
    """Emergency response."""
    crisis_id: str
    title: str
    crisis_type: str
    crisis_level: str
    status: str
    affected_pathways: int
    involves_tk: bool


class InitiateResponseRequest(BaseModel):
    """Request to initiate emergency response."""
    incident_commander: str
    response_team: List[str] = Field(default_factory=list)


class TKConsultationRequest(BaseModel):
    """Request to record TK consultation."""
    guidance: str
    concerns: List[str] = Field(default_factory=list)
    recommendations: List[str] = Field(default_factory=list)
    required_actions: List[str] = Field(default_factory=list)


class ResolveCrisisRequest(BaseModel):
    """Request to resolve crisis."""
    root_cause: str
    lessons_learned: List[str] = Field(default_factory=list)
    preventive_measures: List[str] = Field(default_factory=list)


# ==============================================================================
# Orchestrator Endpoints
# ==============================================================================

@router.get("/pathways", summary="List all 31 pathways")
async def list_pathways():
    """Get all pathway configurations."""
    pathways = _orchestrator.get_all_pathways()
    return {
        "total_pathways": len(pathways),
        "pathways": [
            {
                "pathway_id": p.pathway_id.value,
                "layer": p.layer.value,
                "baseline_qubits": p.baseline_qubits,
                "baseline_memory_gb": p.baseline_memory_gb,
                "baseline_storage_tb": p.baseline_storage_tb,
                "is_active": p.is_active,
                "health_score": p.health_score,
            }
            for p in pathways
        ],
    }


@router.get("/pathways/{pathway_id}", summary="Get pathway details")
async def get_pathway(
    pathway_id: str = Path(..., description="Pathway ID"),
):
    """Get specific pathway configuration."""
    try:
        pid = PathwayID(pathway_id)
    except ValueError:
        raise HTTPException(status_code=400, detail=f"Invalid pathway ID: {pathway_id}")
    
    config = _orchestrator.get_pathway_status(pid)
    if not config:
        raise HTTPException(status_code=404, detail="Pathway not found")
    
    dependencies = _orchestrator.get_pathway_dependencies(pid)
    dependents = _orchestrator.get_dependent_pathways(pid)
    
    return {
        "pathway": {
            "pathway_id": config.pathway_id.value,
            "layer": config.layer.value,
            "baseline_qubits": config.baseline_qubits,
            "baseline_memory_gb": config.baseline_memory_gb,
            "baseline_storage_tb": config.baseline_storage_tb,
            "is_active": config.is_active,
            "health_score": config.health_score,
            "current_load": config.current_load,
        },
        "dependencies": {
            "required": [d.value for d in dependencies["required"]],
            "optional": [d.value for d in dependencies["optional"]],
        },
        "dependent_pathways": [d.value for d in dependents],
    }


@router.post("/workflows", summary="Create workflow", response_model=WorkflowResponse)
async def create_workflow(request: CreateWorkflowRequest):
    """Create a new multi-pathway workflow."""
    try:
        priority = PriorityLevel(request.priority)
    except ValueError:
        raise HTTPException(status_code=400, detail=f"Invalid priority: {request.priority}")
    
    try:
        steps = [
            {
                "pathway": step.pathway,
                "operation": step.operation,
                "parameters": step.parameters,
                "depends_on": step.depends_on,
                "involves_tk": step.involves_tk,
                "tk_communities": step.tk_communities,
            }
            for step in request.steps
        ]
        
        workflow = _orchestrator.create_workflow(
            name=request.name,
            steps=steps,
            priority=priority,
            description=request.description,
            involves_tk=request.involves_tk,
            created_by=request.created_by,
            dao_proposal_id=request.dao_proposal_id,
        )
        
        return WorkflowResponse(
            workflow_id=workflow.workflow_id,
            name=workflow.name,
            status=workflow.status.value,
            priority=workflow.priority.value,
            priority_score=workflow.priority_score,
            steps_count=len(workflow.steps),
            involves_tk=workflow.involves_tk,
            created_at=workflow.created_at.isoformat(),
            message="Workflow created successfully",
        )
    
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))


@router.get("/workflows/{workflow_id}", summary="Get workflow")
async def get_workflow(workflow_id: str):
    """Get workflow details."""
    workflow = _orchestrator.get_workflow(workflow_id)
    if not workflow:
        raise HTTPException(status_code=404, detail="Workflow not found")
    
    return {
        "workflow_id": workflow.workflow_id,
        "name": workflow.name,
        "description": workflow.description,
        "status": workflow.status.value,
        "priority": workflow.priority.value,
        "priority_score": workflow.priority_score,
        "involves_tk": workflow.involves_tk,
        "created_at": workflow.created_at.isoformat(),
        "started_at": workflow.started_at.isoformat() if workflow.started_at else None,
        "completed_at": workflow.completed_at.isoformat() if workflow.completed_at else None,
        "steps": [
            {
                "step_id": s.step_id,
                "pathway_id": s.pathway_id.value,
                "operation": s.operation,
                "status": s.status.value,
                "depends_on": s.depends_on,
            }
            for s in workflow.steps
        ],
        "resource_allocation": {
            "qubits": workflow.resource_allocation.qubits,
            "memory_gb": workflow.resource_allocation.memory_gb,
            "storage_tb": workflow.resource_allocation.storage_tb,
        } if workflow.resource_allocation else None,
    }


@router.post("/workflows/{workflow_id}/start", summary="Start workflow")
async def start_workflow(workflow_id: str):
    """Start workflow execution."""
    try:
        workflow = _orchestrator.start_workflow(workflow_id)
        return {
            "workflow_id": workflow.workflow_id,
            "status": workflow.status.value,
            "started_at": workflow.started_at.isoformat(),
            "resource_allocation": {
                "qubits": workflow.resource_allocation.qubits,
                "memory_gb": workflow.resource_allocation.memory_gb,
            } if workflow.resource_allocation else None,
            "message": "Workflow started",
        }
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))


@router.get("/workflows/{workflow_id}/execution-order", summary="Get execution order")
async def get_execution_order(workflow_id: str):
    """Get optimized step execution order."""
    try:
        steps = _orchestrator.get_execution_order(workflow_id)
        return {
            "workflow_id": workflow_id,
            "execution_order": [
                {
                    "step_id": s.step_id,
                    "pathway_id": s.pathway_id.value,
                    "operation": s.operation,
                }
                for s in steps
            ],
        }
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))


@router.get("/metrics", summary="Get coordination metrics")
async def get_metrics():
    """Get current coordination metrics."""
    metrics = _orchestrator.get_metrics()
    return {
        "timestamp": metrics.timestamp.isoformat(),
        "performance": {
            "orchestration_efficiency_seconds": metrics.orchestration_efficiency_seconds,
            "coordination_accuracy": metrics.coordination_accuracy,
            "tk_preservation_accuracy": metrics.tk_preservation_accuracy,
            "decentralization_coefficient": metrics.decentralization_coefficient,
        },
        "resources": {
            "total_qubits_allocated": metrics.total_qubits_allocated,
            "total_memory_allocated_gb": metrics.total_memory_allocated_gb,
            "total_storage_allocated_tb": metrics.total_storage_allocated_tb,
        },
        "workflows": {
            "active": metrics.active_workflows,
            "queued": metrics.queued_workflows,
            "completed_today": metrics.completed_today,
            "failed_today": metrics.failed_today,
        },
        "pathway_health": {
            "healthy": metrics.healthy_pathways,
            "degraded": metrics.degraded_pathways,
            "offline": metrics.offline_pathways,
        },
    }


@router.get("/resources", summary="Get resource utilization")
async def get_resource_utilization():
    """Get current resource utilization."""
    return _orchestrator.get_resource_utilization()


# ==============================================================================
# Governance Endpoints
# ==============================================================================

@router.post("/governance/members", summary="Register DAO member", response_model=MemberResponse)
async def register_member(request: RegisterMemberRequest):
    """Register a new DAO member."""
    try:
        roles = [MemberRole(r) for r in request.roles] if request.roles else None
    except ValueError as e:
        raise HTTPException(status_code=400, detail=f"Invalid role: {e}")
    
    member = _governance.register_member(
        wallet_address=request.wallet_address,
        display_name=request.display_name,
        governance_tokens=request.governance_tokens,
        roles=roles,
        tk_communities=request.tk_communities,
        is_validator=request.is_validator,
    )
    
    voting_power = _governance.get_member_voting_power(member.member_id)
    
    return MemberResponse(
        member_id=member.member_id,
        wallet_address=member.wallet_address,
        display_name=member.display_name,
        roles=[r.value for r in member.roles],
        governance_tokens=member.governance_tokens,
        is_tk_holder=member.is_tk_holder,
        voting_power=voting_power,
    )


@router.get("/governance/members/{member_id}", summary="Get member")
async def get_member(member_id: str):
    """Get DAO member details."""
    member = _governance.get_member(member_id)
    if not member:
        raise HTTPException(status_code=404, detail="Member not found")
    
    voting_power = _governance.get_member_voting_power(member_id)
    
    return {
        "member": {
            "member_id": member.member_id,
            "wallet_address": member.wallet_address,
            "display_name": member.display_name,
            "roles": [r.value for r in member.roles],
            "governance_tokens": member.governance_tokens,
            "staked_tokens": member.staked_tokens,
            "is_tk_holder": member.is_tk_holder,
            "tk_communities": member.tk_communities,
            "is_validator": member.is_validator,
            "proposals_created": member.proposals_created,
            "votes_cast": member.votes_cast,
            "joined_at": member.joined_at.isoformat(),
        },
        "voting_power": voting_power,
    }


@router.post("/governance/proposals", summary="Create proposal", response_model=ProposalResponse)
async def create_proposal(request: CreateProposalRequest):
    """Create a governance proposal."""
    try:
        proposal_type = ProposalType(request.proposal_type)
    except ValueError:
        raise HTTPException(status_code=400, detail=f"Invalid proposal type: {request.proposal_type}")
    
    try:
        proposal = _governance.create_proposal(
            creator_id=request.creator_id,
            title=request.title,
            description=request.description,
            proposal_type=proposal_type,
            involves_tk=request.involves_tk,
            affected_communities=request.affected_communities,
            affected_pathways=request.affected_pathways,
            execution_data=request.execution_data,
        )
        
        return ProposalResponse(
            proposal_id=proposal.proposal_id,
            title=proposal.title,
            status=proposal.status.value,
            proposal_type=proposal.proposal_type.value,
            impact_level=proposal.impact_level.value,
            involves_tk=proposal.involves_tk,
            requires_elder_endorsement=proposal.requires_elder_endorsement,
            total_for_power=proposal.total_for_power,
            total_against_power=proposal.total_against_power,
            quorum_reached=proposal.quorum_reached,
        )
    
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))


@router.get("/governance/proposals/{proposal_id}", summary="Get proposal")
async def get_proposal(proposal_id: str):
    """Get proposal details."""
    proposal = _governance.get_proposal(proposal_id)
    if not proposal:
        raise HTTPException(status_code=404, detail="Proposal not found")
    
    return {
        "proposal_id": proposal.proposal_id,
        "title": proposal.title,
        "description": proposal.description,
        "proposal_type": proposal.proposal_type.value,
        "status": proposal.status.value,
        "impact_level": proposal.impact_level.value,
        "involves_tk": proposal.involves_tk,
        "affected_communities": proposal.affected_communities,
        "affected_pathways": proposal.affected_pathways,
        "requires_elder_endorsement": proposal.requires_elder_endorsement,
        "elder_endorsements_count": len(proposal.elder_endorsements),
        "minimum_endorsements_required": proposal.minimum_endorsements_required,
        "voting_starts_at": proposal.voting_starts_at.isoformat() if proposal.voting_starts_at else None,
        "voting_ends_at": proposal.voting_ends_at.isoformat() if proposal.voting_ends_at else None,
        "total_for_power": proposal.total_for_power,
        "total_against_power": proposal.total_against_power,
        "votes_count": len(proposal.votes),
        "quorum_reached": proposal.quorum_reached,
        "created_at": proposal.created_at.isoformat(),
    }


@router.post("/governance/proposals/{proposal_id}/vote", summary="Cast vote", response_model=VoteResponse)
async def cast_vote(proposal_id: str, request: CastVoteRequest):
    """Cast vote on a proposal."""
    try:
        vote = _governance.cast_vote(
            proposal_id=proposal_id,
            member_id=request.member_id,
            support=request.support,
            tokens_to_use=request.tokens_to_use,
            rationale=request.rationale,
        )
        
        return VoteResponse(
            vote_id=vote.vote_id,
            proposal_id=vote.proposal_id,
            support=vote.support,
            base_vote_power=vote.base_vote_power,
            multiplier=vote.multiplier,
            final_vote_power=vote.final_vote_power,
        )
    
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))


@router.post("/governance/proposals/{proposal_id}/open-voting", summary="Open voting")
async def open_voting(proposal_id: str):
    """Open voting on a proposal."""
    try:
        proposal = _governance.open_voting(proposal_id)
        return {
            "proposal_id": proposal.proposal_id,
            "status": proposal.status.value,
            "voting_starts_at": proposal.voting_starts_at.isoformat(),
            "voting_ends_at": proposal.voting_ends_at.isoformat(),
        }
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))


@router.post("/governance/proposals/{proposal_id}/finalize", summary="Finalize proposal")
async def finalize_proposal(proposal_id: str):
    """Finalize proposal voting."""
    try:
        proposal = _governance.finalize_proposal(proposal_id)
        return {
            "proposal_id": proposal.proposal_id,
            "status": proposal.status.value,
            "total_for_power": proposal.total_for_power,
            "total_against_power": proposal.total_against_power,
            "quorum_reached": proposal.quorum_reached,
        }
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))


@router.get("/governance/stats", summary="Get governance statistics")
async def get_governance_stats():
    """Get overall governance statistics."""
    return _governance.get_governance_stats()


# ==============================================================================
# Emergency Endpoints
# ==============================================================================

@router.post("/emergency/responses", summary="Create emergency response", response_model=EmergencyResponse)
async def create_emergency_response(request: CreateEmergencyRequest):
    """Create a new emergency response."""
    try:
        crisis_type = CrisisType(request.crisis_type)
        crisis_level = CrisisLevel(request.crisis_level)
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
    
    response = _emergency.create_emergency_response(
        crisis_type=crisis_type,
        crisis_level=crisis_level,
        title=request.title,
        description=request.description,
        involves_tk=request.involves_tk,
        affected_communities=request.affected_communities,
    )
    
    return EmergencyResponse(
        crisis_id=response.crisis_id,
        title=response.title,
        crisis_type=response.crisis_type.value,
        crisis_level=response.crisis_level.value,
        status=response.status.value,
        affected_pathways=len(response.affected_pathways),
        involves_tk=response.involves_tk,
    )


@router.get("/emergency/responses/{crisis_id}", summary="Get emergency response")
async def get_emergency_response(crisis_id: str):
    """Get emergency response details."""
    response = _emergency.get_response(crisis_id)
    if not response:
        raise HTTPException(status_code=404, detail="Crisis not found")
    
    return {
        "crisis_id": response.crisis_id,
        "title": response.title,
        "description": response.description,
        "crisis_type": response.crisis_type.value,
        "crisis_level": response.crisis_level.value,
        "status": response.status.value,
        "detection_time": response.detection_time.isoformat(),
        "response_initiated_at": response.response_initiated_at.isoformat() if response.response_initiated_at else None,
        "contained_at": response.contained_at.isoformat() if response.contained_at else None,
        "resolved_at": response.resolved_at.isoformat() if response.resolved_at else None,
        "incident_commander": response.incident_commander,
        "response_team": response.response_team,
        "affected_pathways": [
            {
                "pathway_name": p.pathway_name,
                "mode": p.mode.value,
                "health_score": p.health_score,
            }
            for p in response.affected_pathways
        ],
        "involves_tk": response.involves_tk,
        "tk_consultations_count": len(response.tk_consultations),
        "actions_count": len(response.actions_log),
    }


@router.post("/emergency/responses/{crisis_id}/initiate", summary="Initiate response")
async def initiate_response(crisis_id: str, request: InitiateResponseRequest):
    """Initiate emergency response."""
    try:
        response = _emergency.initiate_response(
            crisis_id=crisis_id,
            incident_commander=request.incident_commander,
            response_team=request.response_team,
        )
        return {
            "crisis_id": response.crisis_id,
            "status": response.status.value,
            "response_initiated_at": response.response_initiated_at.isoformat(),
            "message": "Response initiated",
        }
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))


@router.post("/emergency/responses/{crisis_id}/allocate-resources", summary="Allocate emergency resources")
async def allocate_emergency_resources(
    crisis_id: str,
    qubits: int = Query(1000, ge=0),
    memory_gb: int = Query(100, ge=0),
):
    """Allocate emergency resources to crisis."""
    try:
        response = _emergency.allocate_emergency_resources(crisis_id, qubits, memory_gb)
        return {
            "crisis_id": response.crisis_id,
            "emergency_resources_used": response.emergency_resources_used,
            "message": "Resources allocated",
        }
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))


@router.post("/emergency/responses/{crisis_id}/contain", summary="Contain crisis")
async def contain_crisis(crisis_id: str):
    """Mark crisis as contained."""
    try:
        response = _emergency.contain_crisis(crisis_id)
        return {
            "crisis_id": response.crisis_id,
            "status": response.status.value,
            "contained_at": response.contained_at.isoformat(),
        }
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))


@router.post("/emergency/responses/{crisis_id}/resolve", summary="Resolve crisis")
async def resolve_crisis(crisis_id: str, request: ResolveCrisisRequest):
    """Resolve crisis and complete post-analysis."""
    try:
        response = _emergency.resolve_crisis(
            crisis_id=crisis_id,
            root_cause=request.root_cause,
            lessons_learned=request.lessons_learned,
            preventive_measures=request.preventive_measures,
        )
        return {
            "crisis_id": response.crisis_id,
            "status": response.status.value,
            "resolved_at": response.resolved_at.isoformat(),
            "root_cause": response.root_cause,
        }
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))


@router.get("/emergency/responses/{crisis_id}/metrics", summary="Get response metrics")
async def get_response_metrics(crisis_id: str):
    """Get metrics for a crisis response."""
    metrics = _emergency.get_response_metrics(crisis_id)
    if not metrics:
        raise HTTPException(status_code=404, detail="Crisis not found")
    return metrics


@router.get("/emergency/active", summary="Get active emergencies")
async def get_active_emergencies():
    """Get all active emergency responses."""
    responses = _emergency.get_active_responses()
    return {
        "count": len(responses),
        "responses": [
            {
                "crisis_id": r.crisis_id,
                "title": r.title,
                "crisis_level": r.crisis_level.value,
                "status": r.status.value,
            }
            for r in responses
        ],
    }


# ==============================================================================
# Dashboard Endpoints
# ==============================================================================

@router.get("/dashboard/overview", summary="Get system overview")
async def get_dashboard_overview():
    """Get system overview dashboard."""
    return _dashboard.get_system_overview()


@router.get("/dashboard/pathways", summary="Get pathway dashboard")
async def get_pathway_dashboard(
    pathway_id: Optional[str] = Query(None, description="Specific pathway ID"),
):
    """Get pathway performance dashboard."""
    return _dashboard.get_pathway_dashboard(pathway_id)


@router.get("/dashboard/tk-preservation", summary="Get TK preservation dashboard")
async def get_tk_preservation_dashboard(
    privacy_level: str = Query("public", description="Privacy level"),
):
    """Get TK preservation metrics (privacy-filtered)."""
    try:
        level = PrivacyLevel(privacy_level)
    except ValueError:
        raise HTTPException(status_code=400, detail=f"Invalid privacy level: {privacy_level}")
    
    metrics = _dashboard.get_tk_preservation_dashboard(level)
    return metrics.model_dump()


@router.get("/dashboard/governance", summary="Get governance dashboard")
async def get_governance_dashboard():
    """Get governance metrics dashboard."""
    return _dashboard.get_governance_dashboard().model_dump()


@router.get("/dashboard/resources", summary="Get resource dashboard")
async def get_resource_dashboard():
    """Get resource utilization dashboard."""
    return _dashboard.get_resource_dashboard().model_dump()


@router.get("/dashboard/alerts", summary="Get active alerts")
async def get_dashboard_alerts(
    severity: Optional[str] = Query(None, description="Filter by severity"),
    category: Optional[str] = Query(None, description="Filter by category"),
):
    """Get active dashboard alerts."""
    sev = AlertSeverity(severity) if severity else None
    cat = MetricCategory(category) if category else None
    
    alerts = _dashboard.get_active_alerts(sev, cat)
    return {
        "count": len(alerts),
        "alerts": [a.model_dump() for a in alerts],
    }


@router.post("/dashboard/alerts/{alert_id}/acknowledge", summary="Acknowledge alert")
async def acknowledge_alert(
    alert_id: str,
    acknowledged_by: str = Query(..., description="User acknowledging"),
):
    """Acknowledge a dashboard alert."""
    try:
        alert = _dashboard.acknowledge_alert(alert_id, acknowledged_by)
        return {
            "alert_id": alert.alert_id,
            "acknowledged": alert.acknowledged,
            "acknowledged_by": alert.acknowledged_by,
        }
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))


@router.get("/dashboard/export", summary="Export dashboard snapshot")
async def export_dashboard(
    privacy_level: str = Query("public", description="Privacy level"),
):
    """Export complete dashboard snapshot."""
    try:
        level = PrivacyLevel(privacy_level)
    except ValueError:
        raise HTTPException(status_code=400, detail=f"Invalid privacy level: {privacy_level}")
    
    return _dashboard.export_dashboard_snapshot(level)
