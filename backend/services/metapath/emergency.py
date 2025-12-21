"""
MetaPath Emergency Coordination Module

Trade Secret: Crisis detection algorithms, emergency pathway activation,
0.8-hour response time protocols, TK consultation integration during emergencies.

Key Features:
- 0.8-hour crisis activation target
- Automatic pathway rerouting
- TK community emergency notification
- Resource surge allocation
- Post-crisis analysis and reporting
"""

from enum import Enum
from typing import Dict, List, Optional, Set, Any, Tuple
from pydantic import BaseModel, Field
from datetime import datetime, timedelta
import uuid


class CrisisLevel(str, Enum):
    """Crisis severity levels."""
    ADVISORY = "advisory"          # Monitoring, no action
    ELEVATED = "elevated"          # Increased monitoring
    WARNING = "warning"            # Prepare for action
    EMERGENCY = "emergency"        # Immediate action required
    CRITICAL = "critical"          # Full system response


class CrisisType(str, Enum):
    """Types of crisis scenarios."""
    # Technical Crises
    SYSTEM_FAILURE = "system_failure"
    SECURITY_BREACH = "security_breach"
    DATA_CORRUPTION = "data_corruption"
    NETWORK_OUTAGE = "network_outage"
    
    # TK-Related Crises
    TK_MISAPPROPRIATION = "tk_misappropriation"
    COMMUNITY_DISPUTE = "community_dispute"
    ATTRIBUTION_FAILURE = "attribution_failure"
    
    # External Crises
    REGULATORY_ACTION = "regulatory_action"
    LEGAL_CHALLENGE = "legal_challenge"
    MARKET_DISRUPTION = "market_disruption"
    
    # Environmental Crises
    CLIMATE_EVENT = "climate_event"
    SUPPLY_CHAIN_DISRUPTION = "supply_chain_disruption"


class ResponseStatus(str, Enum):
    """Emergency response status."""
    DETECTED = "detected"
    ASSESSING = "assessing"
    RESPONDING = "responding"
    CONTAINED = "contained"
    RECOVERING = "recovering"
    RESOLVED = "resolved"
    POST_ANALYSIS = "post_analysis"


class PathwayOverrideMode(str, Enum):
    """Pathway operation modes during emergency."""
    NORMAL = "normal"
    DEGRADED = "degraded"
    EMERGENCY = "emergency"
    SUSPENDED = "suspended"
    REROUTED = "rerouted"


class EmergencyContact(BaseModel):
    """Emergency contact for notifications."""
    contact_id: str = Field(default_factory=lambda: str(uuid.uuid4()))
    name: str
    role: str
    
    # Contact methods
    email: str = ""
    phone: str = ""
    signal_id: str = ""  # Secure messaging
    
    # Notification preferences
    crisis_levels: List[CrisisLevel] = Field(default_factory=lambda: [CrisisLevel.EMERGENCY, CrisisLevel.CRITICAL])
    crisis_types: List[CrisisType] = Field(default_factory=list)
    
    # Community association
    community_name: Optional[str] = None
    is_tk_community_contact: bool = False
    is_elder: bool = False


class CrisisIndicator(BaseModel):
    """Indicator that triggered crisis detection."""
    indicator_id: str = Field(default_factory=lambda: str(uuid.uuid4()))
    source: str  # System, pathway, external
    metric_name: str
    
    # Values
    threshold_value: float
    actual_value: float
    deviation_percentage: float
    
    # Timing
    detected_at: datetime = Field(default_factory=datetime.utcnow)
    duration_seconds: int = 0


class PathwayStatus(BaseModel):
    """Status of pathway during emergency."""
    pathway_id: str
    pathway_name: str
    mode: PathwayOverrideMode = PathwayOverrideMode.NORMAL
    
    # Status
    is_affected: bool = False
    health_score: float = 1.0
    load_factor: float = 0.0
    
    # Actions
    rerouted_to: Optional[str] = None
    suspended_operations: List[str] = Field(default_factory=list)
    
    # Resources (quantum-ready: compute_units scale to qubits)
    emergency_resources_allocated: bool = False
    additional_compute_units: int = 0
    additional_memory_gb: int = 0


class TKConsultation(BaseModel):
    """TK community consultation during emergency."""
    consultation_id: str = Field(default_factory=lambda: str(uuid.uuid4()))
    crisis_id: str
    community_name: str
    
    # Consultation
    contacted_at: datetime = Field(default_factory=datetime.utcnow)
    responded_at: Optional[datetime] = None
    
    # Response
    consultation_complete: bool = False
    guidance_provided: str = ""
    concerns_raised: List[str] = Field(default_factory=list)
    recommendations: List[str] = Field(default_factory=list)
    
    # Impact on response
    impacts_response: bool = False
    required_actions: List[str] = Field(default_factory=list)


class EmergencyResponse(BaseModel):
    """Active emergency response."""
    crisis_id: str = Field(default_factory=lambda: str(uuid.uuid4()))
    
    # Crisis identification
    crisis_type: CrisisType
    crisis_level: CrisisLevel
    title: str
    description: str
    
    # Status tracking
    status: ResponseStatus = ResponseStatus.DETECTED
    created_at: datetime = Field(default_factory=datetime.utcnow)
    detection_time: datetime = Field(default_factory=datetime.utcnow)
    response_initiated_at: Optional[datetime] = None
    contained_at: Optional[datetime] = None
    resolved_at: Optional[datetime] = None
    
    # Indicators
    indicators: List[CrisisIndicator] = Field(default_factory=list)
    
    # Affected systems
    affected_pathways: List[PathwayStatus] = Field(default_factory=list)
    
    # TK involvement
    involves_tk: bool = False
    affected_communities: List[str] = Field(default_factory=list)
    tk_consultations: List[TKConsultation] = Field(default_factory=list)
    
    # Response team
    incident_commander: Optional[str] = None
    response_team: List[str] = Field(default_factory=list)
    
    # Actions taken
    actions_log: List[Dict[str, Any]] = Field(default_factory=list)
    
    # Resources
    emergency_resources_used: Dict[str, int] = Field(default_factory=dict)
    
    # Impact
    estimated_impact_score: float = 0.0
    actual_impact_score: float = 0.0
    
    # Post-crisis
    root_cause: Optional[str] = None
    lessons_learned: List[str] = Field(default_factory=list)
    preventive_measures: List[str] = Field(default_factory=list)


class EmergencyConfig(BaseModel):
    """Emergency coordination configuration."""
    # Response time targets
    target_detection_seconds: int = 60
    target_response_activation_seconds: int = 300  # 5 minutes
    target_crisis_activation_hours: float = 0.8
    
    # Escalation thresholds
    auto_escalation_enabled: bool = True
    escalation_timeout_minutes: int = 15
    
    # Resource surge limits (quantum-ready: compute_units map to qubits)
    max_emergency_compute_units: int = 5000
    max_emergency_memory_gb: int = 500
    
    # TK consultation requirements
    tk_consultation_required_types: List[CrisisType] = Field(default_factory=lambda: [
        CrisisType.TK_MISAPPROPRIATION,
        CrisisType.COMMUNITY_DISPUTE,
        CrisisType.ATTRIBUTION_FAILURE,
    ])


class EmergencyCoordinator:
    """
    Emergency Coordination System for MetaPath.
    
    Trade Secret: Crisis detection algorithms, response time optimization,
    TK consultation protocols, resource surge allocation.
    
    Performance Targets:
    - Crisis activation: 0.8 hours
    - Detection latency: < 60 seconds
    - Response initiation: < 5 minutes
    """
    
    # ==========================================================================
    # TRADE SECRET: Crisis Detection Thresholds
    # ==========================================================================
    _DETECTION_THRESHOLDS = {
        # System metrics
        "pathway_health": {"warning": 0.7, "emergency": 0.5, "critical": 0.3},
        "error_rate": {"warning": 0.05, "emergency": 0.1, "critical": 0.2},
        "latency_ms": {"warning": 500, "emergency": 1000, "critical": 2000},
        "resource_utilization": {"warning": 0.8, "emergency": 0.9, "critical": 0.95},
        
        # TK metrics
        "attribution_accuracy": {"warning": 0.9, "emergency": 0.8, "critical": 0.7},
        "tk_preservation_score": {"warning": 0.9, "emergency": 0.85, "critical": 0.8},
        
        # Security metrics
        "auth_failure_rate": {"warning": 0.02, "emergency": 0.05, "critical": 0.1},
        "anomaly_score": {"warning": 0.6, "emergency": 0.75, "critical": 0.9},
    }
    
    # ==========================================================================
    # TRADE SECRET: Pathway Rerouting Matrix
    # ==========================================================================
    _REROUTING_MATRIX = {
        "TechPath": ["GenomePath", "ChemPath"],
        "GenomePath": ["ChemPath", "BioPath"],
        "ChemPath": ["GenomePath", "BioPath"],
        "BioPath": ["ChemPath", "ToxPath"],
        "ToxPath": ["BioPath", "ChemPath"],
        "EthnoPath": ["EduPath", "EquiPath"],
        "ClimatePath": ["EcoPath", "GeoPath"],
        "RetroPath": ["TechPath", "EthnoPath"],
        "ClinPath": ["BioPath", "ToxPath"],
        "RegPath": ["PatentPath", "ClinPath"],
        "EquiPath": ["FundPath", "J.E.D.IPath"],
    }
    
    # ==========================================================================
    # TRADE SECRET: Crisis Type to Pathway Impact Mapping
    # ==========================================================================
    _CRISIS_PATHWAY_IMPACT = {
        CrisisType.SYSTEM_FAILURE: ["TechPath", "MetaPath", "GenomePath"],
        CrisisType.SECURITY_BREACH: ["TechPath", "MetaPath", "EthnoPath", "EquiPath"],
        CrisisType.DATA_CORRUPTION: ["GenomePath", "ChemPath", "BioPath", "TechPath"],
        CrisisType.NETWORK_OUTAGE: ["TechPath", "MetaPath", "LogistiPath"],
        CrisisType.TK_MISAPPROPRIATION: ["EthnoPath", "EquiPath", "J.E.D.IPath", "PatentPath"],
        CrisisType.COMMUNITY_DISPUTE: ["EquiPath", "J.E.D.IPath", "PopuliPath", "EthnoPath"],
        CrisisType.ATTRIBUTION_FAILURE: ["EquiPath", "EthnoPath", "PatentPath"],
        CrisisType.REGULATORY_ACTION: ["RegPath", "PatentPath", "PharmaPath", "MarketPath"],
        CrisisType.LEGAL_CHALLENGE: ["PatentPath", "RegPath", "J.E.D.IPath"],
        CrisisType.MARKET_DISRUPTION: ["MarketPath", "PharmaPath", "FundPath"],
        CrisisType.CLIMATE_EVENT: ["ClimatePath", "EcoPath", "ArgoPath", "GeoPath"],
        CrisisType.SUPPLY_CHAIN_DISRUPTION: ["LogistiPath", "PharmaPath", "ArgoPath"],
    }
    
    def __init__(self, config: Optional[EmergencyConfig] = None):
        self._config = config or EmergencyConfig()
        self._contacts: Dict[str, EmergencyContact] = {}
        self._active_responses: Dict[str, EmergencyResponse] = {}
        self._historical_responses: Dict[str, EmergencyResponse] = {}
    
    # ==========================================================================
    # Crisis Detection
    # ==========================================================================
    def detect_crisis(
        self,
        metrics: Dict[str, float],
    ) -> Optional[Tuple[CrisisLevel, List[CrisisIndicator]]]:
        """
        Analyze metrics to detect potential crisis.
        Trade Secret: Multi-metric crisis detection algorithm.
        """
        indicators = []
        max_level = CrisisLevel.ADVISORY
        
        for metric_name, value in metrics.items():
            if metric_name not in self._DETECTION_THRESHOLDS:
                continue
            
            thresholds = self._DETECTION_THRESHOLDS[metric_name]
            
            # Check threshold breaches
            for level_name, threshold in thresholds.items():
                # Handle both high-is-bad and low-is-bad metrics
                is_breach = False
                if metric_name in ["pathway_health", "attribution_accuracy", "tk_preservation_score"]:
                    is_breach = value < threshold
                else:
                    is_breach = value > threshold
                
                if is_breach:
                    indicator = CrisisIndicator(
                        source="system",
                        metric_name=metric_name,
                        threshold_value=threshold,
                        actual_value=value,
                        deviation_percentage=abs(value - threshold) / threshold * 100,
                    )
                    indicators.append(indicator)
                    
                    # Update max level
                    level_map = {
                        "warning": CrisisLevel.WARNING,
                        "emergency": CrisisLevel.EMERGENCY,
                        "critical": CrisisLevel.CRITICAL,
                    }
                    current_level = level_map.get(level_name, CrisisLevel.ADVISORY)
                    if list(CrisisLevel).index(current_level) > list(CrisisLevel).index(max_level):
                        max_level = current_level
        
        if indicators:
            return max_level, indicators
        return None
    
    # ==========================================================================
    # Emergency Response Creation
    # ==========================================================================
    def create_emergency_response(
        self,
        crisis_type: CrisisType,
        crisis_level: CrisisLevel,
        title: str,
        description: str,
        indicators: Optional[List[CrisisIndicator]] = None,
        involves_tk: bool = False,
        affected_communities: Optional[List[str]] = None,
    ) -> EmergencyResponse:
        """
        Create new emergency response.
        """
        # Determine affected pathways
        affected_pathway_names = self._CRISIS_PATHWAY_IMPACT.get(crisis_type, [])
        pathway_statuses = []
        
        for pathway_name in affected_pathway_names:
            status = PathwayStatus(
                pathway_id=pathway_name.lower().replace(" ", "_"),
                pathway_name=pathway_name,
                mode=PathwayOverrideMode.EMERGENCY if crisis_level in [CrisisLevel.EMERGENCY, CrisisLevel.CRITICAL] else PathwayOverrideMode.DEGRADED,
                is_affected=True,
            )
            pathway_statuses.append(status)
        
        # Check if TK consultation required
        if crisis_type in self._config.tk_consultation_required_types:
            involves_tk = True
        
        response = EmergencyResponse(
            crisis_type=crisis_type,
            crisis_level=crisis_level,
            title=title,
            description=description,
            indicators=indicators or [],
            affected_pathways=pathway_statuses,
            involves_tk=involves_tk,
            affected_communities=affected_communities or [],
        )
        
        self._active_responses[response.crisis_id] = response
        
        # Log creation action
        response.actions_log.append({
            "timestamp": datetime.utcnow().isoformat(),
            "action": "response_created",
            "details": f"Emergency response created for {crisis_type.value}",
        })
        
        return response
    
    # ==========================================================================
    # Response Management
    # ==========================================================================
    def initiate_response(
        self,
        crisis_id: str,
        incident_commander: str,
        response_team: Optional[List[str]] = None,
    ) -> EmergencyResponse:
        """
        Initiate emergency response.
        Trade Secret: Response activation protocol.
        """
        response = self._active_responses.get(crisis_id)
        if not response:
            raise ValueError(f"Crisis not found: {crisis_id}")
        
        response.status = ResponseStatus.RESPONDING
        response.response_initiated_at = datetime.utcnow()
        response.incident_commander = incident_commander
        response.response_team = response_team or []
        
        # Log action
        response.actions_log.append({
            "timestamp": datetime.utcnow().isoformat(),
            "action": "response_initiated",
            "details": f"Response initiated by {incident_commander}",
        })
        
        # Calculate response time
        detection_to_response = (
            response.response_initiated_at - response.detection_time
        ).total_seconds()
        
        # Check against target
        target_seconds = self._config.target_response_activation_seconds
        if detection_to_response > target_seconds:
            response.actions_log.append({
                "timestamp": datetime.utcnow().isoformat(),
                "action": "sla_breach",
                "details": f"Response activation exceeded target ({detection_to_response:.1f}s > {target_seconds}s)",
            })
        
        return response
    
    def allocate_emergency_resources(
        self,
        crisis_id: str,
        compute_units: int,
        memory_gb: int,
    ) -> EmergencyResponse:
        """
        Allocate emergency resources to crisis response.
        Trade Secret: Resource surge allocation algorithm.
        
        Quantum-Ready: compute_units map to qubits when quantum backend available.
        """
        response = self._active_responses.get(crisis_id)
        if not response:
            raise ValueError(f"Crisis not found: {crisis_id}")
        
        # Enforce limits
        compute_units = min(compute_units, self._config.max_emergency_compute_units)
        memory_gb = min(memory_gb, self._config.max_emergency_memory_gb)
        
        response.emergency_resources_used = {
            "compute_units": compute_units,
            "memory_gb": memory_gb,
        }
        
        # Distribute to affected pathways
        if response.affected_pathways:
            per_pathway_compute_units = compute_units // len(response.affected_pathways)
            per_pathway_memory = memory_gb // len(response.affected_pathways)
            
            for pathway_status in response.affected_pathways:
                pathway_status.emergency_resources_allocated = True
                pathway_status.additional_compute_units = per_pathway_compute_units
                pathway_status.additional_memory_gb = per_pathway_memory
        
        # Log action
        response.actions_log.append({
            "timestamp": datetime.utcnow().isoformat(),
            "action": "resources_allocated",
            "details": f"Allocated {compute_units} compute_units, {memory_gb}GB memory",
        })
        
        return response
    
    def reroute_pathway(
        self,
        crisis_id: str,
        pathway_name: str,
        target_pathway: Optional[str] = None,
    ) -> EmergencyResponse:
        """
        Reroute pathway operations during emergency.
        Trade Secret: Pathway rerouting algorithm.
        """
        response = self._active_responses.get(crisis_id)
        if not response:
            raise ValueError(f"Crisis not found: {crisis_id}")
        
        # Find pathway status
        pathway_status = next(
            (p for p in response.affected_pathways if p.pathway_name == pathway_name),
            None
        )
        
        if not pathway_status:
            raise ValueError(f"Pathway not affected: {pathway_name}")
        
        # Determine reroute target
        if not target_pathway:
            alternatives = self._REROUTING_MATRIX.get(pathway_name, [])
            target_pathway = alternatives[0] if alternatives else None
        
        if target_pathway:
            pathway_status.mode = PathwayOverrideMode.REROUTED
            pathway_status.rerouted_to = target_pathway
            
            response.actions_log.append({
                "timestamp": datetime.utcnow().isoformat(),
                "action": "pathway_rerouted",
                "details": f"Rerouted {pathway_name} to {target_pathway}",
            })
        
        return response
    
    # ==========================================================================
    # TK Consultation
    # ==========================================================================
    def initiate_tk_consultation(
        self,
        crisis_id: str,
        community_name: str,
    ) -> TKConsultation:
        """
        Initiate TK community consultation.
        Trade Secret: TK consultation integration protocol.
        """
        response = self._active_responses.get(crisis_id)
        if not response:
            raise ValueError(f"Crisis not found: {crisis_id}")
        
        consultation = TKConsultation(
            crisis_id=crisis_id,
            community_name=community_name,
        )
        
        response.tk_consultations.append(consultation)
        
        response.actions_log.append({
            "timestamp": datetime.utcnow().isoformat(),
            "action": "tk_consultation_initiated",
            "details": f"TK consultation initiated with {community_name}",
        })
        
        return consultation
    
    def record_tk_consultation_response(
        self,
        crisis_id: str,
        consultation_id: str,
        guidance: str,
        concerns: Optional[List[str]] = None,
        recommendations: Optional[List[str]] = None,
        required_actions: Optional[List[str]] = None,
    ) -> TKConsultation:
        """
        Record TK community consultation response.
        """
        response = self._active_responses.get(crisis_id)
        if not response:
            raise ValueError(f"Crisis not found: {crisis_id}")
        
        consultation = next(
            (c for c in response.tk_consultations if c.consultation_id == consultation_id),
            None
        )
        
        if not consultation:
            raise ValueError(f"Consultation not found: {consultation_id}")
        
        consultation.responded_at = datetime.utcnow()
        consultation.consultation_complete = True
        consultation.guidance_provided = guidance
        consultation.concerns_raised = concerns or []
        consultation.recommendations = recommendations or []
        consultation.required_actions = required_actions or []
        consultation.impacts_response = bool(required_actions)
        
        response.actions_log.append({
            "timestamp": datetime.utcnow().isoformat(),
            "action": "tk_consultation_completed",
            "details": f"TK consultation response from {consultation.community_name}",
        })
        
        return consultation
    
    # ==========================================================================
    # Crisis Resolution
    # ==========================================================================
    def contain_crisis(
        self,
        crisis_id: str,
    ) -> EmergencyResponse:
        """Mark crisis as contained."""
        response = self._active_responses.get(crisis_id)
        if not response:
            raise ValueError(f"Crisis not found: {crisis_id}")
        
        response.status = ResponseStatus.CONTAINED
        response.contained_at = datetime.utcnow()
        
        response.actions_log.append({
            "timestamp": datetime.utcnow().isoformat(),
            "action": "crisis_contained",
            "details": "Crisis successfully contained",
        })
        
        return response
    
    def resolve_crisis(
        self,
        crisis_id: str,
        root_cause: str,
        lessons_learned: Optional[List[str]] = None,
        preventive_measures: Optional[List[str]] = None,
    ) -> EmergencyResponse:
        """
        Resolve crisis and move to post-analysis.
        """
        response = self._active_responses.get(crisis_id)
        if not response:
            raise ValueError(f"Crisis not found: {crisis_id}")
        
        response.status = ResponseStatus.RESOLVED
        response.resolved_at = datetime.utcnow()
        response.root_cause = root_cause
        response.lessons_learned = lessons_learned or []
        response.preventive_measures = preventive_measures or []
        
        # Calculate total response time
        if response.response_initiated_at:
            total_response_hours = (
                response.resolved_at - response.detection_time
            ).total_seconds() / 3600
            
            response.actions_log.append({
                "timestamp": datetime.utcnow().isoformat(),
                "action": "crisis_resolved",
                "details": f"Crisis resolved in {total_response_hours:.2f} hours",
            })
            
            # Check against 0.8 hour target
            if total_response_hours <= self._config.target_crisis_activation_hours:
                response.actions_log.append({
                    "timestamp": datetime.utcnow().isoformat(),
                    "action": "sla_met",
                    "details": f"Response time ({total_response_hours:.2f}h) met 0.8h target",
                })
        
        # Move to historical
        self._historical_responses[crisis_id] = response
        del self._active_responses[crisis_id]
        
        return response
    
    # ==========================================================================
    # Contact Management
    # ==========================================================================
    def register_emergency_contact(
        self,
        name: str,
        role: str,
        email: str = "",
        phone: str = "",
        crisis_levels: Optional[List[CrisisLevel]] = None,
        crisis_types: Optional[List[CrisisType]] = None,
        community_name: Optional[str] = None,
        is_tk_community_contact: bool = False,
        is_elder: bool = False,
    ) -> EmergencyContact:
        """Register emergency contact."""
        contact = EmergencyContact(
            name=name,
            role=role,
            email=email,
            phone=phone,
            crisis_levels=crisis_levels or [CrisisLevel.EMERGENCY, CrisisLevel.CRITICAL],
            crisis_types=crisis_types or [],
            community_name=community_name,
            is_tk_community_contact=is_tk_community_contact,
            is_elder=is_elder,
        )
        
        self._contacts[contact.contact_id] = contact
        return contact
    
    def get_contacts_for_crisis(
        self,
        crisis_type: CrisisType,
        crisis_level: CrisisLevel,
    ) -> List[EmergencyContact]:
        """Get relevant contacts for a crisis."""
        relevant = []
        
        for contact in self._contacts.values():
            # Check level
            if crisis_level not in contact.crisis_levels:
                continue
            
            # Check type (if specified)
            if contact.crisis_types and crisis_type not in contact.crisis_types:
                continue
            
            relevant.append(contact)
        
        return relevant
    
    # ==========================================================================
    # Reporting
    # ==========================================================================
    def get_active_responses(self) -> List[EmergencyResponse]:
        """Get all active emergency responses."""
        return list(self._active_responses.values())
    
    def get_response(self, crisis_id: str) -> Optional[EmergencyResponse]:
        """Get emergency response by ID."""
        return (
            self._active_responses.get(crisis_id) or
            self._historical_responses.get(crisis_id)
        )
    
    def get_response_metrics(
        self,
        crisis_id: str,
    ) -> Dict[str, Any]:
        """Get metrics for a crisis response."""
        response = self.get_response(crisis_id)
        if not response:
            return {}
        
        metrics = {
            "crisis_id": crisis_id,
            "crisis_type": response.crisis_type.value,
            "crisis_level": response.crisis_level.value,
            "status": response.status.value,
        }
        
        # Calculate timing metrics
        if response.response_initiated_at:
            detection_to_response_s = (
                response.response_initiated_at - response.detection_time
            ).total_seconds()
            metrics["detection_to_response_seconds"] = detection_to_response_s
            metrics["response_target_met"] = detection_to_response_s <= self._config.target_response_activation_seconds
        
        if response.resolved_at:
            total_hours = (
                response.resolved_at - response.detection_time
            ).total_seconds() / 3600
            metrics["total_response_hours"] = total_hours
            metrics["0_8_hour_target_met"] = total_hours <= 0.8
        
        # Resource metrics
        metrics["emergency_resources_used"] = response.emergency_resources_used
        metrics["affected_pathways_count"] = len(response.affected_pathways)
        
        # TK metrics
        metrics["involves_tk"] = response.involves_tk
        metrics["tk_consultations_count"] = len(response.tk_consultations)
        metrics["tk_consultations_completed"] = sum(
            1 for c in response.tk_consultations if c.consultation_complete
        )
        
        return metrics
    
    def get_historical_summary(
        self,
        days: int = 30,
    ) -> Dict[str, Any]:
        """Get summary of historical responses."""
        cutoff = datetime.utcnow() - timedelta(days=days)
        
        recent = [
            r for r in self._historical_responses.values()
            if r.resolved_at and r.resolved_at >= cutoff
        ]
        
        return {
            "period_days": days,
            "total_crises": len(recent),
            "by_type": {
                crisis_type.value: sum(1 for r in recent if r.crisis_type == crisis_type)
                for crisis_type in CrisisType
            },
            "by_level": {
                level.value: sum(1 for r in recent if r.crisis_level == level)
                for level in CrisisLevel
            },
            "avg_response_hours": (
                sum(
                    (r.resolved_at - r.detection_time).total_seconds() / 3600
                    for r in recent if r.resolved_at
                ) / len(recent) if recent else 0
            ),
            "target_met_percentage": (
                sum(
                    1 for r in recent
                    if r.resolved_at and 
                    (r.resolved_at - r.detection_time).total_seconds() / 3600 <= 0.8
                ) / len(recent) * 100 if recent else 0
            ),
            "tk_related_count": sum(1 for r in recent if r.involves_tk),
        }
