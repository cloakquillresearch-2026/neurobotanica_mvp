"""
Health Monitoring Configuration

Configuration for monitoring the NeuroBotanica API health,
ML model status, and trade secret engine availability.

Used for:
- Nevada pilot monitoring
- Production alerting
- SLA compliance
"""

from dataclasses import dataclass
from typing import Dict, List, Optional
from enum import Enum


class AlertSeverity(str, Enum):
    """Alert severity levels."""
    CRITICAL = "critical"  # Immediate action required
    WARNING = "warning"  # Attention needed soon
    INFO = "info"  # Informational only


class ServiceStatus(str, Enum):
    """Service status levels."""
    HEALTHY = "healthy"
    DEGRADED = "degraded"
    DOWN = "down"
    UNKNOWN = "unknown"


@dataclass
class HealthThreshold:
    """Threshold configuration for health checks."""
    metric_name: str
    warning_threshold: float
    critical_threshold: float
    comparison: str  # "gt" (greater than) or "lt" (less than)
    unit: str


@dataclass
class AlertConfig:
    """Alert configuration for a specific check."""
    name: str
    description: str
    severity: AlertSeverity
    check_interval_seconds: int
    notification_channels: List[str]


# Health check thresholds
HEALTH_THRESHOLDS: Dict[str, HealthThreshold] = {
    "database_latency_ms": HealthThreshold(
        metric_name="database_latency_ms",
        warning_threshold=100,
        critical_threshold=500,
        comparison="gt",
        unit="milliseconds"
    ),
    "api_response_time_ms": HealthThreshold(
        metric_name="api_response_time_ms",
        warning_threshold=200,
        critical_threshold=1000,
        comparison="gt",
        unit="milliseconds"
    ),
    "ml_model_load_time_ms": HealthThreshold(
        metric_name="ml_model_load_time_ms",
        warning_threshold=5000,
        critical_threshold=30000,
        comparison="gt",
        unit="milliseconds"
    ),
    "memory_usage_percent": HealthThreshold(
        metric_name="memory_usage_percent",
        warning_threshold=80,
        critical_threshold=95,
        comparison="gt",
        unit="percent"
    ),
    "error_rate_percent": HealthThreshold(
        metric_name="error_rate_percent",
        warning_threshold=1,
        critical_threshold=5,
        comparison="gt",
        unit="percent"
    ),
    "token_validation_time_ms": HealthThreshold(
        metric_name="token_validation_time_ms",
        warning_threshold=50,
        critical_threshold=100,  # Per OmniPath spec
        comparison="gt",
        unit="milliseconds"
    ),
}

# Alert configurations
ALERT_CONFIGS: Dict[str, AlertConfig] = {
    "database_down": AlertConfig(
        name="Database Connection Failed",
        description="Cannot connect to PostgreSQL database",
        severity=AlertSeverity.CRITICAL,
        check_interval_seconds=30,
        notification_channels=["email", "slack", "pagerduty"]
    ),
    "ml_model_missing": AlertConfig(
        name="ML Model Not Loaded",
        description="One or more ML models failed to load",
        severity=AlertSeverity.WARNING,
        check_interval_seconds=60,
        notification_channels=["email", "slack"]
    ),
    "rdkit_unavailable": AlertConfig(
        name="RDKit Not Available",
        description="RDKit chemistry library not available",
        severity=AlertSeverity.WARNING,
        check_interval_seconds=300,
        notification_channels=["email"]
    ),
    "trade_secret_engine_down": AlertConfig(
        name="Trade Secret Engine Unavailable",
        description="One or more trade secret engines not responding",
        severity=AlertSeverity.CRITICAL,
        check_interval_seconds=60,
        notification_channels=["email", "slack", "pagerduty"]
    ),
    "high_error_rate": AlertConfig(
        name="High Error Rate",
        description="API error rate exceeds threshold",
        severity=AlertSeverity.WARNING,
        check_interval_seconds=60,
        notification_channels=["email", "slack"]
    ),
    "slow_response_time": AlertConfig(
        name="Slow API Response",
        description="API response time exceeds threshold",
        severity=AlertSeverity.WARNING,
        check_interval_seconds=60,
        notification_channels=["email"]
    ),
}

# Trade secret engine health checks
TRADE_SECRET_HEALTH_CHECKS = {
    "chempath": {
        "endpoint": "/api/v1/chempath/health",
        "timeout_seconds": 5,
        "required_for_healthy": False
    },
    "toxpath": {
        "endpoint": "/api/v1/toxpath/health",
        "timeout_seconds": 5,
        "required_for_healthy": False
    },
    "regpath": {
        "endpoint": "/api/v1/regpath/health",
        "timeout_seconds": 5,
        "required_for_healthy": False
    },
    "genomepath": {
        "endpoint": "/api/genomepath/health",
        "timeout_seconds": 10,
        "required_for_healthy": True  # Premium tier critical
    },
    "biopath": {
        "endpoint": "/api/biopath/health",
        "timeout_seconds": 10,
        "required_for_healthy": True  # Premium tier critical
    },
    "clinpath": {
        "endpoint": "/api/clinpath/health",
        "timeout_seconds": 10,
        "required_for_healthy": True  # Premium tier critical
    },
}

# ML model health checks
ML_MODEL_HEALTH_CHECKS = {
    "dimer_predictor": {
        "module": "backend.services.dimer_predictor",
        "class": "DimericPredictor",
        "required_for_healthy": True
    },
    "nextgen_dimer_model": {
        "model_path": "models/dimer_potential_v1.joblib",
        "required_for_healthy": False  # Fallback to legacy
    },
    "therapeutic_model": {
        "module": "backend.services.ml_models",
        "class": "TherapeuticPredictionModel",
        "required_for_healthy": False
    },
    "patient_response_model": {
        "module": "backend.services.ml_models",
        "class": "PatientResponseModel",
        "required_for_healthy": False
    },
}

# SLA configuration for Nevada pilot
NEVADA_PILOT_SLA = {
    "uptime_target_percent": 99.5,
    "max_response_time_p95_ms": 500,
    "max_response_time_p99_ms": 1000,
    "max_error_rate_percent": 0.5,
    "required_trade_secrets": ["chempath", "toxpath", "regpath"],
    "monitoring_hours": "24/7",
    "support_response_time_minutes": {
        "critical": 15,
        "warning": 60,
        "info": 480
    }
}

# Notification channel configuration (templates)
NOTIFICATION_CHANNELS = {
    "email": {
        "enabled": True,
        "recipients": [
            "ops@neurobotanica.com",
            "cto@cloakandquill.org"
        ],
        "template": "alert_email"
    },
    "slack": {
        "enabled": True,
        "webhook_url_env": "SLACK_WEBHOOK_URL",
        "channel": "#neurobotanica-alerts"
    },
    "pagerduty": {
        "enabled": False,  # Enable for production
        "integration_key_env": "PAGERDUTY_KEY"
    }
}

# Metrics to collect for dashboard
DASHBOARD_METRICS = [
    "requests_per_minute",
    "average_response_time_ms",
    "error_rate_percent",
    "active_users",
    "database_connections",
    "ml_predictions_per_minute",
    "trade_secret_calls_per_minute",
    "token_validations_per_minute",
    "memory_usage_mb",
    "cpu_usage_percent"
]

# Logging configuration
LOGGING_CONFIG = {
    "level": "INFO",
    "format": "%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    "handlers": ["console", "file"],
    "file_path": "logs/neurobotanica.log",
    "max_bytes": 10_000_000,  # 10 MB
    "backup_count": 5,
    "log_requests": True,
    "log_responses": False,  # Don't log response bodies (may contain sensitive data)
    "log_trade_secret_access": True,  # Audit trail for trade secrets
}
