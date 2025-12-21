# Services Layer
from .efficacy_analyzer import ComparativeEfficacyAnalyzer
from .schedule_iii_docs import ScheduleIIIDocumentationGenerator

# Conformer generator requires RDKit
try:
    from .conformer_generator import (
        ConformerGenerator,
        BatchConformerGenerator,
        ConformerResult,
        generate_conformers_for_smiles,
        calculate_3d_descriptors
    )
    CONFORMER_AVAILABLE = True
except ImportError:
    CONFORMER_AVAILABLE = False

# OmniPath integration
from .omnipath_client import (
    OmniPathClient,
    OmniPathToken,
    ProvenanceManifest,
    ConsentLevel,
    OperationType,
    get_omnipath_client,
    reset_omnipath_client
)

# Provenance tracking
from .provenance_tracker import (
    ProvenanceTracker,
    AuditEntry,
    get_provenance_tracker,
    reset_provenance_tracker,
    track_provenance,
    tracked_transaction
)

__all__ = [
    "ComparativeEfficacyAnalyzer",
    "ScheduleIIIDocumentationGenerator",
    "ConformerGenerator",
    "BatchConformerGenerator",
    "ConformerResult",
    "generate_conformers_for_smiles",
    "calculate_3d_descriptors",
    "CONFORMER_AVAILABLE",
    # OmniPath
    "OmniPathClient",
    "OmniPathToken",
    "ProvenanceManifest",
    "ConsentLevel",
    "OperationType",
    "get_omnipath_client",
    "reset_omnipath_client",
    # Provenance
    "ProvenanceTracker",
    "AuditEntry",
    "get_provenance_tracker",
    "reset_provenance_tracker",
    "track_provenance",
    "tracked_transaction"
]
