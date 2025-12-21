"""
ChemPath - Chemical Characterization Engine

Trade Secret Module for NeuroBotanica.
Exposes capabilities via API while keeping scoring logic/weights/heuristics proprietary.

Features:
- Molecular structure normalization and validation
- 2D/3D descriptor computation
- COA (Certificate of Analysis) QC validation
- Molecular characterization report generation
- Cannabis-specific cannabinoid/terpene analysis
- Entourage effect modeling
"""

from .analyzer import (
    ChemPathAnalyzer,
    CompoundInput,
    COAInput,
    ChemPathRequest,
    ChemPathResponse,
    QCFlag,
    QCSeverity,
    NormalizedStructure,
    DescriptorSet
)

from .report_generator import (
    ChemPathReportGenerator,
    ReportFormat
)

from .cannabinoid_analyzer import (
    CannabinoidAnalyzer,
    CannabinoidType,
    TerpeneType,
    ReceptorTarget,
    TherapeuticCategory,
    CannabinoidProfile,
    TerpeneProfile,
    SynergyInteraction,
    EntourageProfile,
    FormulationRecommendation,
)

__all__ = [
    # Main analyzer
    "ChemPathAnalyzer",
    
    # Input models
    "CompoundInput",
    "COAInput", 
    "ChemPathRequest",
    
    # Output models
    "ChemPathResponse",
    "QCFlag",
    "QCSeverity",
    "NormalizedStructure",
    "DescriptorSet",
    
    # Report generation
    "ChemPathReportGenerator",
    "ReportFormat",
    
    # Cannabis-specific analysis
    "CannabinoidAnalyzer",
    "CannabinoidType",
    "TerpeneType",
    "ReceptorTarget",
    "TherapeuticCategory",
    "CannabinoidProfile",
    "TerpeneProfile",
    "SynergyInteraction",
    "EntourageProfile",
    "FormulationRecommendation",
]
