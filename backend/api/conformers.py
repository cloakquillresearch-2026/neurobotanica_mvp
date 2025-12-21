"""
Conformer API Routes
3D Conformer generation and retrieval endpoints

Patent Claims Support:
- 3D molecular descriptor calculation
- Conformer ensemble generation for receptor docking
- Batch processing for drug discovery pipelines
"""
from fastapi import APIRouter, Depends, HTTPException, Query, BackgroundTasks
from sqlalchemy.orm import Session
from pydantic import BaseModel, Field
from typing import List, Optional, Dict, Any
from datetime import datetime

from ..models.database import get_db
from ..models.compound import Cannabinoid

# Import conformer generator
try:
    from ..services.conformer_generator import (
        ConformerGenerator,
        BatchConformerGenerator,
        ConformerResult
    )
    CONFORMER_AVAILABLE = True
except ImportError:
    CONFORMER_AVAILABLE = False


router = APIRouter(prefix="/conformers", tags=["3D Conformers"])


# Pydantic models for request/response
class ConformerRequest(BaseModel):
    """Request to generate conformers."""
    num_conformers: int = Field(default=50, ge=1, le=500)
    force_regenerate: bool = Field(default=False)


class ConformerFromSMILES(BaseModel):
    """Request to generate conformers from SMILES."""
    smiles: str = Field(..., min_length=1)
    num_conformers: int = Field(default=50, ge=1, le=500)


class Descriptors3DResponse(BaseModel):
    """3D descriptor response."""
    pmi1: Optional[float] = None
    pmi2: Optional[float] = None
    pmi3: Optional[float] = None
    npr1: Optional[float] = None
    npr2: Optional[float] = None
    radius_of_gyration: Optional[float] = None
    inertial_shape_factor: Optional[float] = None
    asphericity: Optional[float] = None
    eccentricity: Optional[float] = None
    spherocity_index: Optional[float] = None
    pbf: Optional[float] = None


class ConformerResponse(BaseModel):
    """Conformer generation response."""
    success: bool
    compound_id: Optional[int] = None
    compound_name: Optional[str] = None
    num_conformers_generated: int = 0
    lowest_energy_kcal_mol: Optional[float] = None
    energy_range_kcal_mol: Optional[float] = None
    num_clusters_rmsd_2A: int = 0
    descriptors_3d: Optional[Dict[str, Any]] = None
    conformer_metadata: Optional[Dict[str, Any]] = None
    generation_time_seconds: float = 0.0
    error: Optional[str] = None


class ConformerStatusResponse(BaseModel):
    """Conformer status for a compound."""
    compound_id: int
    compound_name: str
    has_conformers: bool
    generation_method: Optional[str] = None
    num_conformers: int = 0
    generated_at: Optional[datetime] = None
    descriptors_3d: Optional[Dict[str, Any]] = None
    conformer_metadata: Optional[Dict[str, Any]] = None


class BatchConformerResponse(BaseModel):
    """Batch conformer generation response."""
    total_requested: int
    successful: int
    failed: int
    success_rate: float
    processing_time_seconds: float
    results: List[Dict[str, Any]]


def check_rdkit_available():
    """Dependency to check if RDKit is available."""
    if not CONFORMER_AVAILABLE:
        raise HTTPException(
            status_code=503,
            detail="RDKit not available. Please install RDKit to use conformer generation."
        )


@router.get("/status", response_model=Dict[str, Any])
async def get_conformer_service_status():
    """Get conformer service status and statistics."""
    return {
        "service": "3D Conformer Generator",
        "rdkit_available": CONFORMER_AVAILABLE,
        "generation_method": "ETKDG" if CONFORMER_AVAILABLE else None,
        "optimization_method": "MMFF94" if CONFORMER_AVAILABLE else None,
        "default_conformers": 50,
        "max_conformers": 500,
        "supported_operations": [
            "generate_for_compound",
            "generate_from_smiles",
            "batch_generate",
            "get_3d_descriptors"
        ] if CONFORMER_AVAILABLE else []
    }


@router.post("/generate/{compound_id}", response_model=ConformerResponse)
async def generate_conformers_for_compound(
    compound_id: int,
    request: ConformerRequest = ConformerRequest(),
    db: Session = Depends(get_db),
    _: None = Depends(check_rdkit_available)
):
    """Generate 3D conformers for a specific compound.
    
    This endpoint generates a conformer ensemble using the ETKDG method
    with MMFF94 force field optimization. The lowest energy conformer
    is used for 3D descriptor calculation.
    
    Args:
        compound_id: Database ID of the compound
        request: Generation parameters
        
    Returns:
        ConformerResponse with generation results and 3D descriptors
    """
    # Get compound
    compound = db.query(Cannabinoid).filter(Cannabinoid.id == compound_id).first()
    
    if not compound:
        raise HTTPException(status_code=404, detail=f"Compound with ID {compound_id} not found")
    
    # Check if already has conformers
    if compound.has_conformers and not request.force_regenerate:
        return ConformerResponse(
            success=True,
            compound_id=compound.id,
            compound_name=compound.name,
            num_conformers_generated=compound.num_conformers_generated or 0,
            lowest_energy_kcal_mol=compound.conformer_metadata.get("lowest_energy_kcal_mol") if compound.conformer_metadata else None,
            energy_range_kcal_mol=compound.conformer_metadata.get("energy_range_kcal_mol") if compound.conformer_metadata else None,
            num_clusters_rmsd_2A=compound.conformer_metadata.get("num_clusters_rmsd_2A", 0) if compound.conformer_metadata else 0,
            descriptors_3d=compound.rdkit_descriptors_3d,
            conformer_metadata=compound.conformer_metadata,
            generation_time_seconds=0.0,
            error="Conformers already exist. Use force_regenerate=true to regenerate."
        )
    
    # Generate conformers
    generator = ConformerGenerator(num_conformers=request.num_conformers)
    result = generator.generate_conformers(compound.smiles)
    
    if result.success:
        # Update database
        compound.has_conformers = True
        compound.conformer_generation_method = "ETKDG"
        compound.num_conformers_generated = result.num_conformers_generated
        compound.conformer_metadata = result.conformer_metadata
        compound.rdkit_descriptors_3d = result.descriptors_3d
        compound.conformers_generated_at = datetime.now()
        
        db.commit()
        db.refresh(compound)
        
        return ConformerResponse(
            success=True,
            compound_id=compound.id,
            compound_name=compound.name,
            num_conformers_generated=result.num_conformers_generated,
            lowest_energy_kcal_mol=result.lowest_energy,
            energy_range_kcal_mol=result.energy_range,
            num_clusters_rmsd_2A=result.num_clusters_rmsd_2A,
            descriptors_3d=result.descriptors_3d,
            conformer_metadata=result.conformer_metadata,
            generation_time_seconds=result.generation_time_seconds
        )
    else:
        return ConformerResponse(
            success=False,
            compound_id=compound.id,
            compound_name=compound.name,
            error=result.error,
            generation_time_seconds=result.generation_time_seconds
        )


@router.post("/generate-from-smiles", response_model=ConformerResponse)
async def generate_conformers_from_smiles(
    request: ConformerFromSMILES,
    _: None = Depends(check_rdkit_available)
):
    """Generate 3D conformers from a SMILES string.
    
    This endpoint generates conformers without database storage.
    Useful for ad-hoc analysis of new compounds.
    
    Args:
        request: SMILES string and generation parameters
        
    Returns:
        ConformerResponse with generation results
    """
    generator = ConformerGenerator(num_conformers=request.num_conformers)
    result = generator.generate_conformers(request.smiles)
    
    if result.success:
        return ConformerResponse(
            success=True,
            num_conformers_generated=result.num_conformers_generated,
            lowest_energy_kcal_mol=result.lowest_energy,
            energy_range_kcal_mol=result.energy_range,
            num_clusters_rmsd_2A=result.num_clusters_rmsd_2A,
            descriptors_3d=result.descriptors_3d,
            conformer_metadata=result.conformer_metadata,
            generation_time_seconds=result.generation_time_seconds
        )
    else:
        return ConformerResponse(
            success=False,
            error=result.error,
            generation_time_seconds=result.generation_time_seconds
        )


@router.get("/compound/{compound_id}", response_model=ConformerStatusResponse)
async def get_conformer_status(
    compound_id: int,
    db: Session = Depends(get_db)
):
    """Get conformer status for a specific compound.
    
    Args:
        compound_id: Database ID of the compound
        
    Returns:
        ConformerStatusResponse with conformer data if available
    """
    compound = db.query(Cannabinoid).filter(Cannabinoid.id == compound_id).first()
    
    if not compound:
        raise HTTPException(status_code=404, detail=f"Compound with ID {compound_id} not found")
    
    return ConformerStatusResponse(
        compound_id=compound.id,
        compound_name=compound.name,
        has_conformers=compound.has_conformers or False,
        generation_method=compound.conformer_generation_method,
        num_conformers=compound.num_conformers_generated or 0,
        generated_at=compound.conformers_generated_at,
        descriptors_3d=compound.rdkit_descriptors_3d,
        conformer_metadata=compound.conformer_metadata
    )


@router.get("/compound/{compound_id}/descriptors-3d", response_model=Descriptors3DResponse)
async def get_3d_descriptors(
    compound_id: int,
    db: Session = Depends(get_db)
):
    """Get 3D molecular descriptors for a compound.
    
    Returns PMI, asphericity, eccentricity, and other 3D shape descriptors.
    
    Args:
        compound_id: Database ID of the compound
        
    Returns:
        Descriptors3DResponse with 3D shape descriptors
    """
    compound = db.query(Cannabinoid).filter(Cannabinoid.id == compound_id).first()
    
    if not compound:
        raise HTTPException(status_code=404, detail=f"Compound with ID {compound_id} not found")
    
    if not compound.has_conformers or not compound.rdkit_descriptors_3d:
        raise HTTPException(
            status_code=404,
            detail=f"3D descriptors not available for compound {compound.name}. Generate conformers first."
        )
    
    desc = compound.rdkit_descriptors_3d
    return Descriptors3DResponse(
        pmi1=desc.get("pmi1"),
        pmi2=desc.get("pmi2"),
        pmi3=desc.get("pmi3"),
        npr1=desc.get("npr1"),
        npr2=desc.get("npr2"),
        radius_of_gyration=desc.get("radius_of_gyration"),
        inertial_shape_factor=desc.get("inertial_shape_factor"),
        asphericity=desc.get("asphericity"),
        eccentricity=desc.get("eccentricity"),
        spherocity_index=desc.get("spherocity_index"),
        pbf=desc.get("pbf")
    )


@router.get("/summary", response_model=Dict[str, Any])
async def get_conformer_summary(db: Session = Depends(get_db)):
    """Get summary of conformer generation status across all compounds.
    
    Returns:
        Summary statistics including counts and coverage
    """
    total_compounds = db.query(Cannabinoid).count()
    with_conformers = db.query(Cannabinoid).filter(
        Cannabinoid.has_conformers == True
    ).count()
    without_conformers = total_compounds - with_conformers
    
    # Get average descriptors for compounds with conformers
    compounds_with_conf = db.query(Cannabinoid).filter(
        Cannabinoid.has_conformers == True
    ).all()
    
    avg_conformers = 0
    avg_asphericity = 0
    avg_energy_range = 0
    
    if compounds_with_conf:
        conformer_counts = [c.num_conformers_generated or 0 for c in compounds_with_conf]
        avg_conformers = sum(conformer_counts) / len(conformer_counts)
        
        asphericities = [
            c.rdkit_descriptors_3d.get("asphericity", 0) 
            for c in compounds_with_conf 
            if c.rdkit_descriptors_3d
        ]
        if asphericities:
            avg_asphericity = sum(asphericities) / len(asphericities)
        
        energy_ranges = [
            c.conformer_metadata.get("energy_range_kcal_mol", 0)
            for c in compounds_with_conf
            if c.conformer_metadata
        ]
        if energy_ranges:
            avg_energy_range = sum(energy_ranges) / len(energy_ranges)
    
    return {
        "total_compounds": total_compounds,
        "with_conformers": with_conformers,
        "without_conformers": without_conformers,
        "coverage_percent": round(with_conformers / total_compounds * 100, 1) if total_compounds > 0 else 0,
        "average_conformers_per_compound": round(avg_conformers, 1),
        "average_asphericity": round(avg_asphericity, 4),
        "average_energy_range_kcal_mol": round(avg_energy_range, 2),
        "rdkit_available": CONFORMER_AVAILABLE
    }


@router.get("/pending", response_model=List[Dict[str, Any]])
async def get_pending_compounds(
    limit: int = Query(default=100, le=1000),
    db: Session = Depends(get_db)
):
    """Get list of compounds pending conformer generation.
    
    Args:
        limit: Maximum number of compounds to return
        
    Returns:
        List of compounds without conformers
    """
    compounds = db.query(Cannabinoid).filter(
        Cannabinoid.has_conformers == False
    ).limit(limit).all()
    
    return [
        {
            "id": c.id,
            "name": c.name,
            "abbreviation": c.abbreviation,
            "smiles": c.smiles,
            "molecular_weight": c.molecular_weight,
            "is_dimeric": c.is_dimeric
        }
        for c in compounds
    ]


@router.post("/batch-generate", response_model=BatchConformerResponse)
async def batch_generate_conformers(
    compound_ids: List[int],
    num_conformers: int = Query(default=50, ge=1, le=500),
    db: Session = Depends(get_db),
    _: None = Depends(check_rdkit_available)
):
    """Generate conformers for multiple compounds.
    
    Note: This is a synchronous operation. For large batches,
    consider using the background task or CLI script.
    
    Args:
        compound_ids: List of compound database IDs
        num_conformers: Number of conformers per compound
        
    Returns:
        BatchConformerResponse with results for each compound
    """
    import time
    start_time = time.time()
    
    # Get compounds
    compounds = db.query(Cannabinoid).filter(
        Cannabinoid.id.in_(compound_ids)
    ).all()
    
    if not compounds:
        raise HTTPException(status_code=404, detail="No compounds found with the provided IDs")
    
    generator = ConformerGenerator(num_conformers=num_conformers)
    results = []
    successful = 0
    failed = 0
    
    for compound in compounds:
        result = generator.generate_conformers(compound.smiles)
        
        if result.success:
            # Update database
            compound.has_conformers = True
            compound.conformer_generation_method = "ETKDG"
            compound.num_conformers_generated = result.num_conformers_generated
            compound.conformer_metadata = result.conformer_metadata
            compound.rdkit_descriptors_3d = result.descriptors_3d
            compound.conformers_generated_at = datetime.now()
            
            successful += 1
            results.append({
                "compound_id": compound.id,
                "compound_name": compound.name,
                "success": True,
                "num_conformers": result.num_conformers_generated,
                "time_seconds": result.generation_time_seconds
            })
        else:
            failed += 1
            results.append({
                "compound_id": compound.id,
                "compound_name": compound.name,
                "success": False,
                "error": result.error,
                "time_seconds": result.generation_time_seconds
            })
    
    db.commit()
    
    processing_time = time.time() - start_time
    
    return BatchConformerResponse(
        total_requested=len(compound_ids),
        successful=successful,
        failed=failed,
        success_rate=round(successful / len(compounds) * 100, 1) if compounds else 0,
        processing_time_seconds=round(processing_time, 2),
        results=results
    )
