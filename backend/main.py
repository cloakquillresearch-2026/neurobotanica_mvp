"""
NeuroBotanica API - Main Application
Dimeric Cannabinoid Therapeutic Prediction System

FastAPI application with:
- 320-study evidence database access
- FDA compliance documentation generation
- Pharmacology data package APIs
- Comparative efficacy analysis
- 3D Conformer generation (ETKDG method)
- Patient & Treatment management
- OmniPath integration for provenance tracking
- Clinical Evidence API with confidence weighting
- Receptor Affinity with provenance and heterogeneity analysis
"""
from fastapi import FastAPI, Depends, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from sqlalchemy.orm import Session
from sqlalchemy import func
from contextlib import asynccontextmanager
import os

from backend.models.database import get_db, init_db, engine, Base
from backend.api import studies, compounds, fda_compliance, conformers
from backend.api import omnipath
from backend.api import evidence, receptor_affinity
from backend.api import dimers
from backend.routers import patentpath, terpenes, chempath, toxpath, regpath, security
from backend.middleware.token_validation import TokenValidationMiddleware


@asynccontextmanager
async def lifespan(app: FastAPI):
    """Application lifespan - initialize database on startup."""
    print("🌿 NeuroBotanica API starting...")
    init_db()
    print("✅ Database initialized")
    yield
    print("🛑 NeuroBotanica API shutting down...")


app = FastAPI(
    title="NeuroBotanica API",
    description="""
    ## Dimeric Cannabinoid Therapeutic Prediction System
    
    **Patent Claims Support:**
    - FDA Schedule III compliance documentation (Claim 1j)
    - Pharmacology data packages with mechanism of action
    - Comparative efficacy analysis from 320-study database
    - CMC templates for cannabis pharmaceutical approval
    - 3D conformer generation for molecular modeling
    - Clinical evidence aggregation with confidence weighting
    - Receptor affinity data with full provenance tracking
    
    **Data Assets:**
    - 320 clinical studies across 16 conditions
    - 63 cannabinoid compounds with molecular descriptors
    - 3D conformer ensembles (ETKDG method)
    - Dimeric triangulation predictions
    - FDA-approved drug evidence (Epidiolex, Marinol, Cesamet, Sativex)
    - Receptor affinity database with assay-level provenance
    """,
    version="0.2.0",
    lifespan=lifespan
)

# CORS middleware for frontend access
app.add_middleware(
    CORSMiddleware,
    allow_origins=[
        "http://localhost:3000",
        "http://localhost:5173",  # Vite dev server
        "https://neurobotanica.cloakandquill.org"
    ],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


# Include routers
app.include_router(studies.router, prefix="/api/v1/studies", tags=["Clinical Studies"])
app.include_router(compounds.router, prefix="/api/v1/compounds", tags=["Cannabinoid Compounds"])
app.include_router(fda_compliance.router, prefix="/api/v1/fda", tags=["FDA Compliance"])
app.include_router(conformers.router, prefix="/api/v1", tags=["3D Conformers"])
app.include_router(omnipath.router, prefix="/api/v1", tags=["OmniPath Integration"])
app.include_router(evidence.router, prefix="/api/v1", tags=["Clinical Evidence"])
app.include_router(receptor_affinity.router, prefix="/api/v1", tags=["Receptor Affinity"])
app.include_router(dimers.router, tags=["Dimers"])
app.include_router(patentpath.router, prefix="/api/v1", tags=["PatentPath Lite"])
app.include_router(terpenes.router, prefix="/api/v1", tags=["Terpene Analysis"])
app.include_router(chempath.router, tags=["ChemPath"])
app.include_router(toxpath.router, tags=["ToxPath"])
app.include_router(regpath.router, tags=["RegPath"])
app.include_router(security.router, tags=["Security"])

# Add token validation middleware (disabled by default for development)
# Enable with NEUROBOTANICA_TOKEN_VALIDATION=true environment variable
if os.getenv("NEUROBOTANICA_TOKEN_VALIDATION", "false").lower() == "true":
    app.add_middleware(TokenValidationMiddleware, validate_tokens=True)
else:
    # Add middleware but don't validate (allows testing without auth)
    app.add_middleware(TokenValidationMiddleware, validate_tokens=False)


@app.get("/")
async def root():
    """Root endpoint with API information."""
    return {
        "message": "🌿 NeuroBotanica API",
        "version": "0.1.0",
        "description": "Dimeric Cannabinoid Therapeutic Prediction System",
        "patent_claims": ["1j - FDA Schedule III Compliance Support"],
        "data_assets": {
            "clinical_studies": 320,
            "conditions": 16,
            "compounds": 63,
            "fda_approved_drugs": ["Epidiolex", "Marinol", "Cesamet", "Sativex"]
        },
        "endpoints": {
            "studies": "/api/v1/studies",
            "compounds": "/api/v1/compounds",
            "fda_compliance": "/api/v1/fda",
            "conformers": "/api/v1/conformers",
            "omnipath": "/api/v1/omnipath",
            "docs": "/docs"
        }
    }


@app.get("/health")
async def health_check(db: Session = Depends(get_db)):
    """Health check endpoint with database connectivity test."""
    from sqlalchemy import text
    try:
        # Test database connection
        db.execute(text("SELECT 1"))
        db_status = "connected"
    except Exception as e:
        db_status = f"error: {str(e)}"
    
    return {
        "status": "healthy" if db_status == "connected" else "degraded",
        "database": db_status,
        "version": "0.1.0",
        "features": {
            "fda_compliance_module": True,
            "pharmacology_packages": True,
            "comparative_efficacy": True,
            "schedule_iii_support": True,
            "conformer_generation": True,
            "patient_treatment_models": True,
            "omnipath_integration": True,
            "provenance_tracking": True,
            "token_validation": os.getenv("NEUROBOTANICA_TOKEN_VALIDATION", "false").lower() == "true"
        }
    }


@app.get("/api/v1/stats")
async def get_stats(db: Session = Depends(get_db)):
    """Get database statistics."""
    from backend.models.study import ClinicalStudy
    from backend.models.compound import Cannabinoid, DimericPrediction
    
    try:
        study_count = db.query(ClinicalStudy).count()
        compound_count = db.query(Cannabinoid).count()
        dimer_count = db.query(DimericPrediction).count()
        
        # Get condition breakdown
        conditions = db.query(
            ClinicalStudy.condition,
            func.count(ClinicalStudy.id)
        ).group_by(ClinicalStudy.condition).all()
        
        return {
            "total_studies": study_count,
            "total_compounds": compound_count,
            "dimeric_predictions": dimer_count,
            "conditions": {c[0]: c[1] for c in conditions},
            "fda_approved_coverage": ["Epidiolex", "Marinol", "Cesamet", "Sativex"]
        }
    except Exception as e:
        return {
            "total_studies": 0,
            "total_compounds": 0,
            "dimeric_predictions": 0,
            "message": "Database not yet populated",
            "note": "Run data loading scripts to populate"
        }


if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000, reload=True)
