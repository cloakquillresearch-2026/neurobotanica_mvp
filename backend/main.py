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
from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware


app = FastAPI(
    title="NeuroBotanica API",
    description="""
    ## Dimeric Cannabinoid Therapeutic Prediction System

    **Patent Claims Support:**
    - FDA Schedule III compliance documentation (Claim 1j)
    - Pharmacology data packages with mechanism of action
    - Comparative efficacy analysis from 368-study database
    - CMC templates for cannabis pharmaceutical approval
    - 3D conformer generation for molecular modeling
    - Clinical evidence aggregation with confidence weighting
    - Receptor affinity data with full provenance tracking
    - Trade Secret Engines: ChemPath, ToxPath, RegPath, BioPath, ClinPath, GenomePath

    **Data Assets:**
    - 368 clinical studies across 22+ conditions
    - 63 cannabinoid compounds with molecular descriptors
    - 10,084 dimer/entourage validated entries
    - 3D conformer ensembles (ETKDG method)
    - Dimeric triangulation predictions
    - FDA-approved drug evidence (Epidiolex, Marinol, Cesamet, Sativex)
    - Receptor affinity database with assay-level provenance
    """,
    version="0.4.0",
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


@app.get("/")
async def root():
    """Root endpoint returning API information."""
    return {
        "message": "🌿 NeuroBotanica API - Dimeric Cannabinoid Therapeutic Prediction System",
        "version": "0.4.0",
        "status": "operational",
        "endpoints": {
            "health": "/health",
            "api_docs": "/docs",
            "api_redoc": "/redoc"
        }
    }


@app.get("/health")
async def health_check():
    """Health check endpoint."""
    return {"status": "healthy", "timestamp": "2024-01-01T00:00:00Z"}


@app.get("/api/neurobotanica/health")
async def neurobotanica_health():
    """NeuroBotanica health check endpoint."""
    return {"status": "healthy"}


@app.post("/api/neurobotanica/analyze")
async def neurobotanica_analyze(data: dict):
    """NeuroBotanica analyze endpoint."""
    return {"interactions": []}
# app.include_router(omnipath.router, prefix="/api/v1", tags=["OmniPath Integration"])
# app.include_router(evidence.router, prefix="/api/v1", tags=["Clinical Evidence"])
# app.include_router(receptor_affinity.router, prefix="/api/v1", tags=["Receptor Affinity"])
# app.include_router(dimers.router, tags=["Dimers"])
# app.include_router(patentpath.router, prefix="/api/v1", tags=["PatentPath Lite"])
# app.include_router(terpenes.router, prefix="/api/v1", tags=["Terpene Analysis"])
# app.include_router(chempath.router, tags=["ChemPath"])
# app.include_router(toxpath.router, tags=["ToxPath"])
# app.include_router(regpath.router, tags=["RegPath"])
# app.include_router(genomepath.router, tags=["GenomePath"])  # Router already has /api/genomepath prefix
# app.include_router(biopath.router, prefix="/api/biopath", tags=["BioPath"])
# app.include_router(clinpath.router, prefix="/api/clinpath", tags=["ClinPath"])
# app.include_router(dispensary.router, prefix="/api/dispensary", tags=["Dispensary"])
# app.include_router(security.router, tags=["Security"])
# app.include_router(recommendations_router)

# Ensure database schema matches models at import time for test/dev environments
# try:
#     init_db()
# except Exception as e:
#     import logging
#     logging.getLogger(__name__).warning(f"init_db at import failed: {e}")


# @app.on_event("startup")
# async def _start_monitor():
#     # start background health monitor
#     try:
#         await health_monitor.startup_monitor(app)
#     except Exception as e:
#         import logging
#         logging.getLogger(__name__).warning(f"Failed to start health monitor: {e}")
#     # Optionally start JWKS refresher for Firebase token verification
#     try:
#         import logging
#         if os.getenv("ENABLE_JWKS_REFRESH", "false").lower() == "true":
#             from backend.auth import firebase

#             period = int(os.getenv("JWKS_REFRESH_PERIOD", "3600"))
#             try:
#                 await firebase.start_jwks_refresher(app, period)
#             except Exception as e:
#                 logging.getLogger(__name__).warning(f"Failed to start JWKS refresher: {e}")
#     except Exception:
#         pass


# @app.on_event("shutdown")
# async def _stop_monitor():
#     try:
#         await health_monitor.shutdown_monitor(app)
#     except Exception:
#         pass
#     # Stop JWKS refresher if it was started
#     try:
#         from backend.auth import firebase

#         try:
#             await firebase.stop_jwks_refresher(app)
#         except Exception:
#             pass
#     except Exception:
#         pass

# Add token validation middleware (disabled by default for development)
# Add token validation middleware (disabled by default for development)
# When running tests, force validation off to avoid 401s during collection.
# test_running = (
#     "PYTEST_CURRENT_TEST" in os.environ
#     or "PYTEST_ADDOPTS" in os.environ
#     or any(k.startswith("PYTEST") for k in os.environ)
# )
# if test_running:
#     # app.add_middleware(TokenValidationMiddleware, validate_tokens=False)
#     pass
# else:
#     if os.getenv("NEUROBOTANICA_TOKEN_VALIDATION", "false").lower() == "true":
#         # app.add_middleware(TokenValidationMiddleware, validate_tokens=True)
#         pass
#     else:
#         # Add middleware but don't validate (allows testing without auth)
#         # app.add_middleware(TokenValidationMiddleware, validate_tokens=False)
#         pass


# @app.get("/")
# async def root():
#     """Root endpoint with API information."""
#     return {
#         "message": "🌿 NeuroBotanica API",
#         "version": "0.4.0",
#         "description": "Dimeric Cannabinoid Therapeutic Prediction System"
#     }


# @app.get("/test")
# async def test_endpoint():
#     """Simple test endpoint without database dependency."""
#     try:
#         return {"status": "ok", "message": "Test endpoint working"}
#     except Exception as e:
#         import logging
#         logger = logging.getLogger(__name__)
#         logger.error(f"Error in test endpoint: {e}")
#         import traceback
#         traceback.print_exc()
#         return {"status": "error", "message": str(e)}


# @app.get("/health")
# async def health(refresh: bool = False):
#     """Return cached health data. Use ?refresh=true to force an immediate live check.
#     Adds top-level backward-compatible keys: `ml_models` and `features`.
#     """
#     try:
#         data = await health_monitor.get_cached_health(force_refresh=refresh)
#         # Backwards-compatible top-level keys expected by older tests/clients
#         response = data.copy()
#         response['ml_models'] = data.get('results', {}).get('ml_models')
#         response['features'] = data.get('results', {}).get('features')
#         return JSONResponse(response)
#     except Exception as e:
#         return JSONResponse({'status': 'error', 'detail': str(e)}, status_code=500)


# @app.get("/metrics")
# async def metrics():
#     """Prometheus metrics endpoint. Returns 501 if prometheus_client is not installed."""
#     try:
#         from prometheus_client import generate_latest, CONTENT_TYPE_LATEST
#     except Exception:
#         return JSONResponse({'status': 'unavailable', 'detail': 'prometheus_client not installed'}, status_code=501)

#     payload = generate_latest()
#     return PlainTextResponse(payload.decode('utf-8'), media_type=CONTENT_TYPE_LATEST)


# @app.get("/status")
# async def status():
#     """Simple HTML status UI showing basic health summary."""
#     data = await health_monitor.get_cached_health()
#     status_val = data.get('status', 'unknown')
#     checked_at = data.get('checked_at', '')
#     errors = data.get('errors', [])
#     html = f"""
#     <html>
#       <head><title>NeuroBotanica Status</title></head>
#       <body>
#         <h2>NeuroBotanica Status: {status_val}</h2>
#         <p>Last check: {checked_at}</p>
#         <p>Errors: {', '.join(errors) if errors else 'none'}</p>
#       </body>
#     </html>
#     """
#     return HTMLResponse(html)


# Internal debug endpoints (enable with ENABLE_JWKS_DEBUG=true)
# if os.getenv("ENABLE_JWKS_DEBUG", "false").lower() == "true":
#     from fastapi import APIRouter

#     debug_router = APIRouter(prefix="/internal/debug", tags=["internal"])

#     @debug_router.get("/jwks_cache")
#     async def jwks_cache():
#         try:
#             from backend.auth import firebase

#             snapshot = firebase.get_jwks_cache_snapshot()
#             return JSONResponse({"jwks_cache": snapshot})
#         except Exception as e:
#             return JSONResponse({"error": str(e)}, status_code=500)

#     app.include_router(debug_router)


# @app.get("/api/v1/stats")
# async def get_stats(db: Session = Depends(get_db)):
#     """Get database statistics."""
#     from backend.models.study import ClinicalStudy
#     from backend.models.compound import Cannabinoid, DimericPrediction
    
#     try:
#         study_count = db.query(ClinicalStudy).count()
#         compound_count = db.query(Cannabinoid).count()
#         dimer_count = db.query(DimericPrediction).count()
        
#         # Get condition breakdown
#         conditions = db.query(
#             ClinicalStudy.condition,
#             func.count(ClinicalStudy.id)
#         ).group_by(ClinicalStudy.condition).all()
        
#         return {
#             "total_studies": study_count,
#             "total_compounds": compound_count,
#             "dimeric_predictions": dimer_count,
#             "conditions": {c[0]: c[1] for c in conditions},
#             "fda_approved_coverage": ["Epidiolex", "Marinol", "Cesamet", "Sativex"]
#         }
#     except Exception as e:
#         return {
#             "total_studies": 0,
#             "total_compounds": 0,
#             "dimeric_predictions": 0,
#             "message": "Database not yet populated",
#             "note": "Run data loading scripts to populate"
#         }


if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000, reload=True)
