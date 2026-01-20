"""
Unified API Layer for NeuroBotanica Engines
Combines all components for Budtender integration.
Cultural preservation: TK checks and community attribution.
"""

from fastapi import FastAPI, HTTPException, Depends, Header
from pydantic import BaseModel
from src.engines.interactions import DrugInteractionChecker
from src.engines.bias_correction import DemographicBiasCorrection
from src.engines.synergy import SynergyPredictionSystem
from src.engines.whole_plant import WholePlantAnalysisEngine
from src.engines.polysaccharides import PolysaccharideIntegrationEngine
import time
import json

app = FastAPI(title="NeuroBotanica API", version="1.0.0")

# Mock KV cache (replace with Cloudflare KV in prod)
cache = {}

# Dependency injection for engines
def get_engines():
    # Create engines fresh per request to avoid threading issues
    return {
        "interactions": DrugInteractionChecker(),
        "bias": DemographicBiasCorrection(),
        "synergy": SynergyPredictionSystem(),
        "plant": WholePlantAnalysisEngine(),
        "polysaccharides": PolysaccharideIntegrationEngine()
    }

class AnalysisRequest(BaseModel):
    compound_ids: list[str]
    demographics: dict = {}
    customer_tier: str = "computational_only"
    plant_id: str = None

@app.post("/api/neurobotanica/analyze")
async def analyze(request: AnalysisRequest, consent_header: str = Header(None, alias="X-Consent-ID"), db=None):
    """
    End-to-end analysis combining all engines.
    """
    try:
        start_time = time.time()
        
        # TK Handling: Check consent header
        if request.customer_tier == "tk_enhanced" and not consent_header:
            raise HTTPException(status_code=403, detail="TK access requires consent header")
        
        # Check cache
        cache_key = json.dumps(request.dict(), sort_keys=True)
        if cache_key in cache:
            return cache[cache_key]
        
        # Create engines per request to avoid threading issues
        engines = {
            "interactions": DrugInteractionChecker(db),
            "bias": DemographicBiasCorrection(db),
            "synergy": SynergyPredictionSystem(db),
            "plant": WholePlantAnalysisEngine(db),
            "polysaccharides": PolysaccharideIntegrationEngine(db)
        }
        
        # Phase 1: Interactions and Bias
        interactions = engines["interactions"].check_interactions(request.compound_ids, [], request.customer_tier)
        bias_corrected = engines["bias"].apply_corrections(10.0, request.compound_ids[0], request.demographics, request.customer_tier)
        
        # Phase 2: Synergy, Plant, Polysaccharides
        synergy = engines["synergy"].predict_synergy(request.compound_ids[0], request.compound_ids[1] if len(request.compound_ids) > 1 else request.compound_ids[0], request.customer_tier)
        plant_profile = engines["plant"].analyze_whole_plant(request.plant_id or "cbd", "anxiety", request.customer_tier) if request.plant_id else {}
        microbiome = {"bifidobacteria": 100}  # Mock
        polysaccharide_effects = engines["polysaccharides"].predict_cross_kingdom_effects("beta_glucan_001", microbiome, request.customer_tier)
        
        # Combine results
        result = {
            "interactions": interactions,
            "bias_correction": bias_corrected,
            "synergy": synergy,
            "plant_profile": plant_profile,
            "polysaccharide_effects": polysaccharide_effects,
            "processing_time_ms": (time.time() - start_time) * 1000
        }
        
        # Cache result
        cache[cache_key] = result
        
        # Close engines
        for engine in engines.values():
            engine.close()
        
        return result
    except HTTPException:
        # Re-raise HTTP exceptions (like 403 for consent)
        raise
    except Exception as e:
        # Trigger policy halt on TK failures
        if "consent" in str(e).lower():
            # Log to omnipath_policy_rules (mock)
            print(f"Policy violation: {e}")
        raise HTTPException(status_code=500, detail=str(e))

# Example endpoint for Budtender call
@app.get("/api/neurobotanica/health")
async def health_check():
    print("Health check called", flush=True)
    return {"status": "healthy", "engines": ["interactions", "bias", "synergy", "plant", "polysaccharides"]}

# Cloudflare Worker interface - commented out for Railway deployment
# async def fetch(request, env, ctx):
#     """Cloudflare Worker fetch handler for Python Workers."""
#     from fastapi.responses import JSONResponse
#     import json
    
#     # Get the request path and method
#     url = request.url
#     method = request.method
    
#     # Parse the path
#     path = url.path
    
#     if path == "/api/neurobotanica/health" and method == "GET":
#         return JSONResponse({"status": "healthy", "engines": ["interactions", "bias", "synergy", "plant", "polysaccharides"]})
    
#     elif path == "/api/neurobotanica/analyze" and method == "POST":
#         try:
#             # Parse JSON body
#             body = await request.json()
            
#             # Create request object for FastAPI
#             from pydantic import BaseModel
#             class MockRequest(BaseModel):
#                 compound_ids: list
#                 demographics: dict = {}
#                 customer_tier: str = "computational_only"
#                 plant_id: str = None
            
#             analysis_request = MockRequest(**body)
            
#             # Get headers
#             headers = dict(request.headers)
#             consent_header = headers.get("x-consent-id")
            
#             # Call the analysis function
#             result = await analyze(analysis_request, consent_header, env.DB)
            
#             return JSONResponse(result)
#         except Exception as e:
#             return JSONResponse({"error": str(e)}, status_code=500)
    
#     else:
#         return JSONResponse({"error": "Not found"}, status_code=404)