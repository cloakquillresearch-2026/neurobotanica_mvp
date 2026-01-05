from fastapi import FastAPI

app = FastAPI(
    title="NeuroBotanica API",
    description="Dimeric Cannabinoid Therapeutic Prediction System",
    version="0.4.0"
)

@app.get("/")
async def root():
    """Root endpoint with API information."""
    return {"message": "NeuroBotanica API", "version": "0.4.0"}