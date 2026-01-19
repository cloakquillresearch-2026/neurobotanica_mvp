"""Lightweight recommendations router (placeholder).

This module provides a minimal `router` used by the main app and tests.
It intentionally keeps logic small so tests that import the app during
collection will succeed. More advanced recommendation logic lives in
`backend.routers.dispensary` and related services.
"""
from fastapi import APIRouter

router = APIRouter(prefix="/api/recommendations", tags=["Recommendations"])


@router.post("/generate")
async def generate_recommendations(payload: dict):
    """Return an empty recommendations structure for test/import-time safety.

    Tests that exercise recommendation endpoints will override/mock this
    behavior as needed; this keeps startup safe during collection.
    """
    return {"recommendations": [], "note": "placeholder"}
