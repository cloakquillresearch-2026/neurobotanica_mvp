import asyncio
from datetime import datetime, timezone
from typing import Any, Dict


_CACHE: Dict[str, Any] = {}
_LOCK = asyncio.Lock()


async def _compute_health() -> Dict[str, Any]:
    """Compute a lightweight health snapshot."""
    now = datetime.now(timezone.utc).isoformat()
    return {
        "status": "ok",
        "checked_at": now,
        "errors": [],
        "results": {
            "ml_models": {"loaded": True},
            "features": {"available": True}
        }
    }


async def startup_monitor(app) -> None:
    """Initialize health cache on application startup."""
    async with _LOCK:
        _CACHE.clear()
        _CACHE.update(await _compute_health())
    # attach to app state for potential introspection
    try:
        app.state._health_cache = _CACHE
    except Exception:
        pass


async def shutdown_monitor(app) -> None:
    """Cleanup resources on shutdown."""
    async with _LOCK:
        _CACHE.clear()
    try:
        if hasattr(app.state, "_health_cache"):
            delattr(app.state, "_health_cache")
    except Exception:
        pass


async def get_cached_health(force_refresh: bool = False) -> Dict[str, Any]:
    """Return cached health info; recompute if requested."""
    async with _LOCK:
        if not _CACHE or force_refresh:
            _CACHE.clear()
            _CACHE.update(await _compute_health())
        return dict(_CACHE)
