import os
import time
from typing import Dict, Any, Callable

_JWKS_CACHE: dict[str, Any] = {}
_JWKS_CACHE_TTL = 60 * 60  # seconds

# Optional Prometheus metrics (lazy). Use importlib to avoid static import errors
_AUTH_VERIFY_ATTEMPTS = None
_AUTH_VERIFY_FAILURES = None
try:
    import importlib

    _prom = importlib.import_module("prometheus_client")
    Counter = getattr(_prom, "Counter")
    _AUTH_VERIFY_ATTEMPTS = Counter("neurobotanica_auth_verify_attempts_total", "Auth verify attempts")
    _AUTH_VERIFY_FAILURES = Counter("neurobotanica_auth_verify_failures_total", "Auth verify failures")
except Exception:
    _AUTH_VERIFY_ATTEMPTS = None
    _AUTH_VERIFY_FAILURES = None


def _get_jwks_url(project_id: str) -> str:
    return f"https://www.googleapis.com/service_accounts/v1/jwk/securetoken@system.gserviceaccount.com"


def verify_id_token(token: str) -> Dict[str, Any]:
    """Verify a Firebase ID token.

    Behavior:
    - If `FIREBASE_TEST_MODE=1` accepts `test-token-<role>` for local tests/CI.
    - If `firebase_admin` is available and configured via `FIREBASE_SERVICE_ACCOUNT_JSON`, use it.
    - Otherwise attempt JWKS verification using PyJWT/PyJWKClient and `FIREBASE_PROJECT_ID`.
    """
    # Test mode: lightweight deterministic verifier for unit tests
    if os.getenv("FIREBASE_TEST_MODE") == "1":
        if token.startswith("test-token-"):
            role = token.split("test-token-")[-1]
            return {"uid": f"{role}-uid", "email": f"{role}@example.com", "role": role}
        raise ValueError("Invalid test token")

    # Try firebase_admin if available (import dynamically to avoid static analysis errors)
    try:
        import importlib
        fb_mod = importlib.import_module("firebase_admin")
        fb_auth = getattr(fb_mod, "auth")

        if not getattr(fb_mod, "_apps", None):
            cred_json = os.getenv("FIREBASE_SERVICE_ACCOUNT_JSON")
            if cred_json:
                import json
                creds_mod = importlib.import_module("firebase_admin.credentials")
                cred = creds_mod.Certificate(json.loads(cred_json))
                fb_mod.initialize_app(cred)
            else:
                raise RuntimeError("FIREBASE_SERVICE_ACCOUNT_JSON not set; cannot initialize firebase_admin")

        claims = fb_auth.verify_id_token(token)
        return claims
    except Exception:
        # fallthrough to JWKS verification
        pass

    # JWKS verification path
    project_id = os.getenv("FIREBASE_PROJECT_ID")
    if not project_id:
        raise RuntimeError("FIREBASE_PROJECT_ID is required for JWKS-based token verification")

    try:
        # Lazy dynamic imports to avoid requiring PyJWT for tests using FIREBASE_TEST_MODE
        import importlib

        jwt_mod = importlib.import_module("jwt")
        PyJWKClient = getattr(jwt_mod, "PyJWKClient", None)

        if _AUTH_VERIFY_ATTEMPTS:
            try:
                _AUTH_VERIFY_ATTEMPTS.inc()
            except Exception:
                pass

        if PyJWKClient is None:
            raise RuntimeError("PyJWKClient not available; install PyJWT>=2.0.0")

        jwks_url = _get_jwks_url(project_id)

        # Simple JWKS client caching with TTL to reduce lookups
        entry = _JWKS_CACHE.get(jwks_url)
        now = time.time()
        if not entry or (now - entry.get("fetched_at", 0)) > _JWKS_CACHE_TTL:
            jwk_client = PyJWKClient(jwks_url)
            _JWKS_CACHE[jwks_url] = {"client": jwk_client, "fetched_at": now}
        else:
            jwk_client = entry["client"]

        signing_key = jwk_client.get_signing_key_from_jwt(token)
        public_key = signing_key.key

        claims = jwt_mod.decode(
            token,
            public_key,
            algorithms=["RS256"],
            audience=project_id,
            issuer=f"https://securetoken.google.com/{project_id}",
        )
        return claims
    except Exception as exc:  # pragma: no cover - environment dependent
        # Increment failure metric if available
        if _AUTH_VERIFY_FAILURES:
            try:
                _AUTH_VERIFY_FAILURES.inc()
            except Exception:
                pass

        raise RuntimeError(
            "Token verification failed. Ensure FIREBASE_PROJECT_ID is set or use FIREBASE_TEST_MODE for tests."
        ) from exc


# FastAPI dependency helper
def require_roles(*allowed_roles: str) -> Callable:
    """Return a FastAPI dependency that validates Authorization header and enforces roles.

    Usage:
        current_user = Depends(require_roles('admin', 'ts-ps-001'))
    """
    from fastapi import Request, HTTPException, status

    def _dep(request: Request):
        # If protection is disabled via env var, allow through (useful for tests)
        if os.getenv("TS_PS_001_PROTECT", "false").lower() != "true":
            return {"role": "public"}
        auth = request.headers.get("Authorization", "")
        if not auth.startswith("Bearer "):
            raise HTTPException(status_code=status.HTTP_401_UNAUTHORIZED, detail="Missing Bearer token")
        token = auth.split(" ", 1)[1]
        try:
            claims = verify_id_token(token)
        except Exception:
            raise HTTPException(status_code=status.HTTP_401_UNAUTHORIZED, detail="Invalid token")

        # Roles may be in `role` (str) or `roles` (list)
        role_val = claims.get("role") or claims.get("roles") or claims.get("roles_list")
        roles_list = []
        if isinstance(role_val, str):
            roles_list = [role_val]
        elif isinstance(role_val, (list, tuple)):
            roles_list = list(role_val)

        if not any(r in roles_list for r in allowed_roles):
            raise HTTPException(status_code=status.HTTP_403_FORBIDDEN, detail="Insufficient role")

        # Attach the claims for downstream handlers
        return claims

    return _dep


# JWKS refresher helpers -------------------------------------------------
async def start_jwks_refresher(app, interval_seconds: int = 60 * 60):
    """Start a background task on `app` that refreshes the JWKS cache periodically.

    Controlled by the `ENABLE_JWKS_REFRESH` env var. The task is stored on
    `app.state._jwks_refresher_task` so it can be cancelled on shutdown.
    """
    if os.getenv("ENABLE_JWKS_REFRESH", "false").lower() != "true":
        return

    if getattr(app.state, "_jwks_refresher_task", None):
        return

    import asyncio
    import importlib

    async def _runner():
        import logging
        logger = logging.getLogger("neurobotanica.auth.jwks_refresher")
        while True:
            try:
                project_id = os.getenv("FIREBASE_PROJECT_ID")
                if project_id:
                    jwks_url = _get_jwks_url(project_id)
                    jwt_mod = importlib.import_module("jwt")
                    PyJWKClient = getattr(jwt_mod, "PyJWKClient", None)
                    if PyJWKClient:
                        jwk_client = PyJWKClient(jwks_url)
                        _JWKS_CACHE[jwks_url] = {"client": jwk_client, "fetched_at": time.time()}
                        # log cache update
                        try:
                            logger.info("JWKS refresher updated cache for %s", jwks_url)
                        except Exception:
                            pass
            except Exception as e:
                # don't let refresher crash the loop
                try:
                    logger.exception("JWKS refresher error: %s", e)
                except Exception:
                    pass
                try:
                    if _AUTH_VERIFY_FAILURES:
                        _AUTH_VERIFY_FAILURES.inc()
                except Exception:
                    pass
            await asyncio.sleep(interval_seconds)

    task = asyncio.create_task(_runner())
    app.state._jwks_refresher_task = task


async def stop_jwks_refresher(app):
    """Stop the background JWKS refresher task if present."""
    task = getattr(app.state, "_jwks_refresher_task", None)
    if not task:
        return
    try:
        task.cancel()
        import asyncio

        try:
            await task
        except asyncio.CancelledError:
            pass
    finally:
        try:
            delattr(app.state, "_jwks_refresher_task")
        except Exception:
            pass


def get_jwks_cache_snapshot() -> Dict[str, Any]:
    """Return a serializable snapshot of the JWKS cache suitable for debugging.

    Keys: jwks_url -> {fetched_at: float, client_present: bool}
    """
    snapshot: Dict[str, Any] = {}
    for url, entry in _JWKS_CACHE.items():
        snapshot[url] = {
            "fetched_at": entry.get("fetched_at"),
            "client_present": bool(entry.get("client"))
        }
    return snapshot
