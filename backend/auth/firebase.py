import os
import time
from typing import Dict, Any, Callable

_JWKS_CACHE: dict[str, Any] = {}
_JWKS_CACHE_TTL = 60 * 60  # seconds

# Optional Prometheus metrics (lazy)
_AUTH_VERIFY_ATTEMPTS = None
_AUTH_VERIFY_FAILURES = None
_AUTH_JWKS_LAST_REFRESH = None
_AUTH_JWKS_REFRESH_DURATION = None
try:
    from prometheus_client import Counter, Gauge, Histogram

    _AUTH_VERIFY_ATTEMPTS = Counter("neurobotanica_auth_verify_attempts_total", "Auth verify attempts")
    _AUTH_VERIFY_FAILURES = Counter("neurobotanica_auth_verify_failures_total", "Auth verify failures")
    _AUTH_JWKS_LAST_REFRESH = Gauge("neurobotanica_auth_jwks_last_refresh_timestamp", "Last JWKS refresh timestamp")
    _AUTH_JWKS_REFRESH_DURATION = Histogram("neurobotanica_auth_jwks_refresh_seconds", "JWKS refresh duration seconds")
except Exception:
    _AUTH_VERIFY_ATTEMPTS = None
    _AUTH_VERIFY_FAILURES = None
    _AUTH_JWKS_LAST_REFRESH = None
    _AUTH_JWKS_REFRESH_DURATION = None


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

    # Try firebase_admin if available
    try:
        import firebase_admin  # type: ignore
        from firebase_admin import auth as fb_auth  # type: ignore

        if not firebase_admin._apps:
            cred_json = os.getenv("FIREBASE_SERVICE_ACCOUNT_JSON")
            if cred_json:
                import json
                from firebase_admin import credentials

                cred = credentials.Certificate(json.loads(cred_json))
                firebase_admin.initialize_app(cred)
            else:
                raise RuntimeError("FIREBASE_SERVICE_ACCOUNT_JSON not set; cannot initialize firebase_admin")

        claims = fb_auth.verify_id_token(token)
        return claims
    except Exception:
        # fallthrough to JWKS
        pass

    # JWKS verification path
    project_id = os.getenv("FIREBASE_PROJECT_ID")
    if not project_id:
        raise RuntimeError("FIREBASE_PROJECT_ID is required for JWKS-based token verification")

    try:
        # Lazy imports to avoid requiring PyJWT for tests using FIREBASE_TEST_MODE
        import httpx
        import jwt
        from jwt import PyJWKClient

        # Increment attempt metric if available
        if _AUTH_VERIFY_ATTEMPTS:
            try:
                _AUTH_VERIFY_ATTEMPTS.inc()
            except Exception:
                pass

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

        claims = jwt.decode(
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


    # JWKS refresher background task and lifecycle helpers
    _JWKS_REFRESH_TASK = None


    async def _jwks_refresher_loop(period: int = 3600):
        import asyncio
        import logging
        logger = logging.getLogger(__name__)

        project_id = os.getenv("FIREBASE_PROJECT_ID")
        if not project_id:
            logger.info("JWKS refresher: FIREBASE_PROJECT_ID not set, skipping refresher.")
            return

        jwks_url = _get_jwks_url(project_id)
        while True:
            try:
                start = time.time()
                # Lazy import
                from jwt import PyJWKClient

                # Fetch and replace client
                jwk_client = PyJWKClient(jwks_url)
                _JWKS_CACHE[jwks_url] = {"client": jwk_client, "fetched_at": time.time()}

                duration = time.time() - start
                if _AUTH_JWKS_LAST_REFRESH:
                    try:
                        _AUTH_JWKS_LAST_REFRESH.set(time.time())
                    except Exception:
                        pass
                if _AUTH_JWKS_REFRESH_DURATION:
                    try:
                        _AUTH_JWKS_REFRESH_DURATION.observe(duration)
                    except Exception:
                        pass

                logger.debug("JWKS refreshed successfully in %.3fs", duration)
            except Exception as e:
                logger.warning("JWKS refresh failed: %s", e)
            finally:
                await asyncio.sleep(period)


    async def start_jwks_refresher(app, period: int | None = None):
        """Start JWKS refresher background task and attach to `app.state`."""
        import asyncio
        import logging
        logger = logging.getLogger(__name__)

        if period is None:
            try:
                period = int(os.getenv("FIREBASE_JWKS_REFRESH_SECONDS", "3600"))
            except Exception:
                period = 3600

        global _JWKS_REFRESH_TASK
        if _JWKS_REFRESH_TASK is not None:
            return

        _JWKS_REFRESH_TASK = asyncio.create_task(_jwks_refresher_loop(period))
        try:
            app.state._jwks_refresh_task = _JWKS_REFRESH_TASK
        except Exception:
            logger.debug("Could not attach JWKS task to app.state")


    async def stop_jwks_refresher(app):
        import asyncio
        import logging
        logger = logging.getLogger(__name__)

        global _JWKS_REFRESH_TASK
        task = getattr(app.state, "_jwks_refresh_task", None) or _JWKS_REFRESH_TASK
        if task:
            task.cancel()
            try:
                await task
            except asyncio.CancelledError:
                logger.debug("JWKS refresher task cancelled")
        _JWKS_REFRESH_TASK = None
