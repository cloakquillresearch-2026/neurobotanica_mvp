import os
import time
from typing import Dict, Any

_JWKS_CACHE: dict[str, Any] = {}


def _get_jwks_url(project_id: str) -> str:
    return f"https://www.googleapis.com/service_accounts/v1/jwk/securetoken@system.gserviceaccount.com"


def verify_id_token(token: str) -> Dict[str, Any]:
    """Verify a Firebase ID token.

    - If `FIREBASE_TEST_MODE=1` the token format `test-token-<role>` is accepted (useful for CI).
    - If `FIREBASE_SERVICE_ACCOUNT_JSON` is present the code will attempt to use firebase_admin (legacy path).
    - Otherwise it attempts JWKS verification against Google's public keys for Firebase Auth.
    """
    # Test mode: lightweight deterministic verifier for unit tests
    if os.getenv("FIREBASE_TEST_MODE") == "1":
        if token.startswith("test-token-"):
            role = token.split("test-token-")[-1]
            return {"uid": f"{role}-uid", "email": f"{role}@example.com", "role": role}
        raise ValueError("Invalid test token")

    # If firebase_admin is available and service account provided, prefer it
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

        claims = fb_auth.verify_id_token(token)
        return claims
    except Exception:
        # fall through to JWKS approach
        pass

    # JWKS verification path (no firebase_admin)
    project_id = os.getenv("FIREBASE_PROJECT_ID")
    if not project_id:
        raise RuntimeError("FIREBASE_PROJECT_ID is required for JWKS-based token verification")

    try:
        # Lazy imports so tests which use FIREBASE_TEST_MODE do not require PyJWT/httpx
        import httpx
        import jwt
        from jwt import PyJWKClient

        jwks_url = _get_jwks_url(project_id)
        jwk_client = PyJWKClient(jwks_url)
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
        raise RuntimeError(
            "Token verification failed. Ensure FIREBASE_PROJECT_ID is set or use FIREBASE_TEST_MODE for tests."
        ) from exc
import os
from typing import Dict


def verify_id_token(token: str) -> Dict:
    """Verify a Firebase ID token.

    In CI and local unit tests set `FIREBASE_TEST_MODE=1` and use test tokens
    like `test-token-admin`, `test-token-user`, etc. In production this will
    attempt to use the Firebase Admin SDK; if it's not configured a RuntimeError
    is raised.
    """
    # Test mode: lightweight deterministic verifier for unit tests
    if os.getenv("FIREBASE_TEST_MODE") == "1":
        if token.startswith("test-token-"):
            role = token.split("test-token-")[-1]
            return {"uid": f"{role}-uid", "email": f"{role}@example.com", "role": role}
        raise ValueError("Invalid test token")

    # Production mode: try firebase_admin if available
    try:
        import firebase_admin  # type: ignore
        from firebase_admin import auth as fb_auth  # type: ignore

        # Initialize the admin app lazily if not already
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

    except Exception as exc:  # pragma: no cover - platform/environment dependent
        raise RuntimeError(
            "Firebase not configured or verification failed. "
            "For local tests set FIREBASE_TEST_MODE=1 and use test tokens."
        ) from exc
