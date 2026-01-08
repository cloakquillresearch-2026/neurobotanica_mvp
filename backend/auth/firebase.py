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
