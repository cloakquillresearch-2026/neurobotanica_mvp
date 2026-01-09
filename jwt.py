"""
Lightweight local jwt shim used only for local JWKS refresher debugging.
This file is intentionally minimal and should be removed when testing with real PyJWT.
"""

__version__ = "0.0.1-debug-shim"

class PyJWKClient:
    def __init__(self, jwks_url):
        self.jwks_url = jwks_url

    def get_signing_key_from_jwt(self, token: str):
        class _K:
            def __init__(self):
                self.key = "dummy-public-key"
        return _K()


def decode(token, key, algorithms=None, audience=None, issuer=None):
    # Minimal decode shim for local debugging only
    return {
        "aud": audience,
        "iss": issuer,
        "role": "debug",
        "sub": "debug-user",
    }
