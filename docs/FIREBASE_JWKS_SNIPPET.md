# JWKS refresher

NeuroBotanica supports an optional background JWKS refresher that periodically refreshes the cached JSON Web Key Set used for Firebase JWKS verification. This is useful in staging/production to avoid hitting live JWKS endpoints on every verification and to pick up key rotations quickly.

Environment variables:
- `ENABLE_JWKS_REFRESH` : set to `true` to enable the refresher (default: `false`).
- `JWKS_REFRESH_PERIOD` : refresh interval in seconds (default: `3600`).

Example (bash / staging env):
```bash
export ENABLE_JWKS_REFRESH=true
export JWKS_REFRESH_PERIOD=3600
export FIREBASE_PROJECT_ID=neurobotanica-9ddf3
# Restart or redeploy your staging service after setting these env vars
```

Notes:
- The refresher uses the `FIREBASE_PROJECT_ID` to build the JWKS URL; ensure that secret is set in staging.
- If `PyJWT` (with `PyJWKClient`) is not installed in the runtime, the refresher will be a no-op and token verification will fall back to runtime JWKS fetch or Admin SDK if configured.
- Monitor `neurobotanica_auth_verify_failures_total` and `neurobotanica_auth_verify_attempts_total` (if Prometheus is installed) to validate refresher health.
