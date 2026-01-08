# Authentication & RBAC Plan (Firebase)

Date: 2026-01-08

## Goal
Implement Firebase-based authentication with role-based access control (RBAC) to secure backend endpoints (including TS-PS-001), enable role-aware UI, and provide audit hooks for sensitive operations.

## Scope
- Firebase Authentication (ID token verification) on backend
- Lightweight RBAC with roles: **admin**, **researcher**, **dispensary**, **consultant**, **user**
- Middleware to validate ID tokens and attach `current_user` to requests
- Decorators/guards for endpoint-level role checks
- Tests (unit & integration) and CI checks
- Minimal frontend login flow (Firebase web SDK) and role-based UI variations

## Security & Ops
- Secrets required: `FIREBASE_SERVICE_ACCOUNT_JSON` (or `FIREBASE_PROJECT_ID` + public keys fetch), `FIREBASE_CLIENT_ID` for frontend
- Audit logs: record user id, role, action, timestamp for sensitive endpoints
- Rate-limiting and basic abuse protection to be added after initial auth

## Design (high level)
1. Verify Firebase ID token on each protected request using Firebase Admin SDK or cached public keys.
2. Extract `uid` and claims; map claims to our RBAC roles (store in DB if needed).
3. Attach `User` context to request (id, roles, email, claims).
4. Provide dependency for FastAPI (`get_current_user`, `require_roles` decorator).
5. Add `auth` tests in `tests/test_auth.py` exercising token verification, role checks, and RBAC enforcement.

## Implementation Plan (first sprint)
- Add `backend/auth/firebase.py` with token verification util and cache for keys
- Add `backend/dependencies/auth.py` with FastAPI dependencies (`get_current_user`, `require_role`)
- Add integration tests that mock token verification for CI
- Add simple frontend login page and store ID token in cookie/localStorage (feature-flagged)

## Rollout
- Implement & test locally
- Open PR and run CI (unit tests + integration tests)
- Add repo secrets for e2e (when ready)

## Next actions (now)
- Create branch `feature/auth-firebase-rbac` and open PR skeleton
- Add initial files: `backend/auth/firebase.py`, `backend/dependencies/auth.py`, `tests/test_auth.py` (skeletons)

— Cloak and Quill Research