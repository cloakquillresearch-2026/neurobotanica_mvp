# annotated_doc shim removal — completion report

Date: 2026-01-08

Summary:
- The temporary local `annotated_doc` shim was removed and verified on `main`.
- Verification steps performed:
  - Ran `scripts/check_imports.py` (import-smoke): `annotated_doc` imports OK; a small set of developer-only packages (fastapi, anyio, annotated_types, redis) are not present in a minimal environment — expected and acceptable.
  - Ran full test suite: `pytest` — **248 passed**.
  - Verified Cloudflare Pages preview: https://ffd9b510.neurobotanicabudtender.pages.dev (page content validated).

Actions taken:
- Removed local repo shim and backup folder.
- Added scheduled re-check workflow `.github/workflows/annotated-doc-recheck.yml` to periodically test upstream `annotated_doc` availability.
- Added `scripts/cf_trigger_deploy.py` to trigger Pages deploys for PR/branch verification.
- Documented removal and verification in this file.

If any regressions appear, restore using the recorded backup (if available) or reintroduce a minimal shim until upstream packaging stabilizes.

— Cloak and Quill Research