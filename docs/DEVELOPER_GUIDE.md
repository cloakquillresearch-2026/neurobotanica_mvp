Developer checks and troubleshooting

This project includes a CI smoke test to validate that critical Python modules import successfully.

Quick local checks:
- Activate venv: & ".venv/Scripts/Activate.ps1" (PowerShell) or source .venv/bin/activate
- Install deps: pip install -r requirements.txt
- Run import smoke: python scripts/check_imports.py

If any import fails:
- Re-run `pip install --force-reinstall <pkg>` for the failing package
- Remove stale/broken package dirs in `.venv/Lib/site-packages` and reinstall
- If a PyPI package appears broken, consider adding a local shim for dev/testing and open an issue to fix the package in CI.

Tracked issue: We maintain a temporary shim for `annotated_doc` to keep CI and local development working. See `docs/ISSUES/0001-remove-annotated-doc-shim.md` for the tracked removal plan and acceptance criteria. Remove the shim only after confirming the upstream package installs cleanly and tests pass without the shim.
