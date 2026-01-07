Title: Remove development shim for `annotated_doc` and rely on upstream package

Status: open
Priority: medium
Assignee: (unassigned)

Summary
-------
We added a local shim for `annotated_doc` to keep developer workflows and CI stable while the upstream PyPI distribution is problematic or incomplete. This issue tracks removing the shim and ensuring the real `annotated-doc` package can be installed and used reliably in CI and production.

Background
----------
- During dependency installs we observed broken or partial installations of `annotated-doc` (missing symbols and missing metadata) that caused import errors during test collection.
- To unblock development and tests we added a minimal shim in `annotated_doc/__init__.py` and a site shim in `.venv/Lib/site-packages/annotated_doc/__init__.py` to provide a minimal `Doc` symbol so FastAPI can import.

Goals
-----
- Remove the local shim files from the repository and the venv.
- Ensure `pip install annotated-doc` succeeds and `import annotated_doc; getattr(annotated_doc, 'Doc')` resolves as expected in CI and local dev.
- Update CI to verify real package availability and fail if shim removal breaks tests.

Acceptance criteria
-------------------
- The project builds and test suite passes in CI and locally without the shim in the repo.
- No `annotated_doc` shim file remains in the repo.
- Documentation `docs/DEVELOPER_GUIDE.md` includes notes about the removal and verification steps.

Plan
----
1. Attempt to upgrade/reinstall `annotated-doc` in a clean environment (CI job):
   - `pip install --force-reinstall "annotated-doc"`
   - Run `python -c "import annotated_doc; print(hasattr(annotated_doc, 'Doc'))"`
2. If the import succeeds and tests pass, remove the shim files:
   - Remove `annotated_doc/__init__.py` from repo
   - Remove any shim file created in `.venv/Lib/site-packages/annotated_doc/__init__.py`
3. Run full test suite and `scripts/check_imports.py` in CI.
4. If import fails upstream, open an issue with the package owner and add a note to the DEVELOPER_GUIDE with next steps and workaround instructions.

Notes
-----
- This is a maintenance task and can be scheduled for a low-risk maintenance window. The shim exists to avoid interrupting the active development and should be tracked and removed once upstream is stable.

References
----------
- Local shim file: `annotated_doc/__init__.py`
- Developer guide: `docs/DEVELOPER_GUIDE.md`
- CI smoke script: `scripts/check_imports.py`
