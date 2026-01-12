"""Engine loader for TS-PS-001 and related trade-secret engines.

This module provides a small factory that will attempt to load a production
engine if present (path configurable via `INFLAMMATORY_ENGINE_PATH`), and
otherwise falls back to the test stub included in the repo.
"""
from typing import Any
import importlib
import logging
import os

logger = logging.getLogger(__name__)


def load_inflammatory_engine() -> Any:
    """Return an engine instance implementing `predict_inflammatory_synergy(...)`.

    Behavior:
    - If env `INFLAMMATORY_ENGINE_PATH` is set (e.g. "backend.services.prod_engine"),
      attempt to import that module and call `get_inflammatory_synergy_engine()`.
    - Otherwise fall back to the bundled stub `backend.services.inflammatory_synergy_engine`.
    """
    # Try explicit path from env, then default production engine module
    candidates = []
    env_path = os.getenv("INFLAMMATORY_ENGINE_PATH")
    if env_path:
        candidates.append(env_path)
    default_module = "backend.services.ts_ps_001_engine"
    if default_module not in candidates:
        candidates.append(default_module)

    for path in candidates:
        try:
            mod = importlib.import_module(path)
            factory = getattr(mod, "get_inflammatory_synergy_engine", None)
            if callable(factory):
                logger.info(f"Loading inflammatory engine from {path}")
                return factory()
            logger.warning(f"Module {path} does not expose get_inflammatory_synergy_engine()")
        except Exception as e:
            logger.exception(f"Failed to import engine module '{path}': {e}")

    # Fallback to built-in stub
    try:
        stub_mod = importlib.import_module("backend.services.inflammatory_synergy_engine")
        if hasattr(stub_mod, "get_inflammatory_synergy_engine"):
            logger.info("Falling back to bundled inflammatory_synergy_engine stub")
            return stub_mod.get_inflammatory_synergy_engine()
    except Exception as e:
        logger.exception(f"Failed to load fallback inflammatory_synergy_engine: {e}")

    raise RuntimeError("No inflammatory synergy engine available")
