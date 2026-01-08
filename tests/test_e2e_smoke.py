import os
import requests
import pytest

PREVIEW_URL = os.environ.get("PAGES_PREVIEW_URL") or os.environ.get("PREVIEW_URL")

pytestmark = pytest.mark.skipif(not PREVIEW_URL, reason="No preview URL available in env; skipping E2E smoke tests")


def test_preview_root_returns_200_and_contains_brand():
    assert PREVIEW_URL is not None
    r = requests.get(PREVIEW_URL, timeout=10)
    assert r.status_code == 200
    assert "NeuroBotanica" in r.text


def test_preview_has_consultation_ui():
    r = requests.get(PREVIEW_URL, timeout=10)
    assert "Customer Consultation" in r.text or "Consultation" in r.text
