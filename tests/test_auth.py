import os
import pytest

from backend.auth.firebase import verify_id_token
from backend.dependencies.auth import _extract_roles_from_claims


def test_verify_id_token_testmode():
    os.environ["FIREBASE_TEST_MODE"] = "1"

    claims = verify_id_token("test-token-admin")
    assert claims["uid"] == "admin-uid"
    assert claims["email"] == "admin@example.com"
    assert claims["role"] == "admin"


def test_extract_roles_from_claims():
    assert _extract_roles_from_claims({"role": "admin"}) == ["admin"]
    assert _extract_roles_from_claims({"roles": ["admin", "researcher"]}) == ["admin", "researcher"]
    assert _extract_roles_from_claims({}) == []
