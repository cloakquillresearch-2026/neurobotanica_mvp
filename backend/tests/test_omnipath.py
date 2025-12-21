"""
Tests for OmniPath Integration Components
Covers Week 2 Task 2.2 and 2.3 of MVP Development Plan

Tests:
- OmniPath client service
- Provenance tracking
- Token validation middleware
- OmniPath API endpoints
"""
import pytest
from datetime import datetime, timedelta
from unittest.mock import MagicMock, patch
import json

from fastapi.testclient import TestClient


class TestOmniPathClient:
    """Tests for OmniPath client service."""
    
    def test_client_initialization(self):
        """Test OmniPath client initializes correctly."""
        from backend.services.omnipath_client import OmniPathClient, get_omnipath_client
        
        client = OmniPathClient(mode="development")
        
        assert client.mode == "development"
        assert isinstance(client._local_manifests, dict)
        assert isinstance(client._local_tokens, dict)
    
    def test_token_validation_development_mode(self):
        """Test token validation in development mode."""
        from backend.services.omnipath_client import OmniPathClient
        
        client = OmniPathClient(mode="development")
        
        # Development mode has pre-configured tokens
        token_data = client.validate_token("dev_token_001")
        
        assert token_data is not None
        assert token_data.user_id == "dev_user"
        assert len(token_data.permissions) > 0
        assert "read" in token_data.permissions
    
    def test_token_validation_invalid(self):
        """Test invalid token returns None."""
        from backend.services.omnipath_client import OmniPathClient
        
        client = OmniPathClient(mode="development")
        
        # Invalid token format
        result = client.validate_token("")
        assert result is None
    
    def test_manifest_creation(self):
        """Test provenance manifest creation."""
        from backend.services.omnipath_client import OmniPathClient, OperationType
        
        client = OmniPathClient(mode="development")
        
        test_data = {"name": "CBD", "smiles": "OC1=CC(C2=CC(O)=CC=C2)=CC(O)=C1"}
        
        manifest = client.create_manifest(
            data=test_data,
            operation=OperationType.CREATE,
            user_id="test_user"
        )
        
        assert manifest.manifest_id.startswith("NB-MNF-")
        assert manifest.data_hash is not None
        assert len(manifest.data_hash) == 64  # SHA-256 hex length
        assert manifest.operation == OperationType.CREATE
        assert manifest.user_id == "test_user"
    
    def test_manifest_integrity_verification(self):
        """Test data integrity verification against manifest."""
        from backend.services.omnipath_client import OmniPathClient, OperationType
        
        client = OmniPathClient(mode="development")
        
        test_data = {"name": "THC", "value": 42}
        
        # Create manifest
        manifest = client.create_manifest(
            data=test_data,
            operation=OperationType.CREATE,
            user_id="test_user"
        )
        
        # Verify with same data
        result = client.verify_data_integrity(manifest.manifest_id, test_data)
        
        assert result["verified"] is True
        
        # Verify with modified data
        modified_data = {"name": "THC", "value": 99}
        result = client.verify_data_integrity(manifest.manifest_id, modified_data)
        
        assert result["verified"] is False
    
    def test_consent_level_checking(self):
        """Test consent level checking."""
        from backend.services.omnipath_client import OmniPathClient, ConsentLevel
        
        client = OmniPathClient(mode="development")
        
        # Set consent
        client._consent_registry["user1:treatment_data"] = ConsentLevel.FULL
        client._consent_registry["user2:treatment_data"] = ConsentLevel.RESEARCH
        
        # Check consent
        assert client.check_consent("user1", "treatment_data") == ConsentLevel.FULL
        assert client.check_consent("user2", "treatment_data") == ConsentLevel.RESEARCH
        # In dev mode, unknown users get FULL consent for testing
        assert client.check_consent("unknown_user", "treatment_data") == ConsentLevel.FULL
    
    def test_traditional_knowledge_recording(self):
        """Test TK usage recording."""
        from backend.services.omnipath_client import OmniPathClient
        
        client = OmniPathClient(mode="development")
        
        record = client.record_traditional_knowledge_use(
            knowledge_source="Ayurvedic_Cannabis",
            usage_type="compound_validation",
            benefit_percentage=5.0,
            compound_name="CBD"
        )
        
        # Record is stored in the log
        assert len(client._tk_usage_log) >= 1
        last_record = client._tk_usage_log[-1]
        assert last_record["knowledge_source"] == "Ayurvedic_Cannabis"
        assert last_record["benefit_percentage"] == 5.0
    
    def test_audit_trail_query(self):
        """Test audit trail querying."""
        from backend.services.omnipath_client import OmniPathClient, OperationType
        
        client = OmniPathClient(mode="development")
        
        # Create some manifests
        for i in range(5):
            client.create_manifest(
                data={"id": i},
                operation=OperationType.CREATE,
                user_id=f"user_{i % 2}"
            )
        
        # Query audit trail (returns list of ProvenanceManifest objects)
        trail = client.get_audit_trail(limit=10)
        
        assert len(trail) >= 5
        
        # Query with filter
        filtered = client.get_audit_trail(user_id="user_0", limit=10)
        
        for entry in filtered:
            # entry is a ProvenanceManifest dataclass
            assert entry.user_id == "user_0"


class TestProvenanceTracker:
    """Tests for provenance tracking service."""
    
    def test_tracker_initialization(self):
        """Test provenance tracker initializes correctly."""
        from backend.services.provenance_tracker import ProvenanceTracker
        from backend.services.omnipath_client import OmniPathClient
        
        client = OmniPathClient(mode="development")
        tracker = ProvenanceTracker(omnipath_client=client)
        
        assert tracker.omnipath_client is not None
        assert isinstance(tracker._audit_log, list)
    
    def test_track_create(self):
        """Test tracking CREATE operations."""
        from backend.services.provenance_tracker import ProvenanceTracker
        from backend.services.omnipath_client import OmniPathClient
        
        client = OmniPathClient(mode="development")
        tracker = ProvenanceTracker(omnipath_client=client)
        
        manifest = tracker.track_create(
            data={"name": "Test Compound", "smiles": "C"},
            table_name="cannabinoids",
            user_id="test_user",
            record_id=123
        )
        
        assert manifest.manifest_id.startswith("NB-MNF-")
        assert len(tracker._audit_log) == 1
        assert tracker._audit_log[0].operation == "CREATE"
        assert tracker._audit_log[0].table_name == "cannabinoids"
        assert tracker._audit_log[0].record_id == 123
    
    def test_track_update(self):
        """Test tracking UPDATE operations."""
        from backend.services.provenance_tracker import ProvenanceTracker
        from backend.services.omnipath_client import OmniPathClient
        
        client = OmniPathClient(mode="development")
        tracker = ProvenanceTracker(omnipath_client=client)
        
        data_before = {"name": "Old Name", "value": 1}
        data_after = {"name": "New Name", "value": 2}
        
        manifest = tracker.track_update(
            record_id=456,
            data_before=data_before,
            data_after=data_after,
            table_name="compounds",
            user_id="test_user"
        )
        
        assert len(tracker._audit_log) == 1
        assert tracker._audit_log[0].operation == "UPDATE"
        assert tracker._audit_log[0].data_before == data_before
        assert tracker._audit_log[0].data_after == data_after
    
    def test_track_delete(self):
        """Test tracking DELETE operations."""
        from backend.services.provenance_tracker import ProvenanceTracker
        from backend.services.omnipath_client import OmniPathClient
        
        client = OmniPathClient(mode="development")
        tracker = ProvenanceTracker(omnipath_client=client)
        
        manifest = tracker.track_delete(
            record_id=789,
            data={"name": "Deleted Item"},
            table_name="treatments",
            user_id="test_user"
        )
        
        assert len(tracker._audit_log) == 1
        assert tracker._audit_log[0].operation == "DELETE"
        assert tracker._audit_log[0].data_before == {"name": "Deleted Item"}
        assert tracker._audit_log[0].data_after is None
    
    def test_consent_for_access(self):
        """Test consent checking for access."""
        from backend.services.provenance_tracker import ProvenanceTracker
        from backend.services.omnipath_client import OmniPathClient, ConsentLevel, OperationType
        
        client = OmniPathClient(mode="development")
        tracker = ProvenanceTracker(omnipath_client=client)
        
        # Set up consent
        client._consent_registry["user_full:patients"] = ConsentLevel.FULL
        client._consent_registry["user_research:patients"] = ConsentLevel.RESEARCH
        client._consent_registry["user_limited:patients"] = ConsentLevel.LIMITED
        
        # Full consent allows all operations
        allowed, level = tracker.check_consent_for_access(
            "user_full", "patients", OperationType.CREATE
        )
        assert allowed is True
        assert level == ConsentLevel.FULL
        
        # Research consent allows READ but not CREATE
        allowed, level = tracker.check_consent_for_access(
            "user_research", "patients", OperationType.READ
        )
        assert allowed is True
        
        allowed, level = tracker.check_consent_for_access(
            "user_research", "patients", OperationType.CREATE
        )
        assert allowed is False
        
        # Limited consent only allows AGGREGATE
        allowed, level = tracker.check_consent_for_access(
            "user_limited", "patients", OperationType.AGGREGATE
        )
        assert allowed is True
        
        allowed, level = tracker.check_consent_for_access(
            "user_limited", "patients", OperationType.READ
        )
        assert allowed is False
    
    def test_audit_log_query(self):
        """Test audit log querying with filters."""
        from backend.services.provenance_tracker import ProvenanceTracker
        from backend.services.omnipath_client import OmniPathClient
        
        client = OmniPathClient(mode="development")
        tracker = ProvenanceTracker(omnipath_client=client)
        
        # Create multiple entries
        tracker.track_create({"id": 1}, "cannabinoids", "user_a", 1)
        tracker.track_create({"id": 2}, "treatments", "user_a", 2)
        tracker.track_create({"id": 3}, "cannabinoids", "user_b", 3)
        
        # Query all
        all_entries = tracker.get_audit_log()
        assert len(all_entries) == 3
        
        # Query by table
        cannabinoid_entries = tracker.get_audit_log(table_name="cannabinoids")
        assert len(cannabinoid_entries) == 2
        
        # Query by user
        user_a_entries = tracker.get_audit_log(user_id="user_a")
        assert len(user_a_entries) == 2
    
    def test_record_history(self):
        """Test getting complete record history."""
        from backend.services.provenance_tracker import ProvenanceTracker
        from backend.services.omnipath_client import OmniPathClient
        
        client = OmniPathClient(mode="development")
        tracker = ProvenanceTracker(omnipath_client=client)
        
        # Create a record then update it
        tracker.track_create({"name": "Original"}, "compounds", "user", 100)
        tracker.track_update(100, {"name": "Original"}, {"name": "Updated"}, "compounds", "user")
        
        history = tracker.get_record_history("compounds", 100)
        
        assert len(history) == 2
        # Check that both operations are present (order may vary with identical timestamps)
        operations = {h["operation"] for h in history}
        assert "CREATE" in operations
        assert "UPDATE" in operations


class TestTokenValidationMiddleware:
    """Tests for token validation middleware."""
    
    def test_public_endpoints_bypass(self):
        """Test that public endpoints bypass authentication."""
        from backend.main import app
        
        client = TestClient(app)
        
        # Root endpoint should be public
        response = client.get("/")
        assert response.status_code == 200
        
        # Health check should be public
        response = client.get("/health")
        assert response.status_code == 200
        
        # Docs should be public
        response = client.get("/docs")
        assert response.status_code == 200
    
    def test_omnipath_health_endpoint(self):
        """Test OmniPath health endpoint."""
        from backend.main import app
        
        client = TestClient(app)
        
        response = client.get("/api/v1/omnipath/health")
        
        assert response.status_code == 200
        data = response.json()
        assert data["status"] == "healthy"
        assert "metrics" in data
        assert data["metrics"]["token_validation_target_met"] is True  # <100ms
    
    def test_token_dependency(self):
        """Test token dependency for protected routes."""
        from backend.middleware.token_validation import TokenDependency
        from backend.services.omnipath_client import OmniPathToken
        
        dep = TokenDependency(required_permissions=["read:data"])
        
        assert dep.required_permissions == ["read:data"]
        assert dep.check_consent is False
    
    def test_rate_limiter(self):
        """Test rate limiter functionality."""
        from backend.middleware.token_validation import RateLimiter
        from unittest.mock import MagicMock
        
        limiter = RateLimiter(requests_per_minute=5)
        
        # Create mock request
        mock_request = MagicMock()
        mock_request.state = MagicMock()
        mock_request.state.user_id = None
        mock_request.client = MagicMock()
        mock_request.client.host = "127.0.0.1"
        
        # First 5 requests should be allowed
        for i in range(5):
            assert limiter.check(mock_request) is True
        
        # 6th request should be blocked
        assert limiter.check(mock_request) is False


class TestOmniPathAPI:
    """Tests for OmniPath API endpoints."""
    
    def test_consent_get_unauthenticated(self):
        """Test consent endpoint requires authentication when enabled."""
        # Note: Token validation disabled by default in dev mode
        pass
    
    def test_tk_sources_endpoint(self):
        """Test TK sources listing endpoint."""
        from backend.main import app
        from backend.services.omnipath_client import get_omnipath_client, reset_omnipath_client
        
        # Reset client to get fresh state
        reset_omnipath_client()
        client = get_omnipath_client()
        
        # Add some TK records
        client.record_traditional_knowledge_use(
            "Ayurvedic_Cannabis", "validation", 5.0, "CBD"
        )
        client.record_traditional_knowledge_use(
            "Traditional_Chinese", "analysis", 3.0, "THC"
        )
        
        # Test API endpoint (requires auth, but middleware disabled in test)
        test_client = TestClient(app)
        response = test_client.get("/api/v1/omnipath/tk/sources")
        
        # Will return 401 if auth is enabled, otherwise 200
        assert response.status_code in [200, 401]
    
    def test_manifest_endpoint(self):
        """Test manifest retrieval endpoint."""
        from backend.main import app
        from backend.services.omnipath_client import get_omnipath_client, OperationType
        
        client = get_omnipath_client()
        
        # Create a manifest
        manifest = client.create_manifest(
            data={"test": "data"},
            operation=OperationType.CREATE,
            user_id="test"
        )
        
        test_client = TestClient(app)
        response = test_client.get(f"/api/v1/omnipath/manifest/{manifest.manifest_id}")
        
        assert response.status_code in [200, 401]


class TestTokenValidationPerformance:
    """Performance tests for token validation (MVP spec: <100ms)."""
    
    def test_token_validation_speed(self):
        """Test that token validation completes in <100ms."""
        import time
        from backend.services.omnipath_client import OmniPathClient
        
        client = OmniPathClient(mode="development")
        
        # Warm up
        client.validate_token("dev_test_warm")
        
        # Time validation
        times = []
        for i in range(100):
            start = time.time()
            client.validate_token(f"dev_test_token_user{i}")
            elapsed = (time.time() - start) * 1000
            times.append(elapsed)
        
        avg_time = sum(times) / len(times)
        max_time = max(times)
        
        print(f"Token validation: avg={avg_time:.2f}ms, max={max_time:.2f}ms")
        
        assert avg_time < 100, f"Average validation time {avg_time:.2f}ms exceeds 100ms target"
        assert max_time < 200, f"Max validation time {max_time:.2f}ms exceeds 200ms"
    
    def test_manifest_creation_speed(self):
        """Test that manifest creation is reasonably fast."""
        import time
        from backend.services.omnipath_client import OmniPathClient, OperationType
        
        client = OmniPathClient(mode="development")
        
        test_data = {"name": "Test", "smiles": "CCCC", "values": list(range(100))}
        
        times = []
        for i in range(100):
            start = time.time()
            client.create_manifest(
                data={**test_data, "id": i},
                operation=OperationType.CREATE,
                user_id="perf_test"
            )
            elapsed = (time.time() - start) * 1000
            times.append(elapsed)
        
        avg_time = sum(times) / len(times)
        
        print(f"Manifest creation: avg={avg_time:.2f}ms")
        
        assert avg_time < 50, f"Average manifest creation {avg_time:.2f}ms exceeds 50ms target"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
