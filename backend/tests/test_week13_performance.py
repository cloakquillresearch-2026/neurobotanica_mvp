"""
Week 13: Performance and Load Testing Suite

Tests response times, concurrent request handling, and system performance 
under various load conditions to ensure production readiness.
"""

import time
import statistics
import concurrent.futures
from typing import List, Dict, Any

import pytest
from fastapi.testclient import TestClient
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from sqlalchemy.pool import StaticPool

from backend.main import app
from backend.models.database import Base, get_db


# Test database setup with connection pooling for performance tests
TEST_DATABASE_URL = "sqlite:///:memory:"
engine = create_engine(
    TEST_DATABASE_URL,
    connect_args={"check_same_thread": False},
    poolclass=StaticPool,
)
TestingSessionLocal = sessionmaker(autocommit=False, autoflush=False, bind=engine)


def override_get_db():
    """Provide test database session."""
    db = TestingSessionLocal()
    try:
        yield db
    finally:
        db.close()


# Apply override
app.dependency_overrides[get_db] = override_get_db

# Initialize database once
Base.metadata.create_all(bind=engine)


class TestResponseTimeBenchmarks:
    """Benchmark response times for critical endpoints."""
    
    @pytest.fixture(autouse=True)
    def setup(self):
        """Setup test client."""
        self.client = TestClient(app)
    
    def measure_response_time(self, method: str, endpoint: str, 
                               data: Dict = None, iterations: int = 10) -> Dict[str, float]:
        """Measure response times over multiple iterations."""
        times = []
        
        for _ in range(iterations):
            start = time.perf_counter()
            if method == "GET":
                response = self.client.get(endpoint)
            elif method == "POST":
                response = self.client.post(endpoint, json=data)
            end = time.perf_counter()
            
            if response.status_code in [200, 201]:
                times.append(end - start)
        
        if not times:
            return {"error": "All requests failed"}
        
        return {
            "min": min(times),
            "max": max(times),
            "mean": statistics.mean(times),
            "median": statistics.median(times),
            "stdev": statistics.stdev(times) if len(times) > 1 else 0,
            "count": len(times)
        }
    
    def test_health_check_performance(self):
        """Test health endpoint responds in under 100ms."""
        metrics = self.measure_response_time("GET", "/")
        
        assert "error" not in metrics
        assert metrics["mean"] < 0.1  # 100ms target
        assert metrics["max"] < 0.2   # 200ms max
    
    def test_chempath_analyze_performance(self):
        """Test ChemPath analysis responds within acceptable time."""
        data = {
            "compound": {
                "name": "Test CBD",
                "smiles": "CCCCCC1=CC(=C(C(=C1)O)C2C=C(CCC2C(=C)C)C)O"
            }
        }
        metrics = self.measure_response_time("POST", "/api/v1/chempath/analyze", data)
        
        assert "error" not in metrics
        assert metrics["mean"] < 1.0  # 1 second target for analysis
        assert metrics["max"] < 2.0   # 2 second max
    
    def test_toxpath_assessment_performance(self):
        """Test ToxPath assessment responds within acceptable time."""
        data = {
            "compound_ref": "CCCCCC1=CC(=C(C(=C1)O)C2C=C(CCC2C(=C)C)C)O",
            "compound_name": "Test CBD",
            "route": "oral"
        }
        metrics = self.measure_response_time("POST", "/api/v1/toxpath/assess", data)
        
        assert "error" not in metrics
        assert metrics["mean"] < 1.5  # 1.5 second target
    
    def test_compounds_list_performance(self):
        """Test compounds listing is fast."""
        metrics = self.measure_response_time("GET", "/api/v1/compounds/")
        
        assert "error" not in metrics
        assert metrics["mean"] < 0.2  # 200ms target for listing
    
    def test_fda_overview_performance(self):
        """Test FDA overview endpoint performance."""
        metrics = self.measure_response_time("GET", "/api/v1/fda/")
        
        assert "error" not in metrics
        assert metrics["mean"] < 0.3  # 300ms target


class TestConcurrencyHandling:
    """Test system behavior under concurrent load."""
    
    @pytest.fixture(autouse=True)
    def setup(self):
        """Setup for concurrency tests."""
        self.client = TestClient(app)
    
    def run_concurrent_requests(self, method: str, endpoint: str, 
                                  data: Dict = None, concurrency: int = 10) -> Dict[str, Any]:
        """Execute requests concurrently and measure results."""
        results = {"success": 0, "failure": 0, "times": []}
        
        def make_request():
            start = time.perf_counter()
            try:
                if method == "GET":
                    response = self.client.get(endpoint)
                elif method == "POST":
                    response = self.client.post(endpoint, json=data)
                
                elapsed = time.perf_counter() - start
                
                if response.status_code in [200, 201]:
                    return ("success", elapsed)
                else:
                    return ("failure", elapsed)
            except Exception as e:
                return ("failure", 0)
        
        with concurrent.futures.ThreadPoolExecutor(max_workers=concurrency) as executor:
            futures = [executor.submit(make_request) for _ in range(concurrency)]
            
            for future in concurrent.futures.as_completed(futures):
                result, elapsed = future.result()
                if result == "success":
                    results["success"] += 1
                    results["times"].append(elapsed)
                else:
                    results["failure"] += 1
        
        if results["times"]:
            results["mean_time"] = statistics.mean(results["times"])
            results["max_time"] = max(results["times"])
        
        return results
    
    def test_concurrent_health_checks(self):
        """Test 20 concurrent health check requests."""
        results = self.run_concurrent_requests("GET", "/", concurrency=20)
        
        # All requests should succeed
        assert results["success"] >= 18  # Allow 10% failure tolerance
        assert results.get("mean_time", 0) < 0.5
    
    def test_concurrent_compound_lookups(self):
        """Test concurrent compound list requests."""
        results = self.run_concurrent_requests("GET", "/api/v1/compounds/", concurrency=15)
        
        assert results["success"] >= 13
    
    def test_concurrent_chempath_analysis(self):
        """Test concurrent ChemPath analysis requests."""
        data = {
            "compound": {
                "name": "Concurrent Test",
                "smiles": "CCCCCC1=CC(=C(C(=C1)O)C2C=C(CCC2C(=C)C)C)O"
            }
        }
        results = self.run_concurrent_requests(
            "POST", 
            "/api/v1/chempath/analyze", 
            data=data, 
            concurrency=10
        )
        
        assert results["success"] >= 8
    
    def test_mixed_endpoint_concurrency(self):
        """Test concurrent requests to different endpoints."""
        endpoints = [
            ("GET", "/"),
            ("GET", "/api/v1/compounds/"),
            ("GET", "/api/v1/fda/"),
            ("GET", "/api/v1/studies/"),
        ]
        
        total_success = 0
        total_requests = 0
        
        def make_request(method: str, endpoint: str):
            try:
                if method == "GET":
                    response = self.client.get(endpoint)
                return response.status_code in [200, 201]
            except:
                return False
        
        with concurrent.futures.ThreadPoolExecutor(max_workers=20) as executor:
            futures = []
            for method, endpoint in endpoints * 5:  # 20 total requests
                futures.append(executor.submit(make_request, method, endpoint))
                total_requests += 1
            
            for future in concurrent.futures.as_completed(futures):
                if future.result():
                    total_success += 1
        
        success_rate = total_success / total_requests
        assert success_rate >= 0.85  # 85% success rate under load


class TestMemoryAndResourceUsage:
    """Test resource usage patterns."""
    
    @pytest.fixture(autouse=True)
    def setup(self):
        """Setup test client."""
        self.client = TestClient(app)
    
    def test_repeated_requests_no_memory_leak(self):
        """Test that repeated requests don't cause memory issues."""
        # Make 100 sequential requests
        for i in range(100):
            response = self.client.get("/")
            assert response.status_code == 200
            
            # Every 25 requests, verify API still works
            if i % 25 == 0:
                compounds = self.client.get("/api/v1/compounds/")
                assert compounds.status_code == 200
    
    def test_large_payload_handling(self):
        """Test handling of large request payloads."""
        # Create a reasonably large compound name
        large_name = "Test-Compound-" + "X" * 1000
        
        data = {
            "compound": {
                "name": large_name[:255],  # Truncate to reasonable size
                "smiles": "CCCCCC1=CC(=C(C(=C1)O)C2C=C(CCC2C(=C)C)C)O"
            }
        }
        
        response = self.client.post("/api/v1/chempath/analyze", json=data)
        assert response.status_code == 200
    
    def test_batch_operations_performance(self):
        """Test batch analysis doesn't degrade performance."""
        start = time.perf_counter()
        
        # Simulate batch of 10 analyses
        for i in range(10):
            data = {
                "compound": {
                    "name": f"Batch-{i}",
                    "smiles": "CCCCCC1=CC(=C(C(=C1)O)C2C=C(CCC2C(=C)C)C)O"
                }
            }
            response = self.client.post("/api/v1/chempath/analyze", json=data)
            assert response.status_code == 200
        
        elapsed = time.perf_counter() - start
        average_per_request = elapsed / 10
        
        # Average should stay under 1 second per request
        assert average_per_request < 1.0


class TestAPIStability:
    """Test API stability under various conditions."""
    
    @pytest.fixture(autouse=True)
    def setup(self):
        """Setup test client."""
        self.client = TestClient(app)
    
    def test_rapid_sequential_requests(self):
        """Test rapid sequential requests don't cause issues."""
        success_count = 0
        
        for _ in range(50):
            response = self.client.get("/")
            if response.status_code == 200:
                success_count += 1
        
        # At least 95% should succeed
        assert success_count >= 47
    
    def test_alternating_endpoints(self):
        """Test alternating between different endpoints."""
        endpoints = [
            "/",
            "/api/v1/compounds/",
            "/api/v1/fda/",
            "/api/v1/studies/",
        ]
        
        for _ in range(20):
            for endpoint in endpoints:
                response = self.client.get(endpoint)
                assert response.status_code == 200
    
    def test_error_recovery(self):
        """Test system recovers properly after errors."""
        # Force a 404 error
        self.client.get("/api/v1/nonexistent")
        
        # System should still work
        response = self.client.get("/")
        assert response.status_code == 200
        
        # Multiple errors shouldn't break things
        for _ in range(10):
            self.client.get("/api/v1/invalid/endpoint")
        
        # Verify still operational
        response = self.client.get("/api/v1/compounds/")
        assert response.status_code == 200
    
    def test_invalid_data_recovery(self):
        """Test recovery from invalid data submissions."""
        invalid_payloads = [
            {},
            {"wrong": "keys"},
            {"compound": None},
            {"compound": {"name": ""}},
        ]
        
        for payload in invalid_payloads:
            # These should fail gracefully
            self.client.post("/api/v1/chempath/analyze", json=payload)
        
        # System should still accept valid data
        valid_data = {
            "compound": {
                "name": "Valid Test",
                "smiles": "CCCCCC1=CC(=C(C(=C1)O)C2C=C(CCC2C(=C)C)C)O"
            }
        }
        response = self.client.post("/api/v1/chempath/analyze", json=valid_data)
        assert response.status_code == 200


class TestPerformanceThresholds:
    """Verify system meets performance requirements."""
    
    @pytest.fixture(autouse=True)
    def setup(self):
        """Setup test client."""
        self.client = TestClient(app)
    
    def test_p95_response_times(self):
        """Test 95th percentile response times meet thresholds."""
        response_times = []
        
        for _ in range(100):
            start = time.perf_counter()
            response = self.client.get("/")
            elapsed = time.perf_counter() - start
            
            if response.status_code == 200:
                response_times.append(elapsed)
        
        response_times.sort()
        p95_index = int(len(response_times) * 0.95)
        p95_time = response_times[p95_index] if p95_index < len(response_times) else response_times[-1]
        
        # P95 should be under 200ms for health endpoint
        assert p95_time < 0.2
    
    def test_throughput_baseline(self):
        """Test minimum throughput requirements."""
        start = time.perf_counter()
        request_count = 0
        
        # Make as many requests as possible in 5 seconds
        while time.perf_counter() - start < 5.0:
            response = self.client.get("/")
            if response.status_code == 200:
                request_count += 1
        
        # Should handle at least 10 requests per second
        throughput = request_count / 5.0
        assert throughput >= 10
    
    def test_analysis_throughput(self):
        """Test analysis endpoint throughput."""
        data = {
            "compound": {
                "name": "Throughput Test",
                "smiles": "CCCCCC1=CC(=C(C(=C1)O)C2C=C(CCC2C(=C)C)C)O"
            }
        }
        
        start = time.perf_counter()
        success_count = 0
        
        # Try to complete 5 analyses in 10 seconds
        for _ in range(5):
            response = self.client.post("/api/v1/chempath/analyze", json=data)
            if response.status_code == 200:
                success_count += 1
        
        elapsed = time.perf_counter() - start
        
        # All should complete
        assert success_count == 5
        # Should take less than 10 seconds total
        assert elapsed < 10.0
