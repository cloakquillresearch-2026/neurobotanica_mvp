#!/usr/bin/env python3
"""
NeuroBotanica API Test Script
Tests the FastAPI server functionality programmatically
"""
import requests
import time
import threading
import uvicorn
import os
import sys
import logging

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def run_server():
    """Run the server in a background thread."""
    try:
        # Set Python path for imports
        project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        sys.path.insert(0, project_root)

        uvicorn.run(
            "backend.main:app",
            host="127.0.0.1",
            port=8007,
            log_level="info",
            access_log=False,  # Disable access logs for cleaner output
            reload=False
        )
    except Exception as e:
        logger.error(f"Server error: {e}")

def test_api():
    """Test the API endpoints."""
    base_url = "http://127.0.0.1:8007"

    # Wait for server to start
    time.sleep(3)

    tests = [
        ("GET", "/", "Root endpoint"),
        ("GET", "/test", "Test endpoint"),
        ("GET", "/docs", "API documentation"),
        ("GET", "/openapi.json", "OpenAPI schema"),
    ]

    results = []

    for method, endpoint, description in tests:
        try:
            if method == "GET":
                response = requests.get(f"{base_url}{endpoint}", timeout=5)
            else:
                continue

            results.append({
                "endpoint": endpoint,
                "description": description,
                "status": response.status_code,
                "success": response.status_code == 200
            })

            logger.info(f"✅ {description}: {response.status_code}")

        except requests.exceptions.RequestException as e:
            results.append({
                "endpoint": endpoint,
                "description": description,
                "status": "ERROR",
                "success": False,
                "error": str(e)
            })
            logger.error(f"❌ {description}: {e}")

    return results

def main():
    """Main test function."""
    logger.info("🧪 Starting NeuroBotanica API Tests...")

    # Start server in background thread
    server_thread = threading.Thread(target=run_server, daemon=True)
    server_thread.start()

    try:
        # Run tests
        results = test_api()

        # Print summary
        logger.info("\n📊 Test Results:")
        successful = sum(1 for r in results if r["success"])
        total = len(results)

        for result in results:
            status = "✅" if result["success"] else "❌"
            logger.info(f"{status} {result['description']}: {result['status']}")

        logger.info(f"\n🎯 Summary: {successful}/{total} tests passed")

        if successful == total:
            logger.info("🎉 All API endpoints are working correctly!")
            return 0
        else:
            logger.error("⚠️ Some tests failed")
            return 1

    except KeyboardInterrupt:
        logger.info("🛑 Tests interrupted by user")
        return 1
    except Exception as e:
        logger.error(f"❌ Test error: {e}")
        import traceback
        traceback.print_exc()
        return 1

if __name__ == "__main__":
    sys.exit(main())