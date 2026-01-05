#!/usr/bin/env python3
"""
NeuroBotanica Server Diagnostic Script
Diagnoses server stability issues and provides detailed error reporting
"""
import sys
import traceback
import logging
import asyncio
from contextlib import asynccontextmanager

# Set up detailed logging
logging.basicConfig(
    level=logging.DEBUG,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler('server_diagnostic.log')
    ]
)
logger = logging.getLogger(__name__)

def diagnose_imports():
    """Test all imports individually to find any issues."""
    print("🔍 Testing imports...")

    try:
        from fastapi import FastAPI, Depends, HTTPException
        print("✅ FastAPI imports OK")
    except Exception as e:
        print(f"❌ FastAPI import failed: {e}")
        return False

    try:
        from fastapi.middleware.cors import CORSMiddleware
        print("✅ CORS middleware OK")
    except Exception as e:
        print(f"❌ CORS middleware failed: {e}")
        return False

    try:
        from sqlalchemy.orm import Session
        from sqlalchemy import func
        print("✅ SQLAlchemy imports OK")
    except Exception as e:
        print(f"❌ SQLAlchemy import failed: {e}")
        return False

    try:
        from backend.models.database import get_db, init_db, engine, Base
        print("✅ Database models OK")
    except Exception as e:
        print(f"❌ Database models failed: {e}")
        traceback.print_exc()
        return False

    # Test router imports individually
    routers_to_test = [
        ('backend.api.studies', 'studies'),
        ('backend.api.compounds', 'compounds'),
        ('backend.api.fda_compliance', 'fda_compliance'),
        ('backend.api.conformers', 'conformers'),
        ('backend.api.omnipath', 'omnipath'),
        ('backend.api.evidence', 'evidence'),
        ('backend.api.receptor_affinity', 'receptor_affinity'),
        ('backend.api.dimers', 'dimers'),
        ('backend.routers.patentpath', 'patentpath'),
        ('backend.routers.terpenes', 'terpenes'),
        ('backend.routers.chempath', 'chempath'),
        ('backend.routers.toxpath', 'toxpath'),
        ('backend.routers.regpath', 'regpath'),
        ('backend.routers.security', 'security'),
        ('backend.routers.genomepath', 'genomepath'),
        ('backend.routers.biopath', 'biopath'),
        ('backend.routers.clinpath', 'clinpath'),
        ('backend.routers.dispensary', 'dispensary'),
    ]

    for module_name, router_name in routers_to_test:
        try:
            module = __import__(module_name, fromlist=[router_name])
            router = getattr(module, router_name)
            print(f"✅ {router_name} router OK")
        except Exception as e:
            print(f"❌ {router_name} router failed: {e}")
            # Continue testing other routers

    try:
        from backend.middleware.token_validation import TokenValidationMiddleware
        print("✅ Token validation middleware OK")
    except Exception as e:
        print(f"❌ Token validation middleware failed: {e}")
        return False

    return True

def create_minimal_app():
    """Create a minimal FastAPI app for testing."""
    print("🏗️ Creating minimal test app...")

    from fastapi import FastAPI
    from fastapi.middleware.cors import CORSMiddleware
    from contextlib import asynccontextmanager

    @asynccontextmanager
    async def lifespan(app: FastAPI):
        print("🌿 Minimal app starting...")
        yield
        print("🛑 Minimal app shutting down...")

    app = FastAPI(
        title="NeuroBotanica Diagnostic",
        description="Minimal app for server stability testing",
        lifespan=lifespan
    )

    # Add CORS
    app.add_middleware(
        CORSMiddleware,
        allow_origins=["*"],
        allow_credentials=True,
        allow_methods=["*"],
        allow_headers=["*"],
    )

    @app.get("/")
    async def root():
        return {
            "message": "NeuroBotanica Diagnostic Server",
            "status": "running",
            "version": "diagnostic"
        }

    @app.get("/health")
    async def health():
        return {"status": "healthy"}

    return app

async def test_server_lifecycle():
    """Test the server lifecycle manually."""
    print("🔄 Testing server lifecycle...")

    from fastapi.testclient import TestClient
    app = create_minimal_app()

    # Test with TestClient
    with TestClient(app) as client:
        try:
            response = client.get("/")
            print(f"✅ Root endpoint: {response.status_code}")
            print(f"   Response: {response.json()}")

            response = client.get("/health")
            print(f"✅ Health endpoint: {response.status_code}")
            print(f"   Response: {response.json()}")

        except Exception as e:
            print(f"❌ Test client failed: {e}")
            traceback.print_exc()
            return False

    return True

def main():
    print("🚀 NeuroBotanica Server Diagnostic")
    print("=" * 50)

    # Test imports
    if not diagnose_imports():
        print("❌ Import diagnosis failed. Cannot proceed.")
        return 1

    # Test minimal app
    try:
        app = create_minimal_app()
        print(f"✅ Minimal app created with {len(app.routes)} routes")
    except Exception as e:
        print(f"❌ Minimal app creation failed: {e}")
        traceback.print_exc()
        return 1

    # Test server lifecycle
    try:
        result = asyncio.run(test_server_lifecycle())
        if not result:
            print("❌ Server lifecycle test failed")
            return 1
    except Exception as e:
        print(f"❌ Server lifecycle test exception: {e}")
        traceback.print_exc()
        return 1

    # Try to create full app
    print("🔄 Testing full NeuroBotanica app...")
    try:
        from backend.main import app as full_app
        print(f"✅ Full app imported with {len(full_app.routes)} routes")

        # Test a few routes
        from fastapi.testclient import TestClient
        with TestClient(full_app) as client:
            response = client.get("/")
            print(f"✅ Full app root endpoint: {response.status_code}")

    except Exception as e:
        print(f"❌ Full app test failed: {e}")
        traceback.print_exc()
        return 1

    print("✅ All diagnostic tests passed!")
    print("\n📋 Recommendations:")
    print("1. Server stability issue may be environment-specific")
    print("2. Try running with different uvicorn configurations")
    print("3. Check for background process conflicts")
    print("4. Consider using gunicorn with uvicorn workers")

    return 0

if __name__ == "__main__":
    sys.exit(main())