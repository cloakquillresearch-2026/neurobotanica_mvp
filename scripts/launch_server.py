#!/usr/bin/env python3
"""
NeuroBotanica Server Launcher
Programmatic server startup for Windows compatibility
"""
import uvicorn
import os
import sys
import logging

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def main():
    """Launch the NeuroBotanica FastAPI server."""
    try:
        # Set Python path for imports
        project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        sys.path.insert(0, project_root)

        logger.info("🚀 Starting NeuroBotanica API Server...")
        logger.info(f"📁 Project root: {project_root}")

        # Run the server
        uvicorn.run(
            "backend.main:app",
            host="127.0.0.1",
            port=8006,
            log_level="debug",
            access_log=True,
            reload=False  # Disable reload for production stability
        )

    except KeyboardInterrupt:
        logger.info("🛑 Server stopped by user")
    except Exception as e:
        logger.error(f"❌ Server error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()