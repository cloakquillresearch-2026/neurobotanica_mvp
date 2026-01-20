from fastapi import FastAPI
import uvicorn
import traceback
import sys

try:
    app = FastAPI(title="Test API", version="1.0.0")

    @app.get("/")
    async def root():
        return {"message": "Test API working", "status": "ok"}

    @app.get("/health")
    async def health():
        return {"status": "healthy"}

    print("FastAPI app created successfully", flush=True)

except Exception as e:
    print(f"Error creating FastAPI app: {e}", flush=True)
    traceback.print_exc()
    sys.exit(1)

if __name__ == "__main__":
    try:
        print("Starting uvicorn server...", flush=True)
        uvicorn.run(app, host="127.0.0.1", port=8000, log_level="info")
    except Exception as e:
        print(f"Error starting server: {e}", flush=True)
        traceback.print_exc()
        sys.exit(1)