"""
Token Validation Middleware for FastAPI
Implements OmniPath token-gated access control

Supports NeuroBotanica patent claims:
- Sub-100ms token validation
- Permission-based access control
- Traditional knowledge consent checking
- Audit trail integration

Reference: NeuroBotanica MVP Development Plan - Week 2 Task 2.2
"""
from typing import Optional, Callable, List, Any
from datetime import datetime
import time
import logging
from functools import wraps

from fastapi import Request, HTTPException, status, Depends
from fastapi.security import HTTPBearer, HTTPAuthorizationCredentials
from starlette.middleware.base import BaseHTTPMiddleware
from starlette.responses import JSONResponse

from ..services.omnipath_client import (
    get_omnipath_client,
    OmniPathToken,
    OperationType
)
from ..services.provenance_tracker import get_provenance_tracker

logger = logging.getLogger(__name__)

# Security scheme for OpenAPI docs
security = HTTPBearer(auto_error=False)


class TokenValidationMiddleware(BaseHTTPMiddleware):
    """Middleware for validating OmniPath tokens on all requests.
    
    Features:
    - Sub-100ms token validation (target per MVP spec)
    - Permission checking based on endpoint
    - Audit logging for all authenticated requests
    - Public endpoint bypass
    """
    
    # Endpoints that don't require authentication
    PUBLIC_ENDPOINTS = [
        "/",
        "/health",
        "/docs",
        "/openapi.json",
        "/redoc",
        "/api/v1/health",
        "/api/v1/auth/login",
        "/api/v1/auth/register",
    ]
    
    # Endpoints that only require valid token (no specific permissions)
    AUTH_ONLY_ENDPOINTS = [
        "/api/v1/auth/refresh",
        "/api/v1/auth/logout",
        "/api/v1/user/profile",
    ]
    
    def __init__(self, app, validate_tokens: bool = True):
        """Initialize middleware.
        
        Args:
            app: FastAPI application
            validate_tokens: Whether to validate tokens (disable for testing)
        """
        super().__init__(app)
        self.validate_tokens = validate_tokens
        self.omnipath_client = get_omnipath_client()
        
        logger.info(f"TokenValidationMiddleware initialized (validate={validate_tokens})")
    
    async def dispatch(self, request: Request, call_next):
        """Process each request through token validation."""
        start_time = time.time()
        
        # Skip validation for public endpoints
        if self._is_public_endpoint(request.url.path):
            return await call_next(request)
        
        # Skip validation if disabled (testing mode)
        if not self.validate_tokens:
            return await call_next(request)
        
        # Extract token from Authorization header
        auth_header = request.headers.get("Authorization")
        
        if not auth_header:
            return JSONResponse(
                status_code=status.HTTP_401_UNAUTHORIZED,
                content={"detail": "Authorization header required"}
            )
        
        # Parse Bearer token
        if not auth_header.startswith("Bearer "):
            return JSONResponse(
                status_code=status.HTTP_401_UNAUTHORIZED,
                content={"detail": "Invalid authorization scheme. Use Bearer token."}
            )
        
        token = auth_header[7:]  # Remove "Bearer " prefix
        
        # Validate token (must be <100ms per spec)
        validation_start = time.time()
        token_data = self.omnipath_client.validate_token(token)
        validation_time = (time.time() - validation_start) * 1000
        
        if validation_time > 100:
            logger.warning(f"Token validation exceeded 100ms: {validation_time:.2f}ms")
        
        if not token_data:
            return JSONResponse(
                status_code=status.HTTP_401_UNAUTHORIZED,
                content={"detail": "Invalid or expired token"}
            )
        
        # Check if token has required permissions for this endpoint
        required_permission = self._get_required_permission(
            request.url.path,
            request.method
        )
        
        if required_permission and required_permission not in token_data.permissions:
            return JSONResponse(
                status_code=status.HTTP_403_FORBIDDEN,
                content={
                    "detail": "Insufficient permissions",
                    "required": required_permission,
                    "your_permissions": token_data.permissions
                }
            )
        
        # Attach token data to request state for downstream use
        request.state.token = token_data
        request.state.user_id = token_data.user_id
        
        # Process request
        response = await call_next(request)
        
        # Log request for audit trail
        total_time = (time.time() - start_time) * 1000
        logger.debug(
            f"Authenticated request: user={token_data.user_id}, "
            f"path={request.url.path}, method={request.method}, "
            f"time={total_time:.2f}ms"
        )
        
        return response
    
    def _is_public_endpoint(self, path: str) -> bool:
        """Check if endpoint is public (no auth required)."""
        # Exact match
        if path in self.PUBLIC_ENDPOINTS:
            return True
        
        # Prefix match for docs
        if path.startswith("/docs") or path.startswith("/redoc"):
            return True
        
        return False
    
    def _get_required_permission(self, path: str, method: str) -> Optional[str]:
        """Get required permission for endpoint based on path and method.
        
        Permission mapping follows REST conventions:
        - GET -> read:resource
        - POST -> write:resource
        - PUT/PATCH -> write:resource
        - DELETE -> delete:resource
        """
        # Auth-only endpoints don't need specific permissions
        if path in self.AUTH_ONLY_ENDPOINTS:
            return None
        
        # Extract resource from path
        parts = path.strip("/").split("/")
        
        # Find the resource name (usually after "api/v1/")
        if len(parts) >= 3 and parts[0] == "api" and parts[1].startswith("v"):
            resource = parts[2]
        elif len(parts) >= 1:
            resource = parts[0]
        else:
            resource = "unknown"
        
        # Map method to permission prefix
        method_map = {
            "GET": "read",
            "HEAD": "read",
            "POST": "write",
            "PUT": "write",
            "PATCH": "write",
            "DELETE": "delete"
        }
        
        prefix = method_map.get(method.upper(), "read")
        
        return f"{prefix}:{resource}"


class TokenDependency:
    """FastAPI dependency for token validation in route handlers.
    
    Usage:
        from backend.middleware.token_validation import get_current_user
        
        @app.get("/api/v1/protected")
        async def protected_route(user: OmniPathToken = Depends(get_current_user)):
            return {"user_id": user.user_id}
    """
    
    def __init__(
        self,
        required_permissions: Optional[List[str]] = None,
        check_consent: bool = False
    ):
        """Initialize dependency.
        
        Args:
            required_permissions: List of required permissions
            check_consent: Whether to check user consent
        """
        self.required_permissions = required_permissions or []
        self.check_consent = check_consent
        self.omnipath_client = get_omnipath_client()
    
    async def __call__(
        self,
        request: Request,
        credentials: Optional[HTTPAuthorizationCredentials] = Depends(security)
    ) -> OmniPathToken:
        """Validate token and return user info."""
        # Try to get token from request state (set by middleware)
        if hasattr(request.state, "token"):
            token_data = request.state.token
        elif credentials:
            # Fallback: validate token directly
            token_data = self.omnipath_client.validate_token(credentials.credentials)
        else:
            raise HTTPException(
                status_code=status.HTTP_401_UNAUTHORIZED,
                detail="Not authenticated"
            )
        
        if not token_data:
            raise HTTPException(
                status_code=status.HTTP_401_UNAUTHORIZED,
                detail="Invalid or expired token"
            )
        
        # Check required permissions
        for perm in self.required_permissions:
            if perm not in token_data.permissions:
                raise HTTPException(
                    status_code=status.HTTP_403_FORBIDDEN,
                    detail=f"Missing required permission: {perm}"
                )
        
        return token_data


# Pre-configured dependencies for common use cases
async def get_current_user(
    request: Request,
    credentials: Optional[HTTPAuthorizationCredentials] = Depends(security)
) -> OmniPathToken:
    """Get current authenticated user.
    
    Usage:
        @app.get("/profile")
        async def get_profile(user: OmniPathToken = Depends(get_current_user)):
            return {"user_id": user.user_id}
    """
    dep = TokenDependency()
    return await dep(request, credentials)


async def get_admin_user(
    request: Request,
    credentials: Optional[HTTPAuthorizationCredentials] = Depends(security)
) -> OmniPathToken:
    """Get current user, requiring admin permission.
    
    Usage:
        @app.get("/admin/users")
        async def list_users(user: OmniPathToken = Depends(get_admin_user)):
            # Only admins can access this
            ...
    """
    dep = TokenDependency(required_permissions=["admin:all"])
    return await dep(request, credentials)


async def get_research_user(
    request: Request,
    credentials: Optional[HTTPAuthorizationCredentials] = Depends(security)
) -> OmniPathToken:
    """Get current user, requiring research permission.
    
    Usage:
        @app.get("/api/v1/research/data")
        async def get_research_data(user: OmniPathToken = Depends(get_research_user)):
            # Only researchers can access this
            ...
    """
    dep = TokenDependency(required_permissions=["read:research_data"])
    return await dep(request, credentials)


def require_permissions(*permissions: str):
    """Decorator to require specific permissions for a route.
    
    Usage:
        @app.get("/api/v1/data")
        @require_permissions("read:data", "access:terpene_analysis")
        async def get_data():
            ...
    """
    def decorator(func: Callable):
        @wraps(func)
        async def wrapper(*args, **kwargs):
            # Extract request from args
            request = None
            for arg in args:
                if isinstance(arg, Request):
                    request = arg
                    break
            
            if not request or not hasattr(request.state, "token"):
                raise HTTPException(
                    status_code=status.HTTP_401_UNAUTHORIZED,
                    detail="Not authenticated"
                )
            
            token_data = request.state.token
            
            for perm in permissions:
                if perm not in token_data.permissions:
                    raise HTTPException(
                        status_code=status.HTTP_403_FORBIDDEN,
                        detail=f"Missing required permission: {perm}"
                    )
            
            return await func(*args, **kwargs)
        return wrapper
    return decorator


def require_consent(table_name: str, operation: OperationType):
    """Decorator to require user consent for a route.
    
    Usage:
        @app.get("/api/v1/patient/{patient_id}")
        @require_consent("patients", OperationType.READ)
        async def get_patient(patient_id: int):
            ...
    """
    def decorator(func: Callable):
        @wraps(func)
        async def wrapper(*args, **kwargs):
            # Extract request
            request = None
            for arg in args:
                if isinstance(arg, Request):
                    request = arg
                    break
            
            if not request or not hasattr(request.state, "user_id"):
                raise HTTPException(
                    status_code=status.HTTP_401_UNAUTHORIZED,
                    detail="Not authenticated"
                )
            
            tracker = get_provenance_tracker()
            allowed, consent_level = tracker.check_consent_for_access(
                user_id=request.state.user_id,
                table_name=table_name,
                operation=operation
            )
            
            if not allowed:
                raise HTTPException(
                    status_code=status.HTTP_403_FORBIDDEN,
                    detail=f"Consent required for {operation.value} operation on {table_name}"
                )
            
            return await func(*args, **kwargs)
        return wrapper
    return decorator


# Rate limiting helper (basic implementation)
class RateLimiter:
    """Simple rate limiter for API endpoints.
    
    Usage:
        limiter = RateLimiter(requests_per_minute=60)
        
        @app.get("/api/v1/data")
        async def get_data(request: Request):
            if not limiter.check(request):
                raise HTTPException(429, "Rate limit exceeded")
    """
    
    def __init__(self, requests_per_minute: int = 60):
        self.requests_per_minute = requests_per_minute
        self._request_counts: dict = {}
        self._last_reset = time.time()
    
    def check(self, request: Request) -> bool:
        """Check if request is within rate limit."""
        # Reset counts every minute
        now = time.time()
        if now - self._last_reset > 60:
            self._request_counts.clear()
            self._last_reset = now
        
        # Get client identifier (user_id or IP)
        client_id = getattr(request.state, "user_id", None)
        if not client_id:
            client_id = request.client.host if request.client else "unknown"
        
        # Increment and check
        current = self._request_counts.get(client_id, 0)
        
        if current >= self.requests_per_minute:
            return False
        
        self._request_counts[client_id] = current + 1
        return True


# Global rate limiter instance
default_rate_limiter = RateLimiter(requests_per_minute=120)


def rate_limit(requests_per_minute: int = 60):
    """Decorator for rate limiting routes.
    
    Usage:
        @app.get("/api/v1/expensive-operation")
        @rate_limit(10)  # 10 requests per minute
        async def expensive_operation():
            ...
    """
    limiter = RateLimiter(requests_per_minute)
    
    def decorator(func: Callable):
        @wraps(func)
        async def wrapper(*args, **kwargs):
            request = None
            for arg in args:
                if isinstance(arg, Request):
                    request = arg
                    break
            
            if request and not limiter.check(request):
                raise HTTPException(
                    status_code=status.HTTP_429_TOO_MANY_REQUESTS,
                    detail=f"Rate limit exceeded: {requests_per_minute} requests per minute"
                )
            
            return await func(*args, **kwargs)
        return wrapper
    return decorator
