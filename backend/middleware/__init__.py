"""
Backend middleware package for NeuroBotanica.

Provides:
- Token validation middleware
- Rate limiting
- Request logging
"""
from .token_validation import (
    TokenValidationMiddleware,
    TokenDependency,
    get_current_user,
    get_admin_user,
    get_research_user,
    require_permissions,
    require_consent,
    rate_limit,
    RateLimiter,
    default_rate_limiter
)

__all__ = [
    "TokenValidationMiddleware",
    "TokenDependency",
    "get_current_user",
    "get_admin_user",
    "get_research_user",
    "require_permissions",
    "require_consent",
    "rate_limit",
    "RateLimiter",
    "default_rate_limiter"
]
