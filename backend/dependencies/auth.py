from typing import List
from pydantic import BaseModel
from fastapi import Depends, HTTPException, status
from fastapi.security import HTTPBearer, HTTPAuthorizationCredentials

from backend.auth.firebase import verify_id_token

security = HTTPBearer()


class User(BaseModel):
    uid: str
    email: str | None = None
    roles: List[str] = []


def _extract_roles_from_claims(claims: dict) -> List[str]:
    r = []
    maybe = claims.get("role") or claims.get("roles") or []
    if isinstance(maybe, str):
        r = [maybe]
    elif isinstance(maybe, (list, tuple)):
        r = list(maybe)
    return r


async def get_current_user(credentials: HTTPAuthorizationCredentials = Depends(security)) -> User:
    token = credentials.credentials
    try:
        claims = verify_id_token(token)
    except Exception as exc:
        raise HTTPException(status_code=status.HTTP_401_UNAUTHORIZED, detail="Invalid or missing auth token")

    roles = _extract_roles_from_claims(claims)
    return User(uid=claims.get("uid"), email=claims.get("email"), roles=roles)


def require_roles(*required: str):
    async def _require(user: User = Depends(get_current_user)):
        if not any(r in user.roles for r in required):
            raise HTTPException(status_code=status.HTTP_403_FORBIDDEN, detail="Insufficient role")
        return user

    return _require
