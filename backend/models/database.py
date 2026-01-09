"""
NeuroBotanica Database Configuration
SQLAlchemy setup for PostgreSQL/SQLite
"""
from sqlalchemy import create_engine
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker
from pydantic_settings import BaseSettings
import os


class Settings(BaseSettings):
    """Application settings loaded from environment."""
    database_url: str = "sqlite:///./neurobotanica_dev.db"
    secret_key: str = "development_secret_key_change_in_production"
    app_debug: bool = True  # Renamed from 'debug' to avoid env var conflicts
    
    class Config:
        env_file = ".env"
        extra = "ignore"  # Ignore extra env vars like DEBUG=WARN


settings = Settings()

# Create engine - SQLite for development, PostgreSQL for production
engine = create_engine(
    settings.database_url,
    connect_args={"check_same_thread": False} if "sqlite" in settings.database_url else {}
)

SessionLocal = sessionmaker(autocommit=False, autoflush=False, bind=engine)

Base = declarative_base()


def get_db():
    """Dependency for database session."""
    db = SessionLocal()
    try:
        yield db
    finally:
        db.close()


def init_db():
    """Initialize database tables.

    For local development using SQLite, drop and recreate tables to ensure
    the schema matches current models (useful during rapid iteration).
    In production (Postgres), this function intentionally only creates
    missing tables and does not drop existing data.
    """
    # If using local sqlite dev DB, recreate schema to pick up model changes
    if "sqlite" in settings.database_url and settings.database_url.startswith("sqlite:///./"):
        Base.metadata.drop_all(bind=engine)
        Base.metadata.create_all(bind=engine)
    else:
        Base.metadata.create_all(bind=engine)
