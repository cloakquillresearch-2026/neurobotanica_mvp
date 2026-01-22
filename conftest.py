"""
Pytest configuration and fixtures for NeuroBotanica tests.
This file is automatically loaded by pytest.
"""
import os
import sys
import subprocess

# Ensure database exists before running tests
def pytest_configure(config):
    """Setup test database before running tests."""
    db_path = "neurobotanica.db"

    # Only create database if it doesn't exist
    if not os.path.exists(db_path):
        print(f"\n🔧 Creating test database: {db_path}")
        try:
            # Run the setup script
            result = subprocess.run(
                [sys.executable, "scripts/setup_test_database.py", db_path],
                capture_output=True,
                text=True,
                check=True
            )
            print(result.stdout)
        except subprocess.CalledProcessError as e:
            print(f"❌ Failed to create test database: {e}")
            print(f"STDOUT: {e.stdout}")
            print(f"STDERR: {e.stderr}")
            sys.exit(1)
    else:
        print(f"\n✓ Test database already exists: {db_path}")
