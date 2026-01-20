"""
ETL Validation Script for Week 1 Checkpoint
Validates data ingestion: checks foreign keys, runs basic queries, ensures <200ms times.
"""

import os
import time
import sqlite3
import psycopg2

def get_connection():
    db_url = os.getenv('DATABASE_URL')
    if db_url:
        return psycopg2.connect(db_url)
    else:
        return sqlite3.connect('neurobotanica.db')

def validate_etl():
    conn = get_connection()
    cursor = conn.cursor()
    
    # Check record counts
    tables = ["neurobotanica_compounds", "neurobotanica_drug_interactions", "omnipath_consent_artifacts"]
    for table in tables:
        cursor.execute(f"SELECT COUNT(*) FROM {table}")
        count = cursor.fetchone()[0]
        print(f"{table}: {count} records")
    
    # Foreign key check example
    cursor.execute("""
    SELECT COUNT(*) FROM neurobotanica_compounds c
    LEFT JOIN omnipath_consent_artifacts a ON c.consent_id = a.consent_id
    WHERE c.consent_id IS NOT NULL AND a.consent_id IS NULL
    """)
    orphans = cursor.fetchone()[0]
    print(f"Orphaned consent references: {orphans}")
    
    # Performance test
    start = time.time()
    cursor.execute("SELECT * FROM neurobotanica_compounds LIMIT 100")
    results = cursor.fetchall()
    elapsed = (time.time() - start) * 1000
    print(f"Query time: {elapsed:.2f} ms")
    
    conn.close()
    return orphans == 0 and elapsed < 200
    
    # Check record counts
    tables = ["neurobotanica_compounds", "neurobotanica_drug_interactions", "omnipath_consent_artifacts"]
    for table in tables:
        cursor.execute(f"SELECT COUNT(*) FROM {table}")
        count = cursor.fetchone()[0]
        print(f"{table}: {count} records")
    
    # Foreign key check example
    cursor.execute("""
    SELECT COUNT(*) FROM neurobotanica_compounds c
    LEFT JOIN omnipath_consent_artifacts a ON c.consent_id = a.consent_id
    WHERE c.consent_id IS NOT NULL AND a.consent_id IS NULL
    """)
    orphans = cursor.fetchone()[0]
    print(f"Orphaned consent references: {orphans}")
    
    # Performance test
    start = time.time()
    cursor.execute("SELECT * FROM neurobotanica_compounds LIMIT 100")
    results = cursor.fetchall()
    elapsed = (time.time() - start) * 1000
    print(f"Query time: {elapsed:.2f} ms")
    
    conn.close()
    return orphans == 0 and elapsed < 200

if __name__ == "__main__":
    if validate_etl():
        print("ETL validation passed!")
    else:
        print("ETL validation failed.")