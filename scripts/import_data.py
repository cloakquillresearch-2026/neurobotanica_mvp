import os
import psycopg2
from psycopg2 import sql

def import_data():
    database_url = os.getenv('DATABASE_URL')
    if not database_url:
        print("DATABASE_URL not set. Exiting.")
        return

    conn = psycopg2.connect(database_url)
    cursor = conn.cursor()

    # Read SQL file
    with open('import_data.sql', 'r') as f:
        sql_content = f.read()

    # Split into individual statements (basic split by semicolon)
    statements = [stmt.strip() for stmt in sql_content.split(';') if stmt.strip()]

    for statement in statements:
        if statement:
            try:
                cursor.execute(statement)
                print(f"Executed: {statement[:50]}...")
            except Exception as e:
                print(f"Error executing: {statement[:50]}... - {e}")

    conn.commit()
    cursor.close()
    conn.close()
    print("Data import completed.")

if __name__ == "__main__":
    import_data()