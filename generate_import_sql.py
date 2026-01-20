import json
import sqlite3

# Load the data from the dump
with open('neurobotanica_dump.json', 'r') as f:
    data = json.load(f)

# Generate SQL INSERT statements
sql_statements = []

for table_name, table_data in data.items():
    columns = table_data['columns']
    rows = table_data['rows']

    for row in rows:
        # Escape single quotes in values
        escaped_values = []
        for value in row:
            if value is None:
                escaped_values.append('NULL')
            elif isinstance(value, str):
                escaped_values.append(f"'{value.replace(chr(39), chr(39)+chr(39))}'")
            else:
                escaped_values.append(str(value))

        values_str = ', '.join(escaped_values)
        columns_str = ', '.join(f'`{col}`' for col in columns)

        insert_sql = f"INSERT INTO {table_name} ({columns_str}) VALUES ({values_str});"
        sql_statements.append(insert_sql)

# Write to a SQL file
with open('import_data.sql', 'w', encoding='utf-8') as f:
    f.write('\n'.join(sql_statements))

print(f"Generated {len(sql_statements)} INSERT statements in import_data.sql")