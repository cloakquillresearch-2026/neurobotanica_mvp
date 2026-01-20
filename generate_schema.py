import json

# Load the data from the dump
with open('neurobotanica_dump.json', 'r') as f:
    data = json.load(f)

# Generate CREATE TABLE statements
create_statements = []

for table_name, table_data in data.items():
    columns = table_data['columns']
    rows = table_data['rows']

    if not rows:
        continue

    # Infer column types from first row
    column_defs = []
    for i, col in enumerate(columns):
        value = rows[0][i]
        if value is None:
            col_type = "TEXT"
        elif isinstance(value, int):
            col_type = "INTEGER"
        elif isinstance(value, float):
            col_type = "REAL"
        else:
            col_type = "TEXT"

        column_defs.append(f"    {col} {col_type}")

    create_sql = f"CREATE TABLE IF NOT EXISTS {table_name} (\n" + ",\n".join(column_defs) + "\n);"
    create_statements.append(create_sql)

# Write to schema file
with open('dynamic_schema.sql', 'w', encoding='utf-8') as f:
    f.write('\n\n'.join(create_statements))

print(f"Generated schema for {len(create_statements)} tables")