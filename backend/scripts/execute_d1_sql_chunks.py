"""
Split a large SQL file into smaller chunks and execute each against Cloudflare D1 via Wrangler.

Usage:
  python execute_d1_sql_chunks.py --file C:/tmp/norml_import.sql --batch-size 100

This script writes temporary chunk files to `C:/tmp` and runs `npx wrangler d1 execute` for each.
"""
import argparse
import subprocess
from pathlib import Path
import sys


def split_statements(sql_text: str):
    # Split on statement terminator followed by blank line to preserve internal JSON newlines.
    parts = sql_text.split(');\n\n')
    statements = []
    for i, p in enumerate(parts):
        if i < len(parts) - 1:
            statements.append(p + ');')
        else:
            # Last part may already end with ; or not
            if p.strip():
                statements.append(p)
    return [s.strip() for s in statements if s.strip()]


def run_chunk(file_path: Path, config_path: str):
    cmd = [
        'npx', 'wrangler', 'd1', 'execute', 'neurobotanica-clinical-evidence',
        '--file', str(file_path), '--remote', '--config', config_path
    ]
    cmd_str = ' '.join(f'"{p}"' if ' ' in p else p for p in cmd)
    print(f"Running: {cmd_str}")
    # Use shell=True so Windows correctly resolves npx from user's environment
    # Force UTF-8 decoding and replace errors to avoid decode failures in logs
    result = subprocess.run(cmd_str, capture_output=True, text=True, shell=True, encoding='utf-8', errors='replace')
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--file', required=True, help='Path to large SQL file')
    parser.add_argument('--batch-size', type=int, default=100, help='Statements per chunk')
    parser.add_argument('--tmp-dir', default='C:/tmp', help='Directory for chunk files')
    parser.add_argument('--config', default='workers/api-proxy/wrangler.toml', help='Wrangler config path')
    args = parser.parse_args()

    sql_path = Path(args.file)
    if not sql_path.exists():
        print(f"SQL file not found: {sql_path}")
        sys.exit(2)

    text = sql_path.read_text(encoding='utf-8')
    statements = split_statements(text)
    print(f"Total statements parsed: {len(statements)}")

    batch_size = args.batch_size
    tmp_dir = Path(args.tmp_dir)
    tmp_dir.mkdir(parents=True, exist_ok=True)

    for i in range(0, len(statements), batch_size):
        batch = statements[i:i+batch_size]
        chunk_file = tmp_dir / f'norml_import_chunk_{i//batch_size + 1}.sql'
        chunk_text = '\n\n'.join(batch) + '\n'
        chunk_file.write_text(chunk_text, encoding='utf-8')
        print(f"Executing chunk {i//batch_size + 1} ({len(batch)} statements) -> {chunk_file}")
        res = run_chunk(chunk_file, args.config)
        if res.returncode != 0:
            print(f"Chunk {i//batch_size + 1} failed with return code {res.returncode}")
            print(res.stdout)
            print(res.stderr)
            print("Aborting further chunks.")
            sys.exit(1)
        else:
            print(f"Chunk {i//batch_size + 1} succeeded")
            print(res.stdout)

    print("All chunks executed successfully.")


if __name__ == '__main__':
    main()
