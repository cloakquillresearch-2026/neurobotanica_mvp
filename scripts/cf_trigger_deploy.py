"""Trigger a Cloudflare Pages deployment for a branch and fetch status.

Usage:
  python scripts/cf_trigger_deploy.py <branch>

Behavior:
- Reads CLOUDFLARE_API_TOKEN from environment (must be set).
- Looks for a Pages project that matches the name in frontend/wrangler.toml (or falls back to provided project name).
- Attempts to trigger a deployment for the specified branch and prints deployment info.

This script is defensive and will not print the token itself.
"""
import os
import sys
import json
import requests

BRANCH = sys.argv[1] if len(sys.argv) > 1 else None
if not BRANCH:
    print('Usage: python scripts/cf_trigger_deploy.py <branch>')
    sys.exit(2)

TOKEN = os.environ.get('CLOUDFLARE_API_TOKEN')
if not TOKEN:
    print('No CLOUDFLARE_API_TOKEN in environment; aborting')
    sys.exit(1)

HEADERS = {'Authorization': f'Bearer {TOKEN}', 'Content-Type': 'application/json'}

# Try to infer project name from wrangler.toml
PROJECT_NAME = None
try:
    with open('frontend/wrangler.toml', 'r', encoding='utf-8') as f:
        for ln in f:
            if ln.strip().startswith('name'):
                # naive parse: name = "value"
                parts = ln.split('=', 1)
                if len(parts) == 2:
                    PROJECT_NAME = parts[1].strip().strip('"').strip("'")
                    break
except Exception:
    pass

print('Using branch:', BRANCH)
print('Inferred Pages project name:', PROJECT_NAME)

# List accounts
r = requests.get('https://api.cloudflare.com/client/v4/accounts', headers=HEADERS, timeout=15)
if r.status_code != 200:
    print('Failed to list accounts', r.status_code, r.text[:500])
    sys.exit(3)

accounts = r.json().get('result', [])
if not accounts:
    print('No Cloudflare accounts found for token')
    sys.exit(4)

# For each account, list pages projects and try to find match
for acct in accounts:
    acct_id = acct.get('id')
    acct_name = acct.get('name')
    print('\nChecking account', acct_name, acct_id)
    r = requests.get(f'https://api.cloudflare.com/client/v4/accounts/{acct_id}/pages/projects', headers=HEADERS, timeout=15)
    if r.status_code != 200:
        print('  failed to list pages projects for account', acct_id, r.status_code)
        continue
    projects = r.json().get('result', [])
    for p in projects:
        name = p.get('name')
        print('  found project', name)
        if PROJECT_NAME and name != PROJECT_NAME:
            continue
        # Attempt to create a deployment for the branch
        print('  attempting to create deployment for branch', BRANCH, 'on project', name)
        body = {
            'deployment_trigger': {
                'type': 'github',
                'metadata': {
                    'branch': BRANCH
                }
            }
        }
        deploy_url = f'https://api.cloudflare.com/client/v4/accounts/{acct_id}/pages/projects/{name}/deployments'
        r2 = requests.post(deploy_url, headers=HEADERS, json=body, timeout=30)
        print('  status', r2.status_code)
        try:
            dd = r2.json()
        except Exception:
            dd = {'text': r2.text}
        print('  response keys:', list(dd.keys()))
        if r2.status_code in (200, 201) and dd.get('result'):
            res = dd['result']
            print('  deployment id:', res.get('id'))
            print('  url:', res.get('url'))
            print('  phase:', res.get('phase'))
            sys.exit(0)
        else:
            print('  create deployment failed:', dd.get('errors') or dd.get('message') or dd)

print('\nNo matching Pages project/deployment succeeded')
sys.exit(5)
