# Budtender Application — Complete Workflow Specification

## Purpose

This document provides **complete, unambiguous specifications** for the NeuroBotanica Budtender Education Assistant. Any coding agent working on this application must implement these workflows exactly as described.

> **For Coding Agents:** Do not claim features are "working" unless you have verified them against this specification. Do not fabricate test results. If something doesn't work, report it honestly.
> 

---

## 🚨 CURRENT SITUATION (Feb 2, 2026)

<aside>
🟢

**STATUS (Feb 2, 2026):** All critical backend fixes deployed. Worker version `ecfb5b65-2986-4990-8e1b-f3af9fa8ba88` is live.

</aside>

### ✅ WORKING (Verified Feb 2, 2026)

| Feature | Evidence |
| --- | --- |
| **Profile saving** | Created "Benjamin Workman" → `profile_1770057470543_8f78j848` |
| **Customer search** | Search for "Benjamin" returns the saved customer |
| **Analysis endpoint** | Returns `polysaccharide_effects`, `bias_correction`, etc. |
| **Inflammatory synergy** | Returns 200 OK with `confidence: 0.766` |

### 🔴 STILL BROKEN (Verified Feb 2, 2026)

| Issue | Error | Fix Needed |
| --- | --- | --- |
| **Transaction table schema wrong** | `D1_ERROR: table dispensary_transactions has no column named customer_id: SQLITE_ERROR` | The table was created but with wrong columns. Need to run migration to add `customer_id` column. |
| **Synergy score defaults to 0.5** | `synergy_score: 0.5, evidence: "Query failed, using default prediction"` | The synergy query is failing internally. Need to check what query is being attempted and fix it. |
| **"Check your connection" errors** | User sees connection errors on every button click | API calls are failing — likely due to transaction endpoint 500 errors cascading or CORS issues on certain routes. |

---

## 🚨 NEXT STEPS FOR CODING AGENT (Feb 2, 2026)

<aside>
🎯

**AGENT: READ THIS FIRST.** Only 4 things are working. Everything else below is STILL BROKEN and needs your attention. Do not claim anything works until you have tested it and the human has verified.

</aside>

### Priority Order of Fixes

**P0 — Fix these FIRST (blocking all other functionality):**

1. **Transaction table schema** — Table exists but is missing `customer_id` column. This causes 500 errors that cascade into "check your connection" errors throughout the app.
2. **CORS on error responses** — When the Worker returns a 500, it doesn't include CORS headers, so the browser blocks the response and the frontend shows "check your connection" instead of the actual error.

**P1 — Fix after P0:**

1. **Synergy query failure** — The `/api/dispensary/analyze` endpoint returns `synergy_score: 0.5` with `"Query failed, using default prediction"`. The underlying query is failing silently. Debug and fix.
2. **New Client registration flow** — The UI still shows conditions and analysis widgets during registration. It should show ONLY the name form until the client is saved.

**P2 — Fix after P1:**

1. **"undefined" in recommendation titles** — Product name not being passed correctly
2. **THC/CBD shows 0%/0%** — Should show actual values or "N/A"
3. **Notes not persisting** — Notes field doesn't save to database

---

### What You Must Do Right Now

**Step 1: Fix the transaction table**

```sql
-- Run this in D1 console or via wrangler d1 execute:
DROP TABLE IF EXISTS dispensary_transactions;

CREATE TABLE dispensary_transactions (
  id TEXT PRIMARY KEY,
  customer_id TEXT NOT NULL,
  profile_id TEXT,
  consultation_date TEXT NOT NULL,
  conditions TEXT,
  recommendations TEXT,
  notes TEXT,
  created_at TEXT DEFAULT (datetime('now')),
  FOREIGN KEY (customer_id) REFERENCES dispensary_profiles(customer_id)
);
```

**Step 2: Add CORS to ALL error responses**

In `workers/api-proxy/src/index.ts`, find every `catch` block and ensure it returns CORS headers:

```tsx
const corsHeaders = {
  'Access-Control-Allow-Origin': '*',
  'Access-Control-Allow-Methods': 'GET, POST, PUT, DELETE, OPTIONS',
  'Access-Control-Allow-Headers': 'Content-Type, Authorization',
};

// EVERY catch block must do this:
catch (error) {
  return new Response(
    JSON.stringify({ error: error.message }),
    { status: 500, headers: { ...corsHeaders, 'Content-Type': 'application/json' } }
  );
}
```

**Step 3: Debug the synergy query**

Find the code that returns `"Query failed, using default prediction"` and add logging to see what query is failing. Likely a missing table.

**Step 4: Deploy the Worker**

```bash
cd workers/api-proxy
npx wrangler deploy
```

**Step 5: Report back with commit hash**

Do NOT say "it's live" or "changes deployed." Only say:

> ✅ Code pushed to GitHub
> 

> - Commit: `[hash]`
> 

> - Files changed: `[list]`
> 

> ⏳ Awaiting manual deployment
> 

---

## Detailed Fix Instructions

### Fix 1: Transaction Table Schema

**Problem:** The `dispensary_transactions` table exists but is missing required columns.

**Error:**

```
D1_ERROR: table dispensary_transactions has no column named customer_id: SQLITE_ERROR
```

**Required columns for `dispensary_transactions`:**

```sql
CREATE TABLE IF NOT EXISTS dispensary_transactions (
  id TEXT PRIMARY KEY,
  customer_id TEXT NOT NULL,
  profile_id TEXT,
  consultation_date TEXT NOT NULL,
  conditions TEXT,
  recommendations TEXT,
  notes TEXT,
  created_at TEXT DEFAULT (datetime('now')),
  FOREIGN KEY (customer_id) REFERENCES dispensary_profiles(customer_id)
);
```

**To fix:**

1. Check current schema: `SELECT sql FROM sqlite_master WHERE name='dispensary_transactions';`
2. Either DROP and recreate, or ALTER TABLE to add missing columns
3. Redeploy the Worker after migration

**File:** `workers/api-proxy/src/index.ts` — find the transaction INSERT statement and verify column names match

---

### Fix 2: Synergy Query Failure

**Problem:** The analyze endpoint returns `synergy_score: 0.5` with evidence "Query failed, using default prediction"

**This means:**

- The synergy calculation code has a try/catch that falls back to 0.5
- The actual query inside is throwing an error
- Likely a missing table, column, or malformed SQL

**To debug:**

1. Find the synergy calculation in `workers/api-proxy/src/index.ts`
2. Look for the try/catch block that returns `"Query failed, using default prediction"`
3. Add logging to see what query is being attempted
4. Check if the required tables exist (e.g., `terpene_interactions`, `compound_synergies`, etc.)

**Expected behavior:** Synergy should vary based on:

- Selected conditions
- Compound combinations (CBD+THC should score higher than baseline)
- Biomarker inputs

---

### Fix 3: "Check Your Connection" Errors

**Problem:** User sees connection errors on almost every button click.

**Likely causes:**

1. **500 errors from transaction endpoint** — When the transaction save fails with D1_ERROR, the frontend may show a generic connection error
2. **CORS headers missing on error responses** — When the Worker returns a 500, it may not include CORS headers, causing the browser to block the response
3. **Unhandled exceptions** — If the Worker throws before setting headers, CORS fails

**To fix:**

```tsx
// In workers/api-proxy/src/index.ts
// Wrap ALL responses in CORS headers, including errors:

const corsHeaders = {
  'Access-Control-Allow-Origin': '*',
  'Access-Control-Allow-Methods': 'GET, POST, PUT, DELETE, OPTIONS',
  'Access-Control-Allow-Headers': 'Content-Type, Authorization',
};

// For error responses:
catch (error) {
  return new Response(
    JSON.stringify({ error: error.message }),
    { 
      status: 500, 
      headers: { ...corsHeaders, 'Content-Type': 'application/json' } 
    }
  );
}
```

**Verify:** After the fix, error responses should still include CORS headers so the frontend can read the error message instead of showing "check your connection."

---

### Deployment After Fixes

```bash
# After pushing fixes to GitHub:
cd workers/api-proxy
npx wrangler deploy

# If frontend changes were made:
cd ../frontend
npm run build
npx wrangler pages deploy ./out --project-name=neurobotanicabudtender
```

### 🧪 END-TO-END VALIDATION NEEDED

Please test the full flow:

1. **Hard refresh** the app (Cmd+Shift+R / Ctrl+Shift+R)
2. **Create a new client** — Enter first name, last name → Save
3. **Search for that client** — Should now appear in results
4. **Select conditions** → Click "Run Analysis"
5. **Complete the consultation** (generates recommendations)
6. **Check the History tab** — Should show the consultation record

### 🔍 If Issues Persist

Open DevTools (F12) → Console tab and look for:

- `D1_TYPE_ERROR` — means safe extractor isn't working
- `404` on `/api/dispensary/transaction` — means endpoint not deployed
- CORS errors — check which origin is being blocked

Screenshot any errors and share here.

### Recent Commits

| Commit | Description |
| --- | --- |
| `5ee8e53` | "fix: Handle profile id response fields" — adds optional chaining to fix TypeError crash |
| `9e45fb6` | "chore: Clean up recommendation retry handler" — minor cleanup |
| `7262262` | "fix: Improve new client registration flow" — name fields to top, hide conditions, stop auto-compute, add save logs |
| `63e9731` | Debug logging added to CustomerProfile.tsx save handler |
| `725d66c` | "fix: Address frontend issues per spec" |
| `b13778b` | "fix: Improve dynamic calculations for synergy and microbiome scores" |
| `a7f7094` | "fix: Implement correct '+ New Client' workflow with modal" |

### What's Working Now (Verified Feb 1, 2026 7:50 PM)

- ✅ Name fields are now at the **TOP** of the form (before conditions)
- ✅ TS-PS-001 no longer auto-runs on page load (no spinning wheel!)
- ✅ Console shows `Save Client clicked` and `Create profile response` — **API is being called**
- ✅ Condition buttons no longer open registration forms

### What's STILL Wrong (Verified Feb 1, 2026 7:50 PM)

- ❌ **"+ New Client" button still visible in Practice Mode** — should be hidden
- ❌ **Primary Conditions still visible** during registration (should be hidden until after client saved)
- ❌ **API call crashes:** `TypeError: Cannot read properties of undefined (reading 'startsWith')`
- ❌ **customer_id is undefined** — response from API not being parsed correctly
- ❌ **"undefined Evidence-based Recommendation"** still showing

### What Agent Claimed vs Reality (Commit `725d66c`)

| Agent Claimed | Reality (Verified by Human) |
| --- | --- |
| "✅ Hide + New Client in Practice Mode" | ❌ Button still visible in Practice Mode |
| "✅ Stop Auto-Computing in Practice Mode" | ❌ "Computing..." still runs |
| "✅ Name Fields Already at Top" | ❌ Name fields are at BOTTOM |
| "✅ Conditions Already Hidden" | ❌ Conditions visible in registration modal |
| "✅ API Call Already Implemented" | ❌ Network tab shows NO API call on Save |

### What's STILL Wrong (Verified Feb 1, 2026 8:30 PM)

- ❌ **Name fields are at the BOTTOM** — Should be at TOP, before conditions
- ❌ **Conditions visible during registration** — Nausea, Depression, etc. should be HIDDEN until after client saved
- ❌ **"Computing..." auto-runs** — TS-PS-001 starts computing on page load, should wait for "Run Analysis"
- ❌ **Practice Mode has "+ New Client" button** — Button should be hidden in Practice/Sandbox mode
- ❌ **No API call on Save Client** — Network tab shows only CSS files when form submitted
- ❌ **Console errors** — 404 for `/icon-192.png`, TypeError on undefined

### What's Still Broken (Priority Order)

1. **🔴 P0: New clients not saved to database** — API Worker deployed ✅, but frontend still only updates local state (no API call visible in Network tab)
2. **🔴 P1: Customer search not working** — Cannot find previously created customers
3. **🔴 P1: Synergy Score hardcoded** — Fix in `b13778b` deployed to Pages but **Worker not deployed** — still shows 50.0%
4. **🔴 P1: Microbiome Modulation hardcoded** — Fix in `b13778b` deployed to Pages but **Worker not deployed** — still shows ~75%
5. **🟡 P2: Notes not persisting** — Notes field doesn't save to database
6. **🟡 P2: "undefined" in recommendations** — Missing product name

### Your Next Task

**Fix the database save for New Client Registration.**

---

## 🔴 ROOT CAUSE CONFIRMED (Feb 1, 2026 7:00 PM)

### Problem: Save Client button does NOT call the API

**Evidence from Network tab:**

- Only `preflight` requests to `/profile` are visible
- NO `POST` request with `first_name`/`last_name` payload
- The modal closes, but data only goes to local React state

### What's Happening

1. User clicks "+ New Client" → Modal opens ✅
2. User enters name → Input fields work ✅
3. User clicks "Save Client" → Modal closes ✅
4. **API call to `POST /api/dispensary/profile`** → ❌ **NEVER HAPPENS**
5. Name is stored in React state only → Lost on page refresh

### Likely Cause

The "Save Client" button in the New Client modal is probably calling a local state update function instead of the actual API save handler. The code needs to:

```tsx
// WRONG (current behavior)
const handleSave = () => {
  setCustomer({ ...customer, first_name, last_name })  // Local state only
  setShowModal(false)
}

// CORRECT (needed behavior)  
const handleSave = async () => {
  const profileData = { first_name, last_name, /* other fields */ }
  const result = await dispensaryAPI.createProfile(profileData)  // API CALL
  setCustomer({ ...customer, customer_id: result.profile_id })
  setShowModal(false)
}
```

### Files to Check

1. **Where is the New Client modal?**
    - `frontend/src/pages/index.tsx` — likely contains `showNewClientModal` state
    - `frontend/src/components/CustomerSearch.tsx` — may have the modal
    - `frontend/src/components/NewClientModal.tsx` — if this exists
2. **Find the Save button's onClick handler**
    - Look for `onClick` on the "Save Client" button
    - Trace what function it calls
    - Verify it calls `dispensaryAPI.createProfile()`
3. **Check `CustomerProfile.tsx`**
    - The `handleSave` function exists and calls the API
    - But is this component even used in the New Client modal?
    - Or is there a DIFFERENT save handler for new clients?

### Next Step for Coding Agent

**Find the New Client modal code and wire the Save button to the API.**

Don't assume `CustomerProfile.tsx`'s `handleSave` is being called. The New Client modal may have its own separate (broken) save handler that needs to be fixed or replaced.

---

## 🔧 EXACT FIX NEEDED (Feb 1, 2026)

### The Good News

The API call code **already exists** — it just isn't being called from the New Client modal.

**Working code in `api.ts`:**

```tsx
createProfile: (profile) => budtenderApi.post('/api/dispensary/profile', profile)
```

**Working save handler in `CustomerProfile.tsx`:**

```tsx
savedProfile = await dispensaryAPI.createProfile(profileData)
```

### What You Need to Do

1. **Open** `frontend/src/components/CustomerSearch.tsx` in VS Code
2. **Find** the Save Client button's `onClick` handler. It probably looks like this:

```tsx
// WRONG (current behavior - only updates local state)
onClick={() => {
  setCustomer({ ...customer, first_name, last_name })
  setShowModal(false)
}}
```

1. **Change it** to call the API:

```tsx
// CORRECT (calls API to persist to database)
onClick={async () => {
  const result = await dispensaryAPI.createProfile({
    first_name,
    last_name,
    // other fields from the form
  })
  setCustomer({ ...customer, customer_id: result.data.profile_id })
  setShowModal(false)
}}
```

1. **Make sure** `dispensaryAPI` is imported at the top of the file:

```tsx
import { dispensaryAPI } from '@/utils/api'
```

1. **Deploy** with Wrangler after pushing the commit:

```bash
cd frontend && npm run build && npx wrangler pages deploy ./out --project-name=neurobotanicabudtender
```

### Summary

**One file, one change.** The infrastructure is already there — just wire the button to the existing API call.

---

### ✅ API Architecture Fixed (Feb 1, 2026)

**The api-proxy Worker is now deployed with routes to [`budtender.neuro-botanica.com/api/*`](http://budtender.neuro-botanica.com/api/*).**

#### The Problem

The frontend in `frontend/src/utils/api.ts` is configured to call:

```jsx
const BUDTENDER_API_URL = 
  process.env.NEXT_PUBLIC_BUDTENDER_API_URL ||
  'https://budtender.neuro-botanica.com'
```

But [`budtender.neuro-botanica.com`](http://budtender.neuro-botanica.com) is the **Cloudflare Pages frontend** — NOT an API endpoint!

#### Where the APIs Actually Live

| File | Deployed URL | Synergy Value |
| --- | --- | --- |
| `workers/terpene-api/src/index.ts` | [`terpene-api.contessapetrini.workers.dev`](http://terpene-api.contessapetrini.workers.dev) | **Dynamic** (the fix in `b13778b`) |
| `backend/src/api/neurobotanica.js` | FastAPI backend | **0.8** (hardcoded) |
| `worker.js` / `worker.ts` | Unknown Worker | **0.5** (hardcoded base) |

#### What Needs to Happen

**Option A: Update the frontend API URL**

Change `NEXT_PUBLIC_BUDTENDER_API_URL` environment variable to:

```
https://terpene-api.contessapetrini.workers.dev
```

**Option B: Deploy the terpene-api Worker to a better URL**

```bash
cd workers/terpene-api
npx wrangler deploy
```

Then update the frontend to use that URL.

**Option C: Fix the FastAPI backend**

If the frontend should call the Python backend instead, remove hardcoded values from `backend/src/api/neurobotanica.js`.

#### Why Previous Deploys Didn't Work

1. Commit `b13778b` changed `workers/terpene-api/src/index.ts` ✔️
2. Pages deployment uploaded frontend to [`budtender.neuro-botanica.com`](http://budtender.neuro-botanica.com) ✔️
3. Frontend still calls [`budtender.neuro-botanica.com`](http://budtender.neuro-botanica.com) for API ❌
4. That's NOT where the Worker API is deployed ❌
5. Synergy values remain hardcoded ❌

### ⚠️ CI Test Failures (Separate Issue)

The Python tests are failing due to **missing database tables** in the test environment — NOT because of frontend changes. Tables missing:

- `clinical_studies`
- `patients`
- `neurobotanica_demographic_factors`
- `neurobotanica_drug_interactions`
- `neurobotanica_compounds`

**This is a CI configuration issue** that needs the test database to be migrated before pytest runs. It is separate from the Budtender app bugs.

### Deployment Reminder

**You cannot deploy.** Only push code to GitHub and provide the commit hash. The human operator will:

1. Verify your commit exists at https://github.com/cloakquillresearch-2026/neurobotanica_mvp/commits/main
2. Run the manual deploy command
3. Test the changes

Do NOT say "changes are live" or "you can see the fix at [budtender.neuro-botanica.com](http://budtender.neuro-botanica.com)" — you have no way to verify this.

---

## Application Overview

The **Budtender Education Assistant** is a tablet-optimized web application used by dispensary staff to:

1. Create and manage customer profiles
2. Capture customer health information and preferences
3. Generate personalized cannabis product recommendations
4. Educate budtenders on terpene profiles and therapeutic effects

### Tech Stack

- **Frontend:** Next.js 14, TypeScript, Tailwind CSS
- **Backend:** FastAPI (Python), Cloudflare Workers (TypeScript)
- **Database:** Cloudflare D1 (SQLite)
- **Deployment:** Cloudflare Pages

### Repository

`cloakquillresearch-2026/neurobotanica_mvp`

---

## ⚠️ CRITICAL: Deployment Configuration

> **RESOLVED (Jan 31, 2026):** The Budtender app is now deployed via manual Wrangler CLI deployments. After pushing code fixes to GitHub, run the deploy command locally to see changes live.
> 

### Current Cloudflare Pages Projects

| Project | Git Connection | What It Deploys |
| --- | --- | --- |
| `neurobotanica-mvp` | ✅ Connected | Marketing site (auto-deploys) |
| `neurobotanicabudtender` | Manual (Wrangler CLI) | Budtender app ✅ **Now working** |

### Production URLs

| URL | Points To | Status |
| --- | --- | --- |
| [`budtender.neuro-botanica.com`](http://budtender.neuro-botanica.com) | Budtender app | ✅ Working (manual deploy) |
| [`neurobotanicabudtender.pages.dev`](http://neurobotanicabudtender.pages.dev) | Budtender app | ✅ Working (manual deploy) |
| [`production.neurobotanicabudtender.pages.dev`](http://production.neurobotanicabudtender.pages.dev) | Budtender app | ✅ Working (manual deploy) |
| [`neurobotanica-mvp.pages.dev`](http://neurobotanica-mvp.pages.dev) | Marketing site | ✅ Receives updates (wrong app) |

---

## 🚨 DEPLOYMENT WORKFLOW (READ THIS FIRST)

### Why This Matters

The Budtender app uses **manual deployment**. Code pushed to GitHub does **NOT** automatically appear on the live site. A human must run a deploy command after every code change.

### The Complete Workflow

**Step 1: Coding Agent Makes Changes**

- Agent pushes commits to `cloakquillresearch-2026/neurobotanica_mvp`
- Agent provides the **exact commit hash** (e.g., `a1b2c3d4e5f6...`)
- Agent **cannot deploy** — they can only push code

**Step 2: Human Verifies the Commit Exists**

- Go to: https://github.com/cloakquillresearch-2026/neurobotanica_mvp/commits/main
- Find the commit hash the agent claimed
- **If the commit doesn't exist, the agent lied** — ask them to actually push the code

**Step 3: Human Deploys to Cloudflare**

Open terminal and run these commands:

```bash
# Navigate to the repo (clone first if needed)
cd neurobotanica_mvp

# Pull the latest changes
git pull origin main

# Build the frontend
cd frontend
npm install
npm run build

# Deploy to Cloudflare Pages
npx wrangler pages deploy ./out --project-name=neurobotanicabudtender
```

**Step 4: Verify the Deployment**

- Wrangler will output a deployment URL (e.g., [`https://abc123.neurobotanicabudtender.pages.dev`](https://abc123.neurobotanicabudtender.pages.dev))
- Visit that URL to confirm changes are live
- Also check: https://budtender.neuro-botanica.com

---

### Why We Can't Use Automatic Deployment

| Constraint | Explanation |  |
| --- | --- | --- |
| **Cloudflare limit** | One Git connection per repository |  |
| **Already connected** | `neurobotanica_mvp` → `neurobotanica-mvp` project (marketing site) |  |
| **Budtender app** | Uses Direct Upload mode (Wrangler CLI), not Git |  |

**Future fix:** Once the app works, split into a separate repository with automatic deployments.

---

### What Coding Agents CAN Do

- ✅ Read and analyze code in the repository
- ✅ Push commits to GitHub
- ✅ Create and merge pull requests
- ✅ Run tests (if configured in CI)
- ✅ Provide the commit hash for verification

### What Coding Agents CANNOT Do

- ❌ Run `npm run build` on your machine
- ❌ Run `wrangler pages deploy`
- ❌ Access Cloudflare Dashboard
- ❌ Verify changes are live on the actual website
- ❌ **Claim changes are "live" or "deployed"** (only a human can deploy)

---

### For Coding Agents: Honest Reporting Template

When you complete work, report like this:

> ✅ **Code pushed to GitHub**
> 

> - Commit: `[full 40-character hash]`
> 

> - Branch: `main`
> 

> - Files changed: `[list them]`
> 

> 
> 

> ⏳ **Awaiting manual deployment**
> 

> The human operator must run:
> 

> `
> 

> cd frontend && npm run build && npx wrangler pages deploy ./out --project-name=neurobotanicabudtender
> 

> `
> 

> Settings → Builds & deployments → Connect to Git
> 

> 🔍 **What to test after deployment:**
> 

> - [specific feature to verify]
> 

> - [expected behavior]
> 

**Do NOT say:**

- "Changes are now live"
- "You can see the fix at [budtender.neuro-botanica.com](http://budtender.neuro-botanica.com)"
- "Deployment successful"

---

### Quick Reference Commands

**First-time setup:**

```bash
git clone https://github.com/cloakquillresearch-2026/neurobotanica_mvp.git
cd neurobotanica_mvp/frontend
npm install
npm install -g wrangler
wrangler login
```

**Subsequent deployments:**

```bash
cd neurobotanica_mvp
git pull origin main
cd frontend
npm run build
npx wrangler pages deploy ./out --project-name=neurobotanicabudtender
```

**Verify deployment worked:**

- Check Cloudflare Dashboard → Pages → neurobotanicabudtender → Deployments
- Visit https://budtender.neuro-botanica.com

---

## Workflow 1: New Client Registration

### 🚨 CRITICAL: Read This Carefully

The current UI is **completely wrong**. The agent's "fix" made it worse. Here is exactly what must happen:

---

### Correct UI Flow (Step by Step)

**Step 1: User clicks "+ New Client"**

| What SHOULD happen | What currently happens (WRONG) |
| --- | --- |
| A modal or new screen appears with ONLY the registration form | The main consultation view stays visible |
| All analysis widgets are HIDDEN | "TS-PS-001 Cross-Kingdom Insights" starts computing |
| No conditions are selectable yet | Condition buttons are visible and clickable |

**Step 2: User sees the New Client Form**

The form must appear in this EXACT order, top to bottom:

```
┌─────────────────────────────────────┐
│  + NEW CLIENT                       │
├─────────────────────────────────────┤
│  First Name *        [___________]  │  ← REQUIRED, at the TOP
│  Last Name *         [___________]  │  ← REQUIRED
│  Age                 [___________]  │  ← optional
│  Gender              [dropdown___]  │  ← optional
│  Weight (kg)         [___________]  │  ← optional
│  Notes               [             ]│  ← optional, multiline
│                      [             ]│
├─────────────────────────────────────┤
│  [Cancel]              [Save Client]│
└─────────────────────────────────────┘
```

**What must NOT be visible on this screen:**

- ❌ Experience Level selector (First Time, Beginner, etc.)
- ❌ Primary Conditions (Chronic Pain, Anxiety, etc.)
- ❌ "Run Analysis" button
- ❌ TS-PS-001 Cross-Kingdom Insights
- ❌ NeuroBotanica Analysis panel
- ❌ Terpene Education panel
- ❌ Any "Computing..." indicators

**Step 3: User fills in name and clicks "Save Client"**

- Validate that First Name and Last Name are not empty
- Call POST `/api/dispensary/profile` to create customer in database
- Receive back a `customer_id`
- THEN navigate to the consultation view WITH the customer loaded

**Step 4: ONLY NOW show the consultation workflow**

After the customer is saved, the user should see:

- Customer name displayed at the top (e.g., "Consultation for: Jane Doe")
- Experience Level selector
- Primary Conditions selector
- "Run Analysis" button (disabled until at least one condition is selected)

---

### What The Agent Did Wrong

1. **Put name fields at the BOTTOM** — They should be at the TOP, before anything else
2. **Kept analysis widgets visible** — The "TS-PS-001 Computing..." should not appear until analysis is actually requested
3. **Didn't create a separate form screen** — The New Client button should open a dedicated form, not just add fields to the existing view
4. **Made condition buttons open forms** — The condition buttons (🔥 Chronic Pain, 🧘 Anxiety) are for selecting conditions during consultation, NOT for triggering new client registration

---

### The Correct Mental Model

Think of it like a medical office:

1. **Reception desk** (New Client form) — Collect name first
2. **Patient intake** (Experience Level, Conditions) — Collect health info second
3. **Doctor consultation** (Run Analysis) — Generate recommendations third

You cannot skip step 1 and go straight to step 2 or 3.

---

### Files to Modify

- `frontend/src/pages/index.tsx` — Create separate "new client" view state
- `frontend/src/components/CustomerProfile.tsx` — Ensure name fields are at TOP
- `frontend/src/components/NewClientModal.tsx` — May need to create this component
- `backend/routers/[dispensary.py](http://dispensary.py)` — POST `/api/dispensary/profile`

---

## Workflow 2: Returning Customer Lookup

### Expected Behavior

1. **Search box** should query the database for existing customers by:
    - First name
    - Last name
    - Customer ID code
2. **Results should only show customers created by this dispensary** — not test data or mock data
3. **Selecting a customer** should:
    - Load their full profile (name, conditions, preferences, notes)
    - Display their consultation history
    - Pre-populate the recommendation engine with their data

### Current Bug

Search returns hardcoded test data ("John Smith") that was never entered by the user.

### Files to Modify

- `frontend/src/utils/api.ts` — customer search endpoint
- `backend/routers/[dispensary.py](http://dispensary.py)` — GET `/api/dispensary/customers/search`
- Remove any mock/seed data from production database

---

## Workflow 3: Customer Consultation

### Step 3.1: Customer Needs Assessment

After selecting or creating a customer, the budtender should:

1. **Select Experience Level** (required):
    - First Time
    - Beginner
    - Occasional
    - Regular
    - Expert
2. **Select Primary Conditions** (at least one required):
    - Chronic Pain
    - Anxiety
    - Insomnia
    - PTSD
    - Depression
    - Nausea
    - Appetite Loss
    - Inflammation
    - (others as defined in database)
3. **Optional: Inflammatory Biomarkers** (for TS-PS-001 personalization):
    - TNF-α (pg/mL)
    - IL-6 (pg/mL)
    - CRP (mg/L)
    - IL-1β (pg/mL)

### Step 3.2: Run Analysis

When "Run Analysis" is clicked:

1. **Call the NeuroBotanica API** with:
    - Customer profile data
    - Selected conditions
    - Experience level
    - Biomarkers (if provided)
2. **Display Results** including:
    - **Primary Kingdom recommendation** (cannabis, fungal, marine, plant)
    - **Synergy Score** — MUST be dynamically calculated based on:
        - Compound combinations
        - Customer conditions
        - Biomarker data (if available)
        - **NOT hardcoded to 50.0%**
    - **Microbiome Modulation** — MUST be dynamically calculated:
        - Based on selected compounds
        - CBD should show higher scores than baseline
        - **NOT hardcoded to 75.0%**
    - **Recommended Compounds** (e.g., CBD, CBG, THC ratios)
    - **Cannabinoid Profile** with actual THC% and CBD% values
        - **NOT 0%/0%**

---

## Workflow 4: Product Recommendations

### Expected Behavior

1. **Recommendations API** should return:
    - Product name
    - THC/CBD percentages (actual values, not 0%)
    - Dominant terpenes
    - Suggested dosing
    - Evidence tier (clinical, research, traditional knowledge)
2. **Display should include**:
    - Product cards with terpene tags
    - Match score to customer needs
    - Educational talking points for the budtender

### Current Bug

Recommendation title shows "undefined Evidence-based Recommendation"

---

## Workflow 5: Notes Persistence

### Expected Behavior

1. **Notes field** in customer profile should:
    - Accept text input
    - Save to database when profile is saved
    - Persist across sessions
    - Load when customer is selected again
2. **Care Notes** section should:
    - Display previously saved notes
    - Allow editing
    - Save changes to database

### Current Bug

Notes are not persisting to the database.

### Files to Modify

- `frontend/src/components/CustomerProfile.tsx` — notes field onChange handler
- `backend/routers/[dispensary.py](http://dispensary.py)` — include notes in profile CRUD operations
- `unified_schema.sql` — ensure `notes` column exists in customer table

---

## API Endpoints Required

### Customer Management

| Endpoint | Method | Purpose |
| --- | --- | --- |
| `/api/dispensary/profile` | POST | Create new customer |
| `/api/dispensary/profile/{id}` | GET | Get customer by ID |
| `/api/dispensary/profile/{id}` | PUT | Update customer |
| `/api/dispensary/customers/search` | GET | Search customers |

### Analysis & Recommendations

| Endpoint | Method | Purpose |
| --- | --- | --- |
| `/api/dispensary/analyze` | POST | Run full NeuroBotanica analysis |
| `/api/dispensary/recommendations` | POST | Get product recommendations |
| `/api/dispensary/inflammatory-synergy` | POST | TS-PS-001 synergy prediction |

---

## Dynamic Calculation Requirements

### Compound Synergy Score

**MUST NOT** return hardcoded values. Calculate based on:

```
base_score = 0.5
if (compound_a == 'CBD' && compound_b == 'THC'):
    score += 0.3  // CBD-THC entourage effect
if (has_terpenes):
    score += terpene_synergy_modifier
if (biomarkers_elevated):
    score += inflammation_response_modifier
return min(score, 1.0)
```

### Microbiome Modulation

**MUST NOT** return hardcoded 75%. Calculate based on:

```
base_confidence = 0.6
if (compound == 'CBD'):
    confidence += 0.15
if (compound == 'CBG'):
    confidence += 0.10
return min(confidence, 1.0)
```

### THC/CBD Profile

**MUST NOT** show 0%/0%. Either:

- Pull actual values from product database
- Calculate recommended ratios based on condition
- Show "N/A" if no data available (never show 0%/0%)

---

## Testing Checklist

Before claiming any feature is "working", verify:

- [ ]  **New Client** — Opens blank form with name fields, does NOT auto-run analysis
- [ ]  **Customer Search** — Returns only real customers, no test data
- [ ]  **Notes** — Saves to database and loads on next visit
- [ ]  **Synergy Score** — Shows different values for different compound combinations
- [ ]  **Microbiome Modulation** — CBD shows higher than baseline (not always 75%)
- [ ]  **THC/CBD Profile** — Shows actual percentages, not 0%/0%
- [ ]  **Recommendations** — Title shows product name, not "undefined"
- [ ]  **All API endpoints** — Return proper data, not errors

---

## File Structure Reference

```
neurobotanica_mvp/
├── frontend/
│   ├── src/
│   │   ├── components/
│   │   │   ├── CustomerProfile.tsx      # Customer form
│   │   │   ├── RecommendationCard.tsx   # Product display
│   │   │   └── SynergyDisplay.tsx       # Synergy scores
│   │   ├── pages/
│   │   │   └── index.tsx                # Main consultation page
│   │   └── utils/
│   │       └── api.ts                   # API client
├── backend/
│   ├── main.py                          # FastAPI app
│   ├── routers/
│   │   └── dispensary.py                # Dispensary endpoints
│   └── services/
│       └── inflammatory_synergy_engine.py
├── workers/
│   └── terpene-api/
│       └── src/
│           └── index.ts                 # Cloudflare Worker
└── unified_schema.sql                   # Database schema
```

---

## Summary of Current Bugs

| Bug | Location | Priority |
| --- | --- | --- |
| New Client skips form | `frontend/src/pages/index.tsx` | P1 |
| Search returns test data | `backend/routers/[dispensary.py](http://dispensary.py)` | P1 |
| Synergy hardcoded 50% | `workers/terpene-api/src/index.ts` | P1 |
| Microbiome hardcoded 75% | `workers/terpene-api/src/index.ts` | P1 |
| THC/CBD shows 0%/0% | `frontend/src/components/RecommendationCard.tsx` | P1 |
| Notes not persisting | `frontend/src/components/CustomerProfile.tsx` | P2 |
| "undefined" in recommendation | `frontend/src/components/RecommendationCard.tsx` | P2 |

---

## Instructions for Coding Agents

1. **Read this entire document** before making any changes
2. **Do not claim features work** without testing against the checklist above
3. **Report honestly** if something doesn't work
4. **Commit real code** — don't fabricate commit hashes
5. **Test in the actual browser** at [`production.neurobotanicabudtender.pages.dev`](http://production.neurobotanicabudtender.pages.dev)
6. **Screenshot your results** if possible

The goal is a **working application**, not claims of a working application.