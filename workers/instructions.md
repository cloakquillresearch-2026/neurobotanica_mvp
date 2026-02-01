# Budtender Application — Complete Workflow Specification

## Purpose

This document provides **complete, unambiguous specifications** for the NeuroBotanica Budtender Education Assistant. Any coding agent working on this application must implement these workflows exactly as described.

> **For Coding Agents:** Do not claim features are "working" unless you have verified them against this specification. Do not fabricate test results. If something doesn't work, report it honestly.
> 

---

## 🚨 CURRENT SITUATION (Feb 1, 2026)

<aside>
⚠️

**READ THIS FIRST.** A previous coding agent made changes that broke the app further and did not follow the workflow specification. You may need to revert their changes before proceeding.

</aside>

### What the Previous Agent Did Wrong

| Problem | Details |
| --- | --- |
| **Fabricated deployment claims** | Said "Status: ✅ Live" when agents cannot deploy |
| **Misunderstood the UI flow** | Made condition buttons (🔥 Chronic Pain, 🧘 Anxiety) open registration forms — those buttons are for selecting conditions during consultation, NOT for creating new clients |
| **Put form fields in wrong order** | Name fields are at the BOTTOM instead of the TOP |
| **Didn't create a separate form** | The "+ New Client" button should open a modal/separate screen, not add fields to the existing view |
| **Left analysis running on load** | "TS-PS-001 Cross-Kingdom Insights" shows "Computing..." before any customer is selected |
| **Took 3 attempts to actually commit** | Initially claimed commits existed when they didn't |

### Commits That May Need Reverting

- `71a7ed9` — Contains the broken "fixes" described above
- Check the commit history to understand what was changed

### Your First Task

**Fix Workflow 1: New Client Registration** — This is the most broken part.

When a user clicks "+ New Client":

1. A modal or new screen should appear with ONLY the registration form
2. Name fields (First Name, Last Name) must be at the TOP
3. No analysis should run until the user explicitly clicks "Run Analysis" later
4. The consultation view (Experience Level, Conditions, Analysis panels) should be HIDDEN until after the client is saved

Read the full specification in **Workflow 1** below.

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