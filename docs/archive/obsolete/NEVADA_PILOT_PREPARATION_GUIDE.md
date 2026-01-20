# Nevada Pilot Preparation Guide

**Target Launch:** Week of March 26, 2026  
**Organization:** Cloak and Quill Research 501(c)(3)  
**Document Status:** DRAFT - December 23, 2025

---

## Executive Summary

The Nevada pilot represents NeuroBotanica's first production deployment, targeting licensed cannabis dispensaries in Nevada. This guide outlines the regulatory, technical, and partnership requirements for a successful launch.

---

## 1. Nevada Cannabis Regulatory Landscape

### 1.1 Regulatory Authority

**Primary Agency:** Nevada Cannabis Compliance Board (CCB)  
**Website:** https://ccb.nv.gov/  
**Key Regulations:** NAC 678A, NAC 678B, NAC 678D

### 1.2 License Types Relevant to NeuroBotanica

| License Type | Description | NeuroBotanica Relevance |
|--------------|-------------|------------------------|
| Retail Store | Dispensary license | Primary customer target |
| Medical Dispensary | Medical-only retail | Secondary target |
| Consumption Lounge | On-site consumption | Future expansion |
| Testing Laboratory | Product testing | Data partnership opportunity |

### 1.3 Key Compliance Requirements

#### For NeuroBotanica (Software Provider)
- [ ] **No cannabis license required** - We are a software/data provider, not handling cannabis
- [ ] **HIPAA consideration** - If collecting patient data, ensure BAA compliance
- [ ] **Data residency** - Nevada may require patient data stored in-state or US
- [ ] **SOC 2 Type II** - Recommended for enterprise dispensary sales

#### For Dispensary Partners
- [ ] Must maintain valid Nevada retail/medical license
- [ ] Must use state-approved seed-to-sale tracking (Metrc)
- [ ] Patient data handling must comply with NRS 453A (medical)

---

## 2. Partnership Development

### 2.1 Target Dispensary Profiles

**Ideal Pilot Partner Characteristics:**
- Annual revenue: $2M-$10M (mid-size, growth-focused)
- Locations: 2-5 stores (enough data, not too complex)
- Tech-forward: Already using POS analytics, interested in differentiation
- Patient-focused: Medical dispensary or dual-use with medical emphasis
- Geographic spread: Las Vegas metro + one rural location for diversity

**Suggested Initial Targets:**

| Tier | Type | Approach | Timeline |
|------|------|----------|----------|
| Tier 1 | 3 mid-size dispensaries | Direct outreach | Jan 2026 |
| Tier 2 | 2 large chains | Conference/referral | Feb 2026 |
| Tier 3 | 5 boutique/medical | Industry association | Mar 2026 |

### 2.2 Value Proposition for Dispensaries

**For Dispensary Operators:**
```
"NeuroBotanica gives your budtenders science-backed recommendations 
in seconds—no more guesswork. Our AI analyzes 368 clinical studies 
across 22 conditions to match patients with the right products, 
increasing customer satisfaction and repeat visits."
```

**Key Selling Points:**
1. **Liability reduction** - Evidence-backed recommendations, not staff opinions
2. **Staff training shortcut** - New budtenders productive in days, not months
3. **Customer retention** - Personalized recommendations drive loyalty
4. **Premium positioning** - "Science-based dispensary" marketing angle
5. **Compliance documentation** - Audit trail for all recommendations

### 2.3 Partnership Agreement Terms

**Pilot Program Structure (Suggested):**

| Term | Details |
|------|---------|
| Duration | 90-day pilot with 12-month renewal option |
| Pricing | Free during pilot; post-pilot $500-2,500/month based on volume |
| Data sharing | Anonymized outcome data shared with NeuroBotanica for model improvement |
| Exclusivity | Non-exclusive; dispensary can use other tools |
| Exit clause | Either party can exit with 30-day notice |
| Support | Dedicated onboarding + weekly check-ins during pilot |

### 2.4 Outreach Strategy

**Phase 1: Research (January 1-15, 2026)**
- [ ] Build list of 50 licensed Nevada dispensaries from CCB public records
- [ ] Research each: size, ownership, tech stack, patient focus
- [ ] Identify decision-makers (GM, Owner, Compliance Officer)
- [ ] Prioritize top 10 targets

**Phase 2: Initial Contact (January 15-31, 2026)**
- [ ] Send personalized intro emails to top 10
- [ ] Follow up with phone calls
- [ ] Offer free demo/consultation
- [ ] Target: 5 demo meetings scheduled

**Phase 3: Demos & Negotiation (February 2026)**
- [ ] Conduct personalized demos for interested parties
- [ ] Address compliance/integration questions
- [ ] Negotiate pilot terms
- [ ] Target: 3 signed pilot agreements

**Phase 4: Onboarding (March 1-25, 2026)**
- [ ] Technical integration with each pilot partner
- [ ] Staff training sessions
- [ ] Soft launch and testing
- [ ] Go-live: March 26, 2026

---

## 3. Technical Requirements

### 3.1 Infrastructure for Nevada Pilot

**Minimum Production Requirements:**

| Component | Specification | Provider Options |
|-----------|---------------|------------------|
| Hosting | US-based cloud, SOC 2 compliant | AWS us-west-2, GCP, Azure |
| Database | PostgreSQL 15+, encrypted at rest | AWS RDS, Supabase, Railway |
| API | HTTPS only, rate-limited | Current FastAPI setup |
| Auth | OAuth 2.0 / JWT tokens | Auth0, Clerk, or OmniPath |
| Monitoring | APM + error tracking | Datadog, Sentry, New Relic |
| Backup | Daily automated, 30-day retention | AWS S3, Backblaze |

**Architecture Diagram:**

```
┌─────────────────────────────────────────────────────────────────┐
│                     Nevada Production Environment                │
├─────────────────────────────────────────────────────────────────┤
│                                                                  │
│  ┌──────────────┐     ┌──────────────┐     ┌──────────────┐    │
│  │ Dispensary   │────▶│ API Gateway  │────▶│ NeuroBotanica│    │
│  │ POS System   │     │ (Auth + Rate │     │ API (FastAPI)│    │
│  │              │     │  Limiting)   │     │              │    │
│  └──────────────┘     └──────────────┘     └──────────────┘    │
│                                                    │            │
│                                                    ▼            │
│                                            ┌──────────────┐    │
│                                            │  PostgreSQL  │    │
│                                            │  (Encrypted) │    │
│                                            └──────────────┘    │
│                                                    │            │
│                                                    ▼            │
│                                            ┌──────────────┐    │
│                                            │ ML Models    │    │
│                                            │ (Read-only)  │    │
│                                            └──────────────┘    │
│                                                                  │
└─────────────────────────────────────────────────────────────────┘
```

### 3.2 Integration Options for Dispensaries

**Option A: API Integration (Recommended for tech-forward dispensaries)**
- Dispensary's POS/EHR calls NeuroBotanica API directly
- Real-time recommendations at point of sale
- Requires technical staff at dispensary or POS vendor cooperation

**Option B: Web Dashboard (Recommended for pilot)**
- Standalone web app for budtenders
- No POS integration required
- Copy/paste or manual entry workflow
- Lower barrier to adoption

**Option C: Tablet Kiosk (Future)**
- Patient self-service intake
- Recommendations before budtender consultation
- Higher hardware cost

**Pilot Recommendation:** Start with Option B (Web Dashboard) for all pilot partners. Add API integration for partners who request it.

### 3.3 Frontend Requirements for Pilot

**Minimum Viable Dashboard Features:**

| Feature | Priority | Description |
|---------|----------|-------------|
| Patient intake form | P0 | Collect condition, preferences, history |
| Product recommendation | P0 | Display matched products with rationale |
| Recommendation history | P1 | View past recommendations for returning patients |
| Staff login | P0 | Basic auth for budtenders |
| Admin panel | P1 | Manage products, view analytics |
| Export/print | P2 | PDF recommendation summary for patient |

**Tech Stack Recommendation:**
- **Framework:** Next.js 14 (React) or SvelteKit
- **Styling:** Tailwind CSS
- **Auth:** Clerk or Auth0 (fast integration)
- **Hosting:** Vercel or Cloudflare Pages

---

## 4. Data Collection & RWE Framework

### 4.1 Real-World Evidence (RWE) Collection

The Nevada pilot is critical for collecting real-world evidence to improve model accuracy. 

**Data Points to Collect (Anonymized):**

| Category | Data Point | Collection Method |
|----------|------------|-------------------|
| Patient | Age range, gender, conditions | Intake form |
| Product | Product type, cannabinoid profile, terpenes | Dispensary inventory |
| Recommendation | What was recommended, confidence score | System log |
| Outcome | Did patient purchase? Return visit? | POS integration or follow-up |
| Feedback | Patient satisfaction (1-5), symptom improvement | Optional follow-up survey |

### 4.2 Outcome Tracking

**Key Metrics for Pilot Success:**

| Metric | Target | Measurement |
|--------|--------|-------------|
| Recommendation acceptance rate | >60% | Purchases matching recommendations |
| Patient return rate | >40% in 90 days | Repeat visits |
| Staff adoption | >80% | % of transactions using NeuroBotanica |
| NPS score | >50 | Patient survey |
| Model accuracy improvement | +5% | Pre/post pilot comparison |

### 4.3 Privacy & Consent

**Required Consent Language (Example):**

```
"NeuroBotanica uses anonymized health information to provide 
personalized product recommendations. Your data is never sold 
and is protected under our Privacy Policy. By proceeding, you 
consent to the collection and use of your information as described.

□ I consent to NeuroBotanica's data collection and use
□ I would like to receive follow-up surveys about my experience
```

**Data Handling Requirements:**
- [ ] No PII stored without explicit consent
- [ ] All data encrypted in transit (TLS 1.3) and at rest (AES-256)
- [ ] 30-day data retention for non-consented sessions
- [ ] CCPA compliance for California residents visiting Nevada

---

## 5. Timeline & Milestones

### 5.1 Detailed Timeline

```
January 2026
├── Week 1 (Jan 1-7)
│   ├── Finalize dispensary target list
│   ├── Draft partnership agreement template
│   └── Complete production infrastructure setup
│
├── Week 2 (Jan 8-14)
│   ├── Begin outreach to Tier 1 dispensaries
│   ├── Prepare demo environment
│   └── Create dispensary onboarding documentation
│
├── Week 3 (Jan 15-21)
│   ├── First demo meetings
│   ├── Iterate on pitch based on feedback
│   └── Frontend development kickoff
│
└── Week 4 (Jan 22-31)
    ├── Follow-up with interested parties
    ├── Begin pilot agreement negotiations
    └── Continue frontend development

February 2026
├── Week 5-6 (Feb 1-14)
│   ├── Sign first pilot agreements
│   ├── Complete MVP frontend
│   └── Begin integration testing
│
├── Week 7 (Feb 15-21)
│   ├── Staff training content creation
│   ├── Soft launch with first partner
│   └── Bug fixes and iteration
│
└── Week 8 (Feb 22-28)
    ├── Expand soft launch to all pilot partners
    ├── Collect initial feedback
    └── Final production hardening

March 2026
├── Week 9-10 (Mar 1-14)
│   ├── Monitor soft launch performance
│   ├── Address any critical issues
│   └── Prepare marketing materials
│
├── Week 11 (Mar 15-21)
│   ├── Final QA pass
│   ├── Load testing
│   └── Backup/recovery testing
│
└── Week 12 (Mar 22-28)
    ├── LAUNCH: March 26, 2026
    ├── On-site support for pilot partners
    └── Begin outcome tracking
```

### 5.2 Go/No-Go Criteria

**Launch Readiness Checklist (March 25, 2026):**

| Category | Requirement | Status |
|----------|-------------|--------|
| **Technical** | | |
| | API uptime >99.5% over past 14 days | ⬜ |
| | All P0 frontend features complete | ⬜ |
| | SSL/TLS configured, security scan passed | ⬜ |
| | Backup/recovery tested | ⬜ |
| | Monitoring/alerting configured | ⬜ |
| **Business** | | |
| | ≥2 signed pilot agreements | ⬜ |
| | Staff training completed for all partners | ⬜ |
| | Support escalation path defined | ⬜ |
| | Legal review of terms complete | ⬜ |
| **Compliance** | | |
| | Privacy policy published | ⬜ |
| | Consent flow implemented | ⬜ |
| | Data handling documented | ⬜ |

---

## 6. Budget Estimate

### 6.1 Pilot Phase Costs (Jan-Mar 2026)

| Category | Item | Cost |
|----------|------|------|
| **Infrastructure** | | |
| | Cloud hosting (3 months) | $300-600 |
| | Database (managed PostgreSQL) | $100-200 |
| | Domain + SSL | $50 |
| | Monitoring (Sentry/Datadog) | $0-100 |
| **Development** | | |
| | Frontend development (contractor) | $5,000-10,000 |
| | DevOps/deployment (contractor) | $2,000-3,000 |
| **Marketing** | | |
| | Trade show/conference attendance | $500-1,500 |
| | Collateral (brochures, business cards) | $200-500 |
| | Demo video production | $500-1,000 |
| **Legal** | | |
| | Partnership agreement review | $500-1,000 |
| | Privacy policy/Terms of Service | $500-1,000 |
| **Travel** | | |
| | Trips to Nevada for partner meetings | $1,000-2,000 |
| | On-site launch support | $500-1,000 |
| **TOTAL** | | **$11,150-22,400** |

### 6.2 Post-Pilot Revenue Projection

**Conservative Scenario (3 dispensaries, 6 months):**

| Month | Dispensaries | Revenue/Each | Total MRR |
|-------|--------------|--------------|-----------|
| Apr 2026 | 3 | $500 (pilot pricing) | $1,500 |
| May 2026 | 4 | $750 | $3,000 |
| Jun 2026 | 5 | $1,000 | $5,000 |
| Jul 2026 | 7 | $1,000 | $7,000 |
| Aug 2026 | 10 | $1,250 | $12,500 |
| Sep 2026 | 12 | $1,500 | $18,000 |

**6-Month Total:** ~$47,000

---

## 7. Risk Mitigation

### 7.1 Key Risks & Mitigations

| Risk | Likelihood | Impact | Mitigation |
|------|------------|--------|------------|
| No dispensaries sign pilot | Medium | Critical | Start outreach early; offer generous pilot terms |
| Technical issues at launch | Medium | High | Extensive QA; soft launch period; on-call support |
| Regulatory change in Nevada | Low | High | Monitor CCB closely; design for compliance |
| Competitor enters market | Medium | Medium | Move fast; focus on differentiation (TK, 6 engines) |
| Staff adoption resistance | Medium | Medium | Invest in training; make UX dead simple |
| Data breach | Low | Critical | SOC 2 prep; encryption; access controls |

### 7.2 Contingency Plans

**If no dispensaries sign by Feb 15:**
- Pivot to offering free perpetual tier with limited features
- Partner with cannabis industry consultants for warm introductions
- Consider MJBizCon or other trade show presence

**If technical issues at launch:**
- Designated on-call engineer for first 2 weeks
- Rollback plan to previous stable version
- Direct communication channel with pilot partners

---

## 8. Success Metrics

### 8.1 Pilot Success Criteria

**By March 31, 2026 (End of Pilot Month 1):**
- [ ] ≥3 dispensaries actively using platform
- [ ] ≥500 patient recommendations generated
- [ ] ≥60% recommendation acceptance rate
- [ ] Staff NPS ≥40

**By June 30, 2026 (End of Q2):**
- [ ] ≥10 paying customers
- [ ] ≥$5,000 MRR
- [ ] ≥5,000 patient recommendations generated
- [ ] Model accuracy improved by ≥5% from baseline
- [ ] ≥2 community healer validations integrated

---

## 9. Action Items

### Immediate (This Week)

1. [ ] Review and customize this guide for Cloak and Quill's specific situation
2. [ ] Identify point person for dispensary outreach
3. [ ] Begin infrastructure setup (cloud account, database)
4. [ ] Create simple one-pager for dispensary outreach

### Next 30 Days

5. [ ] Build dispensary target list from CCB records
6. [ ] Draft partnership agreement (have legal review)
7. [ ] Design MVP frontend wireframes
8. [ ] Prepare demo environment

### Before Launch

9. [ ] Sign ≥2 pilot agreements
10. [ ] Complete frontend development
11. [ ] Conduct staff training for all partners
12. [ ] Complete go/no-go checklist

---

**Document Owner:** Dr. Contessa Petrini  
**Last Updated:** December 23, 2025  
**Next Review:** January 15, 2026

---

*This document is CONFIDENTIAL and intended for internal use by Cloak and Quill Research.*
