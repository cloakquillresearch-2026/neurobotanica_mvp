# Budtender Education Assistant - Implementation Manual

## Overview
This manual provides instructions for deploying the NeuroBotanica Budtender Education Assistant tablets in dispensaries. The system uses a **cloud-hosted API model** - all AI processing and trade secret algorithms remain on secure Cloak and Quill servers.

**Security Model**: No source code distribution. Dispensaries receive pre-configured tablets that connect to our secure API.

## Architecture

```
┌─────────────────────────────────────────────────────────────────┐
│                    DISPENSARY (Client-Side)                     │
│  ┌─────────────────┐                                            │
│  │  Tablet Device  │  ← Pre-installed app (PWA or native)       │
│  │  (iPad/Android) │                                            │
│  └────────┬────────┘                                            │
└───────────┼─────────────────────────────────────────────────────┘
            │ HTTPS API Calls Only
            ▼
┌─────────────────────────────────────────────────────────────────┐
│              CLOAK AND QUILL SERVERS (Secure)                   │
│  ┌─────────────────────────────────────────────────────────────┐│
│  │  Cloudflare Workers (285+ global edge locations)            ││
│  │  ├── Authentication & Rate Limiting                         ││
│  │  ├── Request Validation                                     ││
│  │  └── Response Caching                                       ││
│  └─────────────────────────────────────────────────────────────┘│
│  ┌─────────────────────────────────────────────────────────────┐│
│  │  Private Backend (NEVER exposed)                            ││
│  │  ├── Trade Secret Engines (BioPath, ClinPath, etc.)         ││
│  │  ├── ML Models (Therapeutic Prediction, Dimer Analysis)     ││
│  │  ├── Clinical Database (505+ studies)                       ││
│  │  └── Adjuvant Optimization Engine                           ││
│  └─────────────────────────────────────────────────────────────┘│
└─────────────────────────────────────────────────────────────────┘
```

## What Dispensaries Receive

| Item | Description |
|------|-------------|
| **Tablet Device** | Pre-configured iPad or Android tablet |
| **App Credentials** | API key + dispensary account login |
| **Training Materials** | Video tutorials + quick reference cards |
| **Support Access** | 24/7 technical support hotline |

**What They DON'T Receive:**
- ❌ Source code access
- ❌ GitHub repository access
- ❌ Direct database access
- ❌ Trade secret documentation
- ❌ ML model files

## Deployment Process

### Step 1: Dispensary Onboarding
1. Execute licensing agreement (includes NDA and trade secret protections)
2. Collect dispensary information (location, staff count, inventory system)
3. Generate unique API credentials for dispensary
4. Ship pre-configured tablet(s)

### Step 2: Device Setup (Performed by Cloak and Quill)
Before shipping, each tablet is configured with:
```
✓ NeuroBotanica PWA installed
✓ API endpoint configured (https://api.neurobotanica.com)
✓ Dispensary-specific API key embedded
✓ Offline cache enabled
✓ Auto-update enabled
✓ Device management enrolled (MDM)
```

### Step 3: Dispensary Receives Tablet
1. **Power on** - Device auto-connects to WiFi
2. **Login** - Staff enters dispensary credentials
3. **Verify** - Run test consultation to confirm API connectivity
4. **Train** - Complete 15-minute interactive tutorial

### Step 4: Go-Live Checklist
- [ ] Tablet connected to dispensary WiFi
- [ ] API connection verified (green status indicator)
- [ ] At least one staff member trained
- [ ] Emergency contact information on file
- [ ] Backup tablet shipped (optional)

## API Endpoints (What the Tablet Accesses)

| Endpoint | Purpose |
|----------|---------|
| `POST /api/dispensary/recommend` | Get AI recommendations |
| `POST /api/dispensary/profile` | Create/update customer profile |
| `POST /api/adjuvant/optimize` | Get adjuvant suggestions |
| `GET /api/dispensary/statistics` | View usage analytics |
| `POST /api/dispensary/feedback` | Submit recommendation feedback |

**Rate Limits:**
- 1000 API calls/day per dispensary (MVP tier)
- 10,000 API calls/day (Professional tier)

## Security Measures

### API Authentication
- Each dispensary has unique API key
- Keys rotated quarterly
- All requests logged with IP and timestamp
- Anomaly detection for suspicious patterns

### Device Security
- Mobile Device Management (MDM) enrolled
- Remote wipe capability
- App cannot be uninstalled without admin PIN
- Screen lock after 5 minutes inactivity

### Data Protection
- All API calls over HTTPS (TLS 1.3)
- No customer PII stored on device (session only)
- Customer data encrypted at rest on servers
- HIPAA-compliant infrastructure

## Troubleshooting

### Common Issues

**1. "API Connection Failed"**
- Check WiFi connection
- Verify tablet date/time is correct
- Restart app
- Contact support if persists

**2. "Recommendations Not Loading"**
- Check API status at status.neurobotanica.com
- Verify API rate limit not exceeded
- Try refreshing customer profile

**3. "App Not Responding"**
- Force close and reopen app
- Restart tablet
- Check for pending updates

### Support Channels
- **Phone**: 1-800-NEUROBOT (24/7)
- **Email**: support@neurobotanica.com
- **Portal**: help.neurobotanica.com
- **Emergency**: Direct line to on-call engineer

## Pricing Tiers

| Tier | Monthly Cost | API Calls | Tablets | Support |
|------|-------------|-----------|---------|---------|
| **Starter** | $199/mo | 1,000/day | 1 | Email |
| **Professional** | $499/mo | 10,000/day | 3 | Phone + Email |
| **Enterprise** | Custom | Unlimited | Unlimited | Dedicated |

## Compliance

### Nevada Regulatory Compliance
- Age verification prompts built-in
- Transaction logging for audit trails
- METRC-compatible data export (Enterprise tier)
- Monthly compliance reports

### HIPAA Compliance
- BAA available for medical dispensaries
- Encrypted data transmission
- Access logging and audit trails
- Annual security assessments

## Intellectual Property Notice

**IMPORTANT**: The NeuroBotanica system, including all AI models, algorithms, and trade secret implementations, is the exclusive property of Cloak and Quill Research 501(c)(3).

- Patent Portfolio: $684M-$1.026B valuation
- Protected under 18 U.S.C. § 1836 (DTSA)
- Unauthorized reverse engineering prohibited
- API access does not grant IP rights

Dispensary partners receive a **limited license** to use the service via API only. No ownership or IP rights are transferred.

---

*Cloak and Quill Research 501(c)(3) | Henderson, Nevada*
*Last Updated: December 24, 2025*
