# NeuroBotanica Deployment Guide

## December 27, 2025 - Production Deployment Strategy

### 🎯 Deployment Overview

**NeuroBotanica MVP** will be deployed using Cloudflare's edge computing infrastructure for optimal performance and global accessibility. This deployment strategy leverages:

- **Cloudflare Workers** for the FastAPI backend (285+ global locations, <200ms response times)
- **Cloudflare Pages** for the Next.js frontend (global CDN, automatic HTTPS)
- **Cloudflare D1** for the SQLite database (serverless, encrypted)
- **Cloudflare KV** for caching and session storage

---

## 🚀 Backend Deployment (Cloudflare Workers)

### Prerequisites
```bash
# Install Wrangler CLI
npm install -g wrangler

# Authenticate with Cloudflare
wrangler auth login

# Install Python dependencies for deployment
pip install -r requirements.txt
```

### Step 1: Create Cloudflare Worker Project
```bash
# Initialize Wrangler project
cd backend
wrangler init neurobotanica-api --type python

# Configure wrangler.toml
cat > wrangler.toml << EOF
name = "neurobotanica-api"
main = "main.py"
compatibility_date = "2024-01-01"
compatibility_flags = ["python_workers"]

[build]
command = "pip install -r requirements.txt"

[[d1_databases]]
binding = "DB"
database_name = "neurobotanica-prod"
database_id = "your-database-id"

[[kv_namespaces]]
binding = "CACHE"
id = "your-kv-namespace-id"
preview_id = "your-kv-preview-id"
EOF
```

### Step 2: Database Setup (Cloudflare D1)
```bash
# Create production database
wrangler d1 create neurobotanica-prod

# Get database ID and update wrangler.toml
wrangler d1 list

# Migrate database schema
wrangler d1 execute neurobotanica-prod --file=../scripts/init_database.sql

# Load clinical data
wrangler d1 execute neurobotanica-prod --file=../data/clinical_evidence/init_data.sql
```

### Step 3: Environment Variables
```bash
# Set production environment variables
wrangler secret put API_KEY
wrangler secret put DATABASE_URL
wrangler secret put CLOUDFLARE_ACCOUNT_ID
wrangler secret put CLOUDFLARE_API_TOKEN
```

### Step 4: Deploy Backend
```bash
# Deploy to production
wrangler deploy

# Get deployment URL
wrangler tail --format=pretty
```

---

## 🎨 Frontend Deployment (Cloudflare Pages)

### Prerequisites
```bash
# Install dependencies
cd frontend
npm install

# Install Wrangler for Pages
npm install -g wrangler
```

### Step 1: Configure Next.js for Static Export
```javascript
// next.config.js
module.exports = {
  output: 'export',
  trailingSlash: true,
  images: {
    unoptimized: true,
  },
  env: {
    API_BASE_URL: process.env.NODE_ENV === 'production'
      ? 'https://neurobotanica-api.your-domain.workers.dev'
      : 'http://localhost:8000'
  }
}
```

### Step 2: Build Static Site
```bash
# Build for production
npm run build:static

# Verify build output
ls -la out/
```

### Step 3: Deploy to Cloudflare Pages
```bash
# Create Pages project
wrangler pages project create neurobotanica-pos

# Deploy static files
wrangler pages deployment create out/ --project-name=neurobotanica-pos

# Set environment variables
wrangler pages secret put API_BASE_URL
```

### Step 4: Configure Custom Domain (Optional)
```bash
# Add custom domain
wrangler pages domain add your-domain.com --project-name=neurobotanica-pos

# Configure DNS (point to Cloudflare)
# Add CNAME record: pos.your-domain.com -> neurobotanica-pos.pages.dev
```

---

## 🎭 Demo Deployment (Cloudflare Pages Dev)

### Create Demo Environment
```bash
# Create demo project
wrangler pages project create neurobotanica-demo

# Deploy demo build
npm run build:static
wrangler pages deployment create out/ --project-name=neurobotanica-demo --environment=demo

# Demo URL will be: https://neurobotanica-demo.pages.dev
```

### Demo Configuration
```javascript
// Demo-specific config
const DEMO_CONFIG = {
  apiBaseUrl: 'https://neurobotanica-demo-api.your-domain.workers.dev',
  features: {
    fullTradeSecrets: false, // Limited demo access
    sampleDataOnly: true,
    watermark: 'DEMO VERSION'
  }
}
```

---

## 🔧 Deployment Scripts

### Automated Deployment Script
```bash
#!/bin/bash
# deploy.sh - Complete deployment automation

echo "🚀 Starting NeuroBotanica Deployment..."

# Backend deployment
echo "📦 Deploying Backend..."
cd backend
wrangler deploy

# Frontend deployment
echo "🎨 Deploying Frontend..."
cd ../frontend
npm run build:static
wrangler pages deployment create out/

# Demo deployment
echo "🎭 Deploying Demo..."
wrangler pages deployment create out/ --project-name=neurobotanica-demo

echo "✅ Deployment Complete!"
echo "🌐 Production: https://neurobotanica-pos.pages.dev"
echo "🎭 Demo: https://neurobotanica-demo.pages.dev"
```

### Rollback Script
```bash
#!/bin/bash
# rollback.sh - Emergency rollback

echo "🔄 Rolling back deployment..."

# Rollback backend
wrangler deployments list
wrangler deployments rollback <previous-deployment-id>

# Rollback frontend
wrangler pages deployment list --project-name=neurobotanica-pos
wrangler pages deployment rollback <previous-deployment-id>

echo "✅ Rollback Complete!"
```

---

## 🔒 Security Configuration

### API Security
```python
# backend/middleware/security.py
from fastapi import Request, HTTPException
import os

API_KEYS = os.getenv("API_KEYS", "").split(",")

async def validate_api_key(request: Request):
    api_key = request.headers.get("X-API-Key")
    if api_key not in API_KEYS:
        raise HTTPException(status_code=401, detail="Invalid API Key")
    return True
```

### CORS Configuration
```python
# backend/main.py
app.add_middleware(
    CORSMiddleware,
    allow_origins=[
        "https://neurobotanica-pos.pages.dev",
        "https://neurobotanica-demo.pages.dev",
        "http://localhost:3000"  # Development
    ],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)
```

---

## 📊 Monitoring & Analytics

### Cloudflare Analytics
```bash
# Enable analytics
wrangler tail --format=json > logs/production.log

# Monitor performance
curl -H "Authorization: Bearer $CLOUDFLARE_API_TOKEN" \
  "https://api.cloudflare.com/client/v4/accounts/$ACCOUNT_ID/workers/scripts/neurobotanica-api/stats"
```

### Error Tracking
```python
# backend/middleware/error_tracking.py
import logging
from fastapi import Request

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

async def log_errors(request: Request, exc: Exception):
    logger.error(f"Error: {exc} | URL: {request.url} | Method: {request.method}")
    # Send to error tracking service
```

---

## 🧪 Testing Deployment

### Pre-Deployment Tests
```bash
# Run integration tests
python -m pytest tests/integration/ -v

# Test API endpoints
curl -H "X-API-Key: $API_KEY" https://neurobotanica-api.your-domain.workers.dev/api/health

# Test frontend build
cd frontend && npm run build && npm run export
```

### Post-Deployment Validation
```bash
# Health check
curl https://neurobotanica-pos.pages.dev/api/health

# API functionality test
curl -X POST https://neurobotanica-api.your-domain.workers.dev/api/dispensary/recommend \
  -H "Content-Type: application/json" \
  -H "X-API-Key: $API_KEY" \
  -d '{"condition": "chronic_pain", "experience_level": "intermediate"}'
```

---

## 🚦 Traffic Management

### Blue-Green Deployment
```bash
# Deploy to staging first
wrangler deploy --env staging

# Test staging
curl https://neurobotanica-api-staging.your-domain.workers.dev/api/health

# Promote to production
wrangler deployments promote <staging-deployment-id>
```

### Load Balancing
```bash
# Cloudflare automatically handles load balancing
# Configure traffic splitting for A/B testing
wrangler pages deployment tail --project-name=neurobotanica-pos --percentage=10
```

---

## 💰 Cost Optimization

### Cloudflare Pricing (Free Tier Maximization)
- **Workers:** 100,000 requests/day free
- **D1 Database:** 500MB free, $0.001/GB after
- **KV Storage:** 30GB free, $0.001/GB after
- **Pages:** Unlimited static sites free

### Cost Monitoring
```bash
# Monitor usage
wrangler tail --format=json | jq '.event' | sort | uniq -c

# Set up billing alerts
# Cloudflare Dashboard > Billing > Alerts
```

---

## 🔄 CI/CD Pipeline

### GitHub Actions Workflow
```yaml
# .github/workflows/deploy.yml
name: Deploy NeuroBotanica
on:
  push:
    branches: [ main ]
  pull_request:
    branches: [ main ]

jobs:
  test:
    needs: test
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v3
      - name: Run tests
        run: |
          python -m pytest tests/ -v

  deploy-backend:
    needs: test
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v3
      - name: Deploy Backend
        run: |
          cd backend
          wrangler deploy

  deploy-frontend:
    needs: test
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v3
      - name: Deploy Frontend
        run: |
          cd frontend
          npm install
          npm run build:static
          wrangler pages deployment create out/
```

---

## 🎯 Deployment Checklist

### Pre-Deployment
- [ ] All tests passing (30/30 integration tests)
- [ ] Database schema migrated
- [ ] Environment variables configured
- [ ] API keys generated and distributed
- [ ] Custom domain configured (optional)
- [ ] SSL certificates valid

### Deployment
- [ ] Backend deployed to Cloudflare Workers
- [ ] Frontend deployed to Cloudflare Pages
- [ ] Demo environment deployed
- [ ] DNS records updated
- [ ] CDN cache purged

### Post-Deployment
- [ ] Health checks passing
- [ ] API endpoints responding
- [ ] Frontend loading correctly
- [ ] Demo environment accessible
- [ ] Monitoring alerts configured
- [ ] Backup procedures tested

---

## 🚨 Emergency Procedures

### Service Outage
1. Check Cloudflare status: https://www.cloudflarestatus.com/
2. Review error logs: `wrangler tail`
3. Rollback to previous deployment
4. Notify stakeholders

### Data Loss Recovery
1. Restore from D1 database backup
2. Verify data integrity
3. Update frontend cache
4. Test all functionality

---

## 📞 Support & Maintenance

### Regular Maintenance
- **Weekly:** Review error logs and performance metrics
- **Monthly:** Update dependencies and security patches
- **Quarterly:** Review and optimize costs

### Support Contacts
- **Technical Support:** dev@cloakandquill.org
- **Business Support:** admin@cloakandquill.org
- **Emergency:** +1 (702) 555-0123

---

## 🎉 Success Metrics

### Deployment Success Criteria
- ✅ **Uptime:** 99.9% availability
- ✅ **Performance:** <200ms API response time
- ✅ **Security:** Zero data breaches
- ✅ **User Adoption:** 5 dispensaries onboarded in Month 1
- ✅ **Revenue:** $11,000 MRR achieved

### Monitoring Dashboard
- **Real-time Metrics:** Cloudflare Analytics
- **Error Tracking:** Custom logging
- **Performance:** Response time monitoring
- **Business:** Revenue and user metrics

---

*Powered by Cloak and Quill Research 501(c)(3) | Henderson, Nevada | Patent Portfolio: $684M-$1.026B Value*
# render.yaml
services:
  - type: web
    name: neurobotanica-api
    runtime: python3
    buildCommand: pip install -r requirements.txt
    startCommand: uvicorn backend.main:app --host 0.0.0.0 --port $PORT
    envVars:
      - key: DATABASE_URL
        value: your_db_url
      - key: API_ENV
        value: production
```

---

## 🎨 Frontend Deployment (Next.js → Cloudflare Pages)

### Prerequisites
```bash
# Install Cloudflare CLI
npm install -g wrangler

# Login to Cloudflare
wrangler auth login
```

### 1. Create Cloudflare Pages Project
```bash
cd frontend

# Create pages project
wrangler pages project create neurobotanica-pos

# Configure build settings
wrangler pages deployment create . \
  --project-name neurobotanica-pos \
  --branch main \
  --commit-hash $(git rev-parse HEAD) \
  --compatibility-date 2023-12-01
```

### 2. Configure Build Settings
Create `wrangler.toml` in frontend directory:
```toml
name = "neurobotanica-pos"
compatibility_date = "2023-12-01"
pages_build_output_dir = "out"

[build]
command = "npm run build:static"
cwd = "."

[build.environment_variables]
NODE_ENV = "production"
API_URL = "https://neurobotanica-api.your-domain.workers.dev"

[[pages_build_config]]
# Build configuration for Pages
build_command = "npm run build:static"
destination_dir = "out"
root_dir = "."
```

### 3. Deploy Frontend
```bash
# Build and deploy
npm run build:static
wrangler pages deploy out --project-name neurobotanica-pos

# Get deployment URL
# https://neurobotanica-pos.pages.dev
```

### 4. Custom Domain (Optional)
```bash
# Add custom domain
wrangler pages domain add neurobotanica-pos your-domain.com

# SSL certificate is automatic
```

---

## 🎭 Demo Deployments on Cloudflare Pages

### Strategy: Multiple Demo Environments

#### 1. Demo Configuration
Create `frontend/demo-config.js`:
```javascript
// Demo environment configuration
export const DEMO_CONFIG = {
  // Use mock data instead of live API
  useMockData: true,

  // Demo-specific settings
  demoMode: true,
  showDebugInfo: true,

  // Mock API responses
  mockResponses: {
    recommendations: [
      {
        rank: 1,
        product_name: "Blue Dream (Demo)",
        match_score: 0.87,
        cannabinoid_profile: { THC: 18, CBD: 1, CBG: 1.2 },
        key_terpenes: [
          { name: "myrcene", percent: 0.8, effects: ["sedative"] }
        ],
        why_recommended: "High THC for pain relief, myrcene for muscle relaxation",
        expected_benefits: ["pain reduction", "muscle relaxation"],
        dosage_guidance: {
          starting_dose: "5-10mg THC",
          wait_time_minutes: 30,
          max_dose_per_session: "20mg THC"
        }
      }
    ]
  }
};
```

#### 2. Demo-Specific Build
Update `package.json` scripts:
```json
{
  "scripts": {
    "build:demo": "NODE_ENV=demo next build && next export",
    "deploy:demo": "npm run build:demo && wrangler pages deploy out --project-name neurobotanica-demo"
  }
}
```

#### 3. Create Demo Pages Project
```bash
# Create separate demo project
wrangler pages project create neurobotanica-demo

# Deploy demo version
npm run deploy:demo
```

#### 4. Preview Deployments
```bash
# Create preview deployment for each PR
wrangler pages deployment create . \
  --project-name neurobotanica-pos \
  --branch feature/new-feature \
  --commit-hash $(git rev-parse HEAD)

# Preview URL: https://abc123.neurobotanica-pos.pages.dev
```

### Demo Environments Available:
- **Production Demo:** `https://neurobotanica-pos.pages.dev`
- **Development Demo:** `https://neurobotanica-demo.pages.dev`
- **Feature Previews:** `https://[commit-hash].neurobotanica-pos.pages.dev`

---

## 🔧 Environment Configuration

### Backend Environment Variables
```bash
# Production
DATABASE_URL=postgresql://user:pass@host:5432/db
API_ENV=production
CLOUDFLARE_API_TOKEN=your_token
OPENAI_API_KEY=your_key  # For advanced features

# Development
DATABASE_URL=sqlite:///neurobotanica_dev.db
API_ENV=development
DEBUG=true
```

### Frontend Environment Variables
```bash
# .env.local (development)
API_URL=http://localhost:8000
NODE_ENV=development

# .env.production (production)
API_URL=https://neurobotanica-api.your-domain.workers.dev
NODE_ENV=production

# .env.demo (demo)
API_URL=https://demo-api.neurobotanica.com
NODE_ENV=demo
USE_MOCK_DATA=true
```

### Cloudflare Pages Environment Variables
```bash
# Set via wrangler
wrangler pages secret put API_URL --project-name neurobotanica-pos
wrangler pages secret put DATABASE_URL --project-name neurobotanica-pos
```

---

## 🔄 CI/CD Pipeline (GitHub Actions)

### 1. Backend CI/CD
Create `.github/workflows/backend-deploy.yml`:
```yaml
name: Backend Deploy
on:
  push:
    branches: [main]
    paths: ['backend/**', 'requirements.txt']

jobs:
  deploy:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v3
      - name: Setup Python
        uses: actions/setup-python@v4
        with:
          python-version: '3.13'
      - name: Install dependencies
        run: pip install -r requirements.txt
      - name: Run tests
        run: python -m pytest tests/ -v
      - name: Deploy to Cloudflare Workers
        run: |
          npm install -g wrangler
          wrangler deploy
        env:
          CLOUDFLARE_API_TOKEN: ${{ secrets.CLOUDFLARE_API_TOKEN }}
```

### 2. Frontend CI/CD
Create `.github/workflows/frontend-deploy.yml`:
```yaml
name: Frontend Deploy
on:
  push:
    branches: [main]
    paths: ['frontend/**']

jobs:
  deploy:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v3
      - name: Setup Node
        uses: actions/setup-node@v3
        with:
          node-version: '18'
      - name: Install dependencies
        run: cd frontend && npm ci
      - name: Run tests
        run: cd frontend && npm run lint
      - name: Build
        run: cd frontend && npm run build:static
      - name: Deploy to Cloudflare Pages
        run: |
          cd frontend
          npm install -g wrangler
          wrangler pages deploy out --project-name neurobotanica-pos
        env:
          CLOUDFLARE_API_TOKEN: ${{ secrets.CLOUDFLARE_API_TOKEN }}
```

### 3. Demo Deployment
Create `.github/workflows/demo-deploy.yml`:
```yaml
name: Demo Deploy
on:
  push:
    branches: [develop]
    paths: ['frontend/**']

jobs:
  demo:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v3
      - name: Deploy Demo
        run: |
          cd frontend
          npm ci
          npm run deploy:demo
        env:
          CLOUDFLARE_API_TOKEN: ${{ secrets.CLOUDFLARE_API_TOKEN }}
```

---

## 📊 Database Setup (Cloudflare D1)

### 1. Create D1 Database
```bash
# Create database
wrangler d1 create neurobotanica-db

# Get database ID and update wrangler.toml
wrangler d1 list
```

### 2. Database Schema Migration
```sql
-- Create tables in D1
CREATE TABLE clinical_studies (
  id INTEGER PRIMARY KEY,
  condition TEXT,
  cannabinoid TEXT,
  efficacy_score REAL,
  study_reference TEXT,
  bias_corrected BOOLEAN
);

CREATE TABLE patient_profiles (
  id TEXT PRIMARY KEY,
  experience_level TEXT,
  conditions TEXT,
  preferences TEXT,
  created_at DATETIME DEFAULT CURRENT_TIMESTAMP
);

CREATE TABLE recommendations (
  id TEXT PRIMARY KEY,
  patient_id TEXT,
  recommendation_data TEXT,
  timestamp DATETIME DEFAULT CURRENT_TIMESTAMP,
  feedback_score INTEGER
);
```

### 3. Seed Database
```bash
# Import clinical data
wrangler d1 execute neurobotanica-db --file=data/clinical_studies.sql

# Import ML models to KV storage
wrangler kv:key put --binding MODELS therapeutic_model @models/therapeutic_prediction_v1.joblib
```

---

## 🧪 Testing Deployments

### 1. End-to-End Testing
```bash
# Test production API
curl -X POST https://neurobotanica-api.workers.dev/api/dispensary/recommend \
  -H "Content-Type: application/json" \
  -d '{"condition": "chronic_pain", "experience_level": "intermediate"}'

# Test frontend deployment
curl https://neurobotanica-pos.pages.dev
```

### 2. Performance Monitoring
```bash
# Cloudflare Analytics
# Monitor response times, error rates, and usage

# Custom monitoring (add to backend)
@app.middleware("http")
async def log_requests(request: Request, call_next):
    start_time = time.time()
    response = await call_next(request)
    process_time = time.time() - start_time
    logger.info(f"Request: {request.url} - Time: {process_time:.3f}s")
    return response
```

---

## 🚀 Deployment Checklist

### Pre-Deployment
- [ ] Domain purchased and configured
- [ ] SSL certificates (automatic with Cloudflare)
- [ ] Environment variables set
- [ ] Database migrated and seeded
- [ ] ML models uploaded to KV storage

### Backend Deployment
- [ ] Cloudflare Workers project created
- [ ] API endpoints deployed and tested
- [ ] CORS configured for frontend domain
- [ ] Rate limiting implemented
- [ ] Error handling and logging configured

### Frontend Deployment
- [ ] Cloudflare Pages project created
- [ ] Build configuration optimized
- [ ] Environment variables configured
- [ ] Custom domain connected (optional)
- [ ] PWA functionality tested

### Demo Environments
- [ ] Demo project created
- [ ] Mock data configured
- [ ] Preview deployments working
- [ ] Demo URLs documented

### Post-Deployment
- [ ] End-to-end testing completed
- [ ] Performance monitoring active
- [ ] Backup procedures documented
- [ ] Rollback procedures tested

---

## 🔗 URLs and Endpoints

### Production
- **Frontend:** `https://neurobotanica-pos.pages.dev`
- **API:** `https://neurobotanica-api.your-domain.workers.dev`
- **Docs:** `https://docs.neurobotanica.com`

### Demo Environments
- **Main Demo:** `https://neurobotanica-demo.pages.dev`
- **Feature Previews:** `https://[commit-hash].neurobotanica-pos.pages.dev`

### Development
- **Local Frontend:** `http://localhost:3000`
- **Local API:** `http://localhost:8000`

---

## 🎯 Go-Live Sequence

1. **Week 1:** Deploy to staging environment, full testing
2. **Week 2:** Deploy demos, gather feedback from pilot dispensaries
3. **Week 3:** Production deployment, monitor performance
4. **Week 4:** Scale to additional dispensaries, achieve $11,000 MRR target

**Status: DEPLOYMENT READY** 🚀

---

*Powered by Cloak and Quill Research 501(c)(3) | Henderson, Nevada | Patent Portfolio: $684M-$1.026B Value*</content>
<parameter name="filePath">c:\Users\Dr. Contessa Petrini\OneDrive\Desktop\neurobotanica_project\docs\NeuroBotanica_Deployment_Guide.md