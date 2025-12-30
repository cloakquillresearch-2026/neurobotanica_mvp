# Quick Deployment Guide - NeuroBotanica MVP

## 🚀 One-Command Deployment

### Prerequisites (Run Once)
```bash
# Install Railway CLI (for backend)
npm install -g @railway/cli

# Install Wrangler CLI (for frontend)
npm install -g wrangler

# Authenticate with services
railway login
wrangler auth login

# Install Python dependencies
pip install -r requirements.txt
```

### Full Deployment (Production + Demo)
```bash
# Run the automated deployment script
./scripts/deploy.bat
```

This will:
1. ✅ Run all tests (30/30 integration tests)
2. 📦 Deploy FastAPI backend to Railway (Python-compatible)
3. 🎨 Deploy Next.js frontend to Cloudflare Pages
4. 🎭 Deploy demo environment to Cloudflare Pages
5. 🔍 Validate deployment

## 🌐 Access URLs

After deployment, your NeuroBotanica MVP will be available at:

### Production Environment
- **Frontend (Budtender POS)**: `https://neurobotanica-pos.pages.dev`
- **Backend API**: `https://neurobotanica-api.up.railway.app` (Railway URL)

### Demo Environment
- **Demo Frontend**: `https://neurobotanica-demo.pages.dev`
- **Demo API**: Same as production (with demo flag)

## 🎭 Demo Features

The demo environment includes:
- 🚫 **Watermark**: "DEMO VERSION" banner
- 📊 **Sample Data**: Pre-loaded with example conditions
- 🔒 **Limited Access**: 5 recommendations max per session
- 📧 **Contact Info**: Demo-specific contact information
- 🎯 **Educational Focus**: Showcases capabilities without full access

## 🧪 Testing Deployment

### Quick Health Check
```bash
# Test backend API
curl -H "X-API-Key: your-api-key" https://neurobotanica-api.up.railway.app/api/health

# Test frontend
curl https://neurobotanica-pos.pages.dev
```

### Functional Testing
```bash
# Test recommendation endpoint
curl -X POST https://neurobotanica-api.up.railway.app/api/dispensary/recommend \
  -H "Content-Type: application/json" \
  -H "X-API-Key: your-api-key" \
  -d '{
    "condition": "chronic_pain",
    "experience_level": "intermediate",
    "administration_method": "flower"
  }'
```

## 🔄 Rollback (If Needed)

If deployment fails or issues arise:
```bash
# Rollback Railway deployment
railway rollback

# Rollback Cloudflare Pages
./scripts/rollback.bat
```

## ⚙️ Configuration

### Environment Variables
Set these in Railway dashboard and Cloudflare:

**Railway (Backend)**:
```
API_KEY=your-production-api-key
SECRET_KEY=your-secret-key
DATABASE_URL=sqlite:///neurobotanica_prod.db
API_ENV=production
```

**Cloudflare Pages (Frontend)**:
```
API_BASE_URL=https://neurobotanica-api.up.railway.app
NODE_ENV=production
```

### Custom Domain (Optional)
```bash
# Add custom domain to Railway
# Go to Railway dashboard > Settings > Domains

# Add custom domain to Cloudflare Pages
wrangler pages domain add pos.neurobotanica.com --project-name=neurobotanica-pos
```

## 📊 Monitoring

### Real-time Logs
```bash
# Railway logs
railway logs

# Cloudflare Pages logs
wrangler pages deployment tail --project-name=neurobotanica-pos
```

### Performance Metrics
- **Railway Dashboard**: Real-time metrics and logs
- **Cloudflare Dashboard**: CDN analytics and performance
- **API Response Times**: Target <200ms globally
- **Uptime**: 99.9% SLA with Railway + Cloudflare

## 🚨 Emergency Contacts

- **Technical Issues**: dev@cloakandquill.org
- **Business Issues**: admin@cloakandquill.org
- **Railway Support**: https://docs.railway.app/
- **Cloudflare Support**: https://support.cloudflare.com/

## 🎯 Success Checklist

- [ ] ✅ Tests passing (30/30)
- [ ] 🌐 Frontend accessible
- [ ] 🔗 API responding
- [ ] 🎭 Demo working
- [ ] 📊 Analytics configured
- [ ] 🔒 Security enabled
- [ ] 📧 Dispensaries notified

## 💰 Cost Estimate

**Free Tier (First 6 months)**:
- Railway: $5/month credit, then $0 for hobby plan
- Cloudflare Pages: Unlimited static sites free
- Total: ~$5/month

**Growth Tier (After free limits)**:
- Railway Pro: $10/month
- Cloudflare: Minimal additional costs
- Estimated cost at scale: <$50/month

---

*Ready for Nevada dispensary deployment targeting $11,000 MRR*

*Powered by Cloak and Quill Research 501(c)(3)*

## 🧪 Testing Deployment

### Quick Health Check
```bash
# Test backend API
curl -H "X-API-Key: your-api-key" https://neurobotanica-api.your-domain.workers.dev/api/health

# Test frontend
curl https://neurobotanica-pos.pages.dev
```

### Functional Testing
```bash
# Test recommendation endpoint
curl -X POST https://neurobotanica-api.your-domain.workers.dev/api/dispensary/recommend \
  -H "Content-Type: application/json" \
  -H "X-API-Key: your-api-key" \
  -d '{
    "condition": "chronic_pain",
    "experience_level": "intermediate",
    "administration_method": "flower"
  }'
```

## 🔄 Rollback (If Needed)

If deployment fails or issues arise:
```bash
# Rollback to previous version
./scripts/rollback.sh
```

## ⚙️ Configuration

### Environment Variables
Set these in Cloudflare Dashboard or via Wrangler:

**Backend (Workers)**:
```
API_KEY=your-production-api-key
DATABASE_URL=your-d1-database-url
CLOUDFLARE_ACCOUNT_ID=your-account-id
```

**Frontend (Pages)**:
```
API_BASE_URL=https://neurobotanica-api.your-domain.workers.dev
NODE_ENV=production
```

### Custom Domain (Optional)
```bash
# Add custom domain to production
wrangler pages domain add pos.neurobotanica.com --project-name=neurobotanica-pos

# Add custom domain to demo
wrangler pages domain add demo.neurobotanica.com --project-name=neurobotanica-demo
```

## 📊 Monitoring

### Real-time Logs
```bash
# Backend logs
wrangler tail --format=pretty

# Frontend deployment logs
wrangler pages deployment tail --project-name=neurobotanica-pos
```

### Performance Metrics
- **Cloudflare Dashboard**: Real-time analytics
- **API Response Times**: Target <200ms globally
- **Uptime**: 99.9% SLA with Cloudflare

## 🚨 Emergency Contacts

- **Technical Issues**: dev@cloakandquill.org
- **Business Issues**: admin@cloakandquill.org
- **Cloudflare Support**: https://support.cloudflare.com/

## 🎯 Success Checklist

- [ ] ✅ Tests passing (30/30)
- [ ] 🌐 Frontend accessible
- [ ] 🔗 API responding
- [ ] 🎭 Demo working
- [ ] 📊 Analytics configured
- [ ] 🔒 Security enabled
- [ ] 📧 Dispensaries notified

## 💰 Cost Estimate

**Free Tier (First 6 months)**:
- Workers: 100,000 requests/day
- D1 Database: 500MB storage
- Pages: Unlimited static sites
- Total: $0/month

**Growth Tier (After free limits)**:
- Additional requests: $0.15/1,000
- Database storage: $0.001/GB
- Estimated cost at scale: <$50/month

---

*Ready for Nevada dispensary deployment targeting $11,000 MRR*

*Powered by Cloak and Quill Research 501(c)(3)*