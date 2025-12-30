# Manual Deployment Guide - NeuroBotanica MVP

## Step-by-Step Deployment (No Automation)

If you prefer manual control over the automated script, follow these steps:

---

## Step 1: Backend Deployment (Railway)

### 1.1 Create Railway Account
1. Go to https://railway.app
2. Sign up/Sign in with GitHub
3. Verify your email

### 1.2 Deploy Backend
```bash
# Navigate to backend directory
cd backend

# Login to Railway
railway login

# Initialize project
railway init neurobotanica-api

# Set environment variables
railway variables set API_ENV=production
railway variables set SECRET_KEY=your-super-secret-key-here
railway variables set API_KEY=your-production-api-key

# Deploy
railway up
```

### 1.3 Get Backend URL
- Go to Railway Dashboard
- Select your project
- Copy the deployment URL (e.g., `https://neurobotanica-api.up.railway.app`)

---

## Step 2: Frontend Deployment (Cloudflare Pages)

### 2.1 Prepare Frontend
```bash
# Navigate to frontend directory
cd ../frontend

# Install dependencies
npm install

# Update API URL in next.config.js
# Replace the API_BASE_URL with your Railway URL
```

### 2.2 Deploy to Cloudflare Pages
```bash
# Login to Cloudflare
npx wrangler auth login

# Build static site
npm run build:static

# Create Pages project
npx wrangler pages project create neurobotanica-pos

# Deploy
npx wrangler pages deployment create out/ --project-name=neurobotanica-pos
```

---

## Step 3: Demo Deployment

### 3.1 Create Demo Build
```bash
# In frontend directory
npm run build:demo  # If you set up demo config

# Or just use regular build for demo
npm run build:static
```

### 3.2 Deploy Demo
```bash
# Create demo project
npx wrangler pages project create neurobotanica-demo

# Deploy demo
npx wrangler pages deployment create out/ --project-name=neurobotanica-demo
```

---

## Step 4: Testing

### 4.1 Test Backend
```bash
# Replace with your actual Railway URL
curl https://your-railway-url.up.railway.app/api/health
```

### 4.2 Test Frontend
- Open `https://neurobotanica-pos.pages.dev`
- Try a sample recommendation request

### 4.3 Test Demo
- Open `https://neurobotanica-demo.pages.dev`
- Verify demo watermark and limitations

---

## Step 5: Configuration

### Environment Variables Checklist
- [ ] Railway: `API_KEY` set
- [ ] Railway: `SECRET_KEY` set
- [ ] Railway: `API_ENV=production`
- [ ] Frontend: API_BASE_URL updated in next.config.js

### Custom Domain (Optional)
1. **Railway**: Dashboard → Settings → Domains → Add domain
2. **Cloudflare Pages**: `wrangler pages domain add yourdomain.com`

---

## Troubleshooting

### Backend Issues
- Check Railway logs: `railway logs`
- Verify environment variables in Railway dashboard
- Ensure requirements.txt is in backend directory

### Frontend Issues
- Check build output in `out/` directory
- Verify API_BASE_URL in next.config.js
- Check Cloudflare Pages deployment logs

### Demo Issues
- Ensure demo config is set up
- Check if demo build script exists
- Verify demo project creation

---

## Quick Commands Reference

```bash
# Backend
cd backend
railway login
railway init neurobotanica-api
railway variables set API_KEY=your-key
railway up

# Frontend
cd ../frontend
npm install
npm run build:static
npx wrangler auth login
npx wrangler pages project create neurobotanica-pos
npx wrangler pages deployment create out/ --project-name=neurobotanica-pos

# Demo
npx wrangler pages project create neurobotanica-demo
npx wrangler pages deployment create out/ --project-name=neurobotanica-demo
```

---

*Manual deployment gives you full control over each step*