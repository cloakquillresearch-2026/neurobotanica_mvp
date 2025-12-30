# Manual Cloudflare Deployment Guide

## Since API Token Issues - Manual Deployment Steps

### Step 1: Prepare Your Files
✅ **Already Done**: Your static build is ready in `frontend/out/`

### Step 2: Cloudflare Dashboard Deployment

1. **Go to Cloudflare Dashboard**
   - Visit: https://dash.cloudflare.com/
   - Login with your account

2. **Navigate to Pages**
   - Click "Pages" in the left sidebar
   - Click "Create a project"

3. **Create NeuroBotanica POS Project**
   - Project name: `neurobotanica-pos`
   - Click "Direct upload"

4. **Upload Your Files**
   - Drag and drop all files from `frontend/out/` folder
   - Or click "Select files" and choose the `out` folder
   - Click "Deploy site"

5. **Set Environment Variables** (Optional)
   - Go to your Pages project settings
   - Add environment variables:
     ```
     API_BASE_URL = https://your-backend-api-url
     NODE_ENV = production
     ```

### Step 3: Custom Domain (Optional)
- In Pages project settings → "Custom domains"
- Add: `pos.neurobotanica.com` (or your preferred domain)

### Step 4: Backend Deployment Options

Since Cloudflare Workers don't support Python, here are your backend options:

#### Option A: Railway (Recommended)
1. Go to https://railway.app
2. Connect your GitHub repo
3. Railway will auto-detect Python/FastAPI
4. Deploy automatically
5. Get your API URL (e.g., `https://neurobotanica-api.up.railway.app`)

#### Option B: Render
1. Go to https://render.com
2. Create "Web Service" from Git
3. Select Python/Python3
4. Set build command: `pip install -r requirements.txt`
5. Set start command: `uvicorn backend.main:app --host 0.0.0.0 --port $PORT`

#### Option C: Heroku
1. Install Heroku CLI: `npm install -g heroku`
2. `heroku create neurobotanica-api`
3. `git push heroku main`

### Step 5: Connect Frontend to Backend

Update your frontend's API calls to point to your backend URL:

```javascript
// In frontend code, update API_BASE_URL to your backend URL
const API_BASE_URL = 'https://your-backend-url-here';
```

### Step 6: Test Deployment

1. **Frontend**: Visit your Pages URL
2. **API**: Test endpoints like `/api/health`
3. **Integration**: Try a dispensary recommendation

---

## Quick Commands Reference

```bash
# Check your build is ready
cd frontend && dir out

# Manual upload instructions above

# Backend deployment (choose one):
# Railway: railway.app (easiest)
# Render: render.com
# Heroku: heroku.com
```

---

## Troubleshooting

### API Token Issues
- Go to Cloudflare Dashboard → Profile → API Tokens
- Create new token with Pages permissions
- Set `CLOUDFLARE_API_TOKEN` environment variable

### Build Issues
- Ensure `npm run build:static` completes successfully
- Check `frontend/out/` contains `index.html` and `_next/` folder

### Backend Issues
- Test locally: `uvicorn backend.main:app --reload`
- Check requirements.txt has all dependencies
- Ensure database file exists or is configured

---

**Your MVP is ready - just needs manual deployment to Cloudflare!** 🚀