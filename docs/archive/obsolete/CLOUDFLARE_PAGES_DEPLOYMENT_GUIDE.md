# NeuroBotanica Cloudflare Pages Deployment Guide

## Overview
This guide provides step-by-step instructions for deploying the complete NeuroBotanica MVP to Cloudflare Pages using Pages Functions for API endpoints.

## Prerequisites
- Cloudflare account with Pages enabled
- Node.js and npm installed locally
- Access to your Cloudflare dashboard

## Project Structure
```
frontend/
├── out/                    # Static build output
├── functions/             # Pages Functions (API endpoints)
│   ├── api-health.js
│   ├── api-dispensary-recommend.js
│   ├── api-chempath-analyze.js
│   ├── api-toxpath-analyze.js
│   ├── api-biopath-analyze.js
│   ├── api-genomepath-analyze.js
│   ├── api-regpath-analyze.js
│   └── api-metapath-analyze.js
├── wrangler.toml         # Cloudflare configuration
└── package.json
```

## Deployment Steps

### Step 1: Build the Static Frontend
```bash
cd frontend
npm install
npm run build:static
```

### Step 2: Verify Build Output
Ensure the `out/` directory contains your built files:
- `index.html`
- `_next/` directory with assets
- Other static assets

### Step 3: Deploy to Cloudflare Pages

#### Option A: Using Wrangler CLI (Recommended)
```bash
# Install Wrangler if not already installed
npm install -g wrangler

# Login to Cloudflare
wrangler auth login

# Deploy to Pages
wrangler pages deploy out --compatibility-date 2024-01-01
```

#### Option B: Manual Upload via Dashboard
1. Go to [Cloudflare Pages Dashboard](https://dash.cloudflare.com/pages)
2. Click "Create a project"
3. Choose "Upload assets" or connect your Git repository
4. Upload the contents of the `out/` directory
5. Configure build settings:
   - Build command: `npm run build:static`
   - Build output directory: `out`
   - Root directory: `frontend`

### Step 4: Configure Pages Functions
After deployment, you need to add the Pages Functions:

1. In your Cloudflare Pages project dashboard, go to "Functions"
2. Upload each function file from the `functions/` directory:
   - `api-health.js` → `/api/health`
   - `api-dispensary-recommend.js` → `/api/dispensary-recommend`
   - `api-chempath-analyze.js` → `/api/chempath-analyze`
   - `api-toxpath-analyze.js` → `/api/toxpath-analyze`
   - `api-biopath-analyze.js` → `/api/biopath-analyze`
   - `api-genomepath-analyze.js` → `/api/genomepath-analyze`
   - `api-regpath-analyze.js` → `/api/regpath-analyze`
   - `api-metapath-analyze.js` → `/api/metapath-analyze`

### Step 5: Configure Environment Variables
In your Pages project settings, add these environment variables:
- `NODE_ENV`: `production`

### Step 6: Set Up Custom Domain (Optional)
If you have a custom domain:
1. Go to "Custom domains" in your Pages project
2. Add your domain (e.g., `neurobotanica.yourdomain.com`)
3. Configure DNS records as instructed

## API Endpoints

After deployment, your API endpoints will be available at:
- Health Check: `https://your-project.pages.dev/api/health`
- Dispensary Recommendations: `https://your-project.pages.dev/api/dispensary-recommend`
- ChemPath Analysis: `https://your-project.pages.dev/api/chempath-analyze`
- ToxPath Analysis: `https://your-project.pages.dev/api/toxpath-analyze`
- BioPath Analysis: `https://your-project.pages.dev/api/biopath-analyze`
- GenomePath Analysis: `https://your-project.pages.dev/api/genomepath-analyze`
- RegPath Analysis: `https://your-project.pages.dev/api/regpath-analyze`
- MetaPath Analysis: `https://your-project.pages.dev/api/metapath-analyze`

## Testing Deployment

### Frontend Testing
1. Visit your deployed URL
2. Test the main interface
3. Verify static assets load correctly

### API Testing
Use curl or a REST client to test endpoints:

```bash
# Health check
curl "https://your-project.pages.dev/api/health"

# Dispensary recommendation
curl "https://your-project.pages.dev/api/dispensary-recommend?condition=anxiety&experience=beginner"

# ChemPath analysis
curl "https://your-project.pages.dev/api/chempath-analyze?compound=CBD"
```

## Troubleshooting

### Common Issues

1. **Functions not working**: Ensure functions are uploaded to the correct routes
2. **CORS errors**: Check that CORS headers are properly configured in functions
3. **Build failures**: Verify `npm run build:static` completes successfully
4. **Assets not loading**: Check that `out/` directory contains all necessary files

### Logs and Monitoring
- Check Cloudflare Pages deployment logs
- Monitor function execution in the Functions tab
- Use browser developer tools to check for console errors

## Security Considerations
- All API endpoints include CORS headers for web access
- Functions run in Cloudflare's edge network
- No sensitive data is stored in the functions
- Mock data is used for demonstration purposes

## Performance
- Static assets served via Cloudflare's global CDN
- Functions execute at the edge for low latency
- Automatic scaling based on request volume

## Next Steps
1. Replace mock data with real NeuroBotanica algorithms
2. Implement proper authentication
3. Add rate limiting and security measures
4. Set up monitoring and analytics
5. Configure backup and disaster recovery

## Support
For issues with this deployment:
1. Check Cloudflare Pages documentation
2. Review function logs in the dashboard
3. Test locally before deploying
4. Ensure all dependencies are properly installed