# NeuroBotanica API Proxy Worker

Cloudflare Worker that keeps `/api/*` calls on the same hostname while serving the `platform-site` static assets.

## Usage
```bash
cd workers/api-proxy
npm install
npm run dev        # local testing
npm run deploy     # deploys to production account
```

### Required variables
- `API_BASE_URL`: Public FastAPI base (e.g. `https://nb-backend-production.up.railway.app`).
- `PLATFORM_PAGES_URL`: Cloudflare Pages URL (e.g. `https://neurobotanica-mvp.pages.dev`).

Set them using `npx wrangler secret put API_BASE_URL` and `npx wrangler secret put PLATFORM_PAGES_URL`.

### Routing
After deployment, attach a route such as `neurobotanica-mvp.pages.dev/*` (or the custom domain) so requests hit the Worker. All non `/api/*` paths stream the Pages response; `/api/*` paths are forwarded to FastAPI with CORS headers.
