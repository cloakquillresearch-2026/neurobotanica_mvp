param(
    [switch]$Preview
)

$root = Split-Path -Parent $PSScriptRoot
$platformPath = Join-Path $root "platform-site"

if (-not (Test-Path $platformPath)) {
    Write-Error "platform-site directory was not found at $platformPath"
    exit 1
}

Push-Location $platformPath
try {
    Write-Host "[platform-site] Installing dependencies..."
    npm install
    if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }

    Write-Host "[platform-site] Building Cloudflare Pages bundle via next-on-pages..."
    npm run cf:build
    if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }

    if ($Preview) {
        Write-Host "[platform-site] Starting local Pages preview (Ctrl+C to stop)..."
        npx wrangler pages dev .vercel/output/static
    }
    else {
        Write-Host "[platform-site] Deploying to Cloudflare Pages project neurobotanica-mvp..."
        npx wrangler pages deploy .vercel/output/static --project-name neurobotanica-mvp --branch main
    }
}
finally {
    Pop-Location
}
