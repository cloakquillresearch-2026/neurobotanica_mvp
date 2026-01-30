$logPath = "$env:USERPROFILE\.wrangler\logs\wrangler-2026-01-30_06-07-43_806.log"
$destPath = ".\.wrangler_log_latest.log"
if (Test-Path $logPath) {
    Copy-Item $logPath $destPath -Force
    Write-Host "Log copied successfully."
} else {
    Write-Host "Log file not found."
}