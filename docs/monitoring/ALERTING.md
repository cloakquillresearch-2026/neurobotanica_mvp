# Monitoring & Alerting for TS-PS-001 (NeuroBotanica)

This document explains sample Prometheus alert rules and Datadog monitor templates you can use to detect regressions and incidents related to the TS-PS-001 Inflammatory Synergy Engine.

## Prometheus Alert Rules
Files:
- `monitoring/prometheus/rules/inflammatory_rules.yml`

Key alerts:
- NeuroBotanicaInflammatoryHighLatency
  - Trigger: 95th percentile latency (5m window) > 20ms
  - Action: Page on-call; check service logs and model loads
- NeuroBotanicaInflammatoryCacheMissRatioHigh
  - Trigger: Cache misses / total > 20% over 5m and > 10 total calls
  - Action: Check Redis connectivity, TTL, and key cardinality
- NeuroBotanicaHealthDegraded
  - Trigger: `neurobotanica_health_status > 0` in Prometheus
  - Action: Check /health endpoint details (db, ml models, coefficients)

> Thresholds are defaults for small deployments. Tune the numeric thresholds to your expected SLOs and latency baselines.

## Alertmanager Example
File:
- `monitoring/alertmanager/alertmanager-example.yml`

This is a minimal Alertmanager config that routes all alerts to an `oncall@example.org` email. Replace receivers with your PagerDuty, Slack, or other integrations and configure routes for maintenance windows and silence rules.

## Datadog Monitor Templates
Files:
- `monitoring/datadog/inflammatory_latency_monitor.json`
- `monitoring/datadog/inflammatory_cache_miss_monitor.json`

These are example monitor payloads you can import via Datadog UI or API. They demonstrate how to alert on the same conditions as Prometheus.

## How to deploy
1. Add `inflammatory_rules.yml` to your Prometheus or Thanos recording rules and reload Prometheus (or apply via the Prometheus Operator).
2. Configure Alertmanager with `alertmanager-example.yml` and connect your notification receivers.
3. Import Datadog monitor JSON files using the Datadog UI or `datadog` CLI/API.

## Runbook tips
- For high latency: check recent deployments, CPU/memory, model file presence, and health monitor outputs.
- For high cache misses: verify `TS_PS_001_REDIS_URL` connectivity, Redis logs, and TTL settings. If Redis is unavailable, the engine falls back to in-memory cache.

## Notes
- These are example templates. Update the alert severity, duration (`for`), and thresholds to match your SLOs.
- Consider adding automatic incident runbooks and escalation policies in your alerting platform.
