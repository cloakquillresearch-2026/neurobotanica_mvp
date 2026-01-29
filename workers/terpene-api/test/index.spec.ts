import { describe, it, expect } from 'vitest';
import { handleRequest } from '../src/index';

// Mock env for tests
const mockEnv = {
  NEUROBOTANICA_DB: {
    prepare: () => ({
      bind: () => ({
        all: async () => ({ results: [] })
      })
    })
  },
  NEUROBOTANICA_CACHE: null
};

describe('Terpene Analysis API', () => {
  it('responds to GET / with API info JSON', async () => {
    const request = new Request('http://localhost/', { method: 'GET' });
    const response = await handleRequest(request, mockEnv);
    expect(response.status).toBe(200);
    const body = await response.json();
    expect(body).toHaveProperty('name', 'NeuroBotanica Terpene Analysis API');
    expect(body).toHaveProperty('version', '1.0.0');
    expect(body).toHaveProperty('endpoints');
  });

  it('responds to GET /health with health check', async () => {
    const request = new Request('http://localhost/health', { method: 'GET' });
    const response = await handleRequest(request, mockEnv);
    expect(response.status).toBe(200);
    const body = await response.json();
    expect(body).toHaveProperty('status', 'healthy');
    expect(body).toHaveProperty('engines');
  });

  it('responds to GET /api/neurobotanica/health with health check', async () => {
    const request = new Request('http://localhost/api/neurobotanica/health', { method: 'GET' });
    const response = await handleRequest(request, mockEnv);
    expect(response.status).toBe(200);
    const body = await response.json();
    expect(body).toHaveProperty('status', 'healthy');
  });

  it('handles CORS preflight requests', async () => {
    const request = new Request('http://localhost/', {
      method: 'OPTIONS',
      headers: { 'Origin': 'https://neurobotanica.com' }
    });
    const response = await handleRequest(request, mockEnv);
    expect(response.status).toBe(200);
    expect(response.headers.get('Access-Control-Allow-Origin')).toBe('*');
  });

  it('processes POST /analyze for terpene analysis', async () => {
    const request = new Request('http://localhost/analyze', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ compound_ids: ['cbd', 'thc'], strain: 'OG Kush' })
    });
    const response = await handleRequest(request, mockEnv);
    expect(response.status).toBe(200);
    const body = await response.json();
    expect(body).toHaveProperty('interactions');
    expect(body).toHaveProperty('bias_correction');
    expect(body).toHaveProperty('synergy');
  });

  it('processes POST /api/neurobotanica/analyze for full analysis', async () => {
    const request = new Request('http://localhost/api/neurobotanica/analyze', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ compound_ids: ['cbd'], demographics: { age: 35 } })
    });
    const response = await handleRequest(request, mockEnv);
    expect(response.status).toBe(200);
    const body = await response.json();
    expect(body).toHaveProperty('processing_time_ms');
  });

  it('rejects POST /analyze with missing compound_ids', async () => {
    const request = new Request('http://localhost/analyze', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ strain: 'OG Kush' })
    });
    const response = await handleRequest(request, mockEnv);
    expect(response.status).toBe(400);
    const body = await response.json();
    expect(body).toHaveProperty('error');
  });

  it('handles invalid POST request body', async () => {
    const request = new Request('http://localhost/analyze', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: 'invalid json'
    });
    const response = await handleRequest(request, mockEnv);
    expect(response.status).toBe(500);
  });

  it('returns 404 for unknown paths', async () => {
    const request = new Request('http://localhost/unknown', { method: 'GET' });
    const response = await handleRequest(request, mockEnv);
    expect(response.status).toBe(404);
    const body = await response.json();
    expect(body).toHaveProperty('error', 'Not Found');
  });

  it('processes POST /api/neurobotanica/bias-correction', async () => {
    const request = new Request('http://localhost/api/neurobotanica/bias-correction', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ compound_id: 'cbd', base_dose: 10, demographics: { age: 45 } })
    });
    const response = await handleRequest(request, mockEnv);
    expect(response.status).toBe(200);
    const body = await response.json();
    expect(body).toHaveProperty('adjusted_dose_mg');
  });
});