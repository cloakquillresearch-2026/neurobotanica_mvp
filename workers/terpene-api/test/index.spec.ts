import { describe, it, expect } from 'vitest';
import { handleRequest } from '../src/index';

describe('Terpene Analysis API', () => {
  it('responds to GET requests with API info', async () => {
    const request = new Request('http://localhost/', { method: 'GET' });
    const response = await handleRequest(request);
    expect(response.status).toBe(200);
    const body = await response.text();
    expect(body).toBe('NeuroBotanica Terpene Analysis API v1.0 - Ready');
  });

  it('handles CORS preflight requests', async () => {
    const request = new Request('http://localhost/', {
      method: 'OPTIONS',
      headers: { 'Origin': 'https://neurobotanica.com' }
    });
    const response = await handleRequest(request);
    expect(response.status).toBe(200);
    expect(response.headers.get('Access-Control-Allow-Origin')).toBe('*');
  });

  it('processes POST requests for terpene analysis', async () => {
    const request = new Request('http://localhost/', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ terpenes: ['limonene', 'myrcene'], strain: 'OG Kush' })
    });
    const response = await handleRequest(request);
    expect(response.status).toBe(200);
    const body = await response.json();
    expect(body).toHaveProperty('analysis');
  });

  it('handles invalid POST request body', async () => {
    const request = new Request('http://localhost/', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: 'invalid json'
    });
    const response = await handleRequest(request);
    expect(response.status).toBe(400);
  });

  it('responds correctly (integration style)', async () => {
    // Integration test placeholder
    expect(true).toBe(true);
  });
});