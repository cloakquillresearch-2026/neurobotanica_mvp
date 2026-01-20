/**
 * NeuroBotanica API Worker
 * Cloudflare Worker for Budtender integration
 */

export default {
  async fetch(request, env) {
    const url = new URL(request.url);

    if (url.pathname === '/api/neurobotanica/health') {
      return new Response(JSON.stringify({
        status: 'healthy',
        engines: ['interactions', 'bias', 'synergy', 'plant', 'polysaccharides']
      }), {
        headers: { 'Content-Type': 'application/json' }
      });
    }

    if (url.pathname === '/api/neurobotanica/analyze' && request.method === 'POST') {
      try {
        const body = await request.json();
        const consentHeader = request.headers.get('X-Consent-ID');

        // TK Check
        if (body.customer_tier === 'tk_enhanced' && !consentHeader) {
          return new Response(JSON.stringify({ detail: 'TK access requires consent header' }), {
            status: 403,
            headers: { 'Content-Type': 'application/json' }
          });
        }

        // Mock analysis (replace with actual engine calls)
        const result = {
          interactions: { warnings: [], severity: 'low' },
          bias_correction: { adjusted_dose_mg: 10.0, factors_applied: [] },
          synergy: { synergy_score: 0.8, confidence: 0.9, tk_enhanced: body.customer_tier === 'tk_enhanced' },
          plant_profile: body.plant_id ? { ratios: { thc_cbd: 1.5 } } : {},
          polysaccharide_effects: { bifidobacteria_increase: 20 },
          processing_time_ms: Math.random() * 100
        };

        return new Response(JSON.stringify(result), {
          headers: { 'Content-Type': 'application/json' }
        });
      } catch (error) {
        return new Response(JSON.stringify({ detail: error.message }), {
          status: 500,
          headers: { 'Content-Type': 'application/json' }
        });
      }
    }

    return new Response('Not Found', { status: 404 });
  }
};