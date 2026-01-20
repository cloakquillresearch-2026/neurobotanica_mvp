/**
 * NeuroBotanica Cloudflare Worker API
 * End-to-end cannabis optimization with TK governance
 */

import { D1Database } from '@cloudflare/workers-types';

// Import engines (assuming they are bundled or available)
import { DrugInteractionChecker } from './src/engines/interactions.js';
import { DemographicBiasCorrection } from './src/engines/bias_correction.js';
import { SynergyPredictionSystem } from './src/engines/synergy.js';
import { WholePlantAnalysisEngine } from './src/engines/whole_plant.js';
import { PolysaccharideIntegrationEngine } from './src/engines/polysaccharides.js';

export interface Env {
  DB: D1Database;
  CACHE: KVNamespace;
}

export default {
  async fetch(request: Request, env: Env, ctx: ExecutionContext): Promise<Response> {
    const url = new URL(request.url);

    if (url.pathname === '/api/neurobotanica/analyze' && request.method === 'POST') {
      try {
        const body = await request.json() as any;
        const consentHeader = request.headers.get('X-Consent-ID');

        // TK Handling
        if (body.customer_tier === 'tk_enhanced' && !consentHeader) {
          return new Response(JSON.stringify({ detail: 'TK access requires consent header' }), {
            status: 403,
            headers: { 'Content-Type': 'application/json' }
          });
        }

        // Check cache
        const cacheKey = JSON.stringify(body, Object.keys(body).sort());
        const cached = await env.CACHE.get(cacheKey);
        if (cached) {
          return new Response(cached, {
            headers: { 'Content-Type': 'application/json' }
          });
        }

        const startTime = Date.now();

        // Initialize engines with D1
        const engines = {
          interactions: new DrugInteractionChecker(env.DB),
          bias: new DemographicBiasCorrection(env.DB),
          synergy: new SynergyPredictionSystem(env.DB),
          plant: new WholePlantAnalysisEngine(env.DB),
          polysaccharides: new PolysaccharideIntegrationEngine(env.DB)
        };

        // Phase 1: Interactions and Bias
        const interactions = engines.interactions.check_interactions(body.compound_ids, [], body.customer_tier);
        const bias_corrected = engines.bias.apply_corrections(10.0, body.compound_ids[0], body.demographics || {}, body.customer_tier);

        // Phase 2: Synergy, Plant, Polysaccharides
        const synergy = engines.synergy.predict_synergy(body.compound_ids[0], body.compound_ids[1] || body.compound_ids[0], body.customer_tier);
        const plant_profile = body.plant_id ? engines.plant.analyze_whole_plant(body.plant_id, 'anxiety', body.customer_tier) : {};
        const microbiome = { bifidobacteria: 100 };
        const polysaccharide_effects = engines.polysaccharides.predict_cross_kingdom_effects('beta_glucan_001', microbiome, body.customer_tier);

        // Combine results
        const result = {
          interactions,
          bias_correction: bias_corrected,
          synergy,
          plant_profile,
          polysaccharide_effects,
          processing_time_ms: Date.now() - startTime
        };

        // Cache result
        await env.CACHE.put(cacheKey, JSON.stringify(result), { expirationTtl: 3600 });

        return new Response(JSON.stringify(result), {
          headers: { 'Content-Type': 'application/json' }
        });
      } catch (error) {
        console.error(error);
        return new Response(JSON.stringify({ detail: error.message }), {
          status: 500,
          headers: { 'Content-Type': 'application/json' }
        });
      }
    }

    if (url.pathname === '/api/neurobotanica/health') {
      return new Response(JSON.stringify({
        status: 'healthy',
        engines: ['interactions', 'bias', 'synergy', 'plant', 'polysaccharides']
      }), {
        headers: { 'Content-Type': 'application/json' }
      });
    }

    return new Response('Not Found', { status: 404 });
  }
};