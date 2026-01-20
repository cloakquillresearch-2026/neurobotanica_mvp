/**
 * NeuroBotanica Cloudflare Worker
 * Edge-deployed API for cannabis optimization with TK preservation
 * Performance Target: <200ms global responses
 */

import { D1Database, KVNamespace } from '@cloudflare/workers-types';

interface Env {
  NEUROBOTANICA_DB: D1Database;
  NEUROBOTANICA_CACHE: KVNamespace;
}

interface AnalysisRequest {
  compound_ids: string[];
  demographics?: Record<string, any>;
  customer_tier?: string;
  plant_id?: string;
}

interface EngineResult {
  interactions?: any;
  bias_correction?: any;
  synergy?: any;
  plant_profile?: any;
  polysaccharide_effects?: any;
  processing_time_ms: number;
}

export default {
  async fetch(request: Request, env: Env): Promise<Response> {
    const url = new URL(request.url);
    
    // Health check endpoint
    if (url.pathname === '/api/neurobotanica/health') {
      return new Response(JSON.stringify({
        status: 'healthy',
        engines: ['interactions', 'bias', 'synergy', 'plant', 'polysaccharides']
      }), {
        headers: { 'Content-Type': 'application/json' }
      });
    }
    
    // Analyze endpoint
    if (url.pathname === '/api/neurobotanica/analyze' && request.method === 'POST') {
      try {
        const reqBody: AnalysisRequest = await request.json();
        const startTime = Date.now();
        
        // Phase 1: Interactions and Bias
        const interactions = await checkInteractions(env, reqBody.compound_ids, reqBody.customer_tier || 'computational_only');
        const biasCorrected = await applyBiasCorrection(env, 10.0, reqBody.compound_ids[0], reqBody.demographics || {}, reqBody.customer_tier || 'computational_only');
        
        // Phase 2: Synergy, Plant, Polysaccharides
        const synergy = await predictSynergy(env, reqBody.compound_ids[0], reqBody.compound_ids[1] || reqBody.compound_ids[0], reqBody.customer_tier || 'computational_only');
        const plantProfile = reqBody.plant_id ? await analyzePlant(env, reqBody.plant_id, reqBody.customer_tier || 'computational_only') : {};
        const polysaccharideEffects = await predictPolysaccharides(env, reqBody.compound_ids, reqBody.customer_tier || 'computational_only');
        
        const result: EngineResult = {
          interactions,
          bias_correction: biasCorrected,
          synergy,
          plant_profile: plantProfile,
          polysaccharide_effects: polysaccharideEffects,
          processing_time_ms: Date.now() - startTime
        };
        
        return new Response(JSON.stringify(result), {
          headers: { 'Content-Type': 'application/json' }
        });
      } catch (error) {
        console.error('Analysis error:', error);
        return new Response(JSON.stringify({ error: 'Analysis failed' }), {
          status: 500,
          headers: { 'Content-Type': 'application/json' }
        });
      }
    }
    
    return new Response('Not Found', { status: 404 });
  }
};

// Helper functions for engines (simplified for edge deployment)
async function checkInteractions(env: Env, compoundIds: string[], tier: string): Promise<any> {
  // TK consent check
  if (tier === 'tk_enhanced') {
    const consentValid = await verifyConsent(env, 'sample_consent_id');
    if (!consentValid) {
      throw new Error('TK consent denied');
    }
  }
  // Mock interaction check
  return { warnings: [], total_warnings: 0 };
}

async function applyBiasCorrection(env: Env, baseDose: number, compoundId: string, demographics: Record<string, any>, tier: string): Promise<any> {
  // Query demographic factors from D1
  const stmt = env.NEUROBOTANICA_DB.prepare('SELECT cyp450_adjustment, dosing_adjustment FROM neurobotanica_demographic_factors WHERE demographic_category = ? AND demographic_value = ?');
  const result = await stmt.bind('age', demographics.age || 'unknown').all();
  const adjustment = result.results[0] || { cyp450_adjustment: 1.0, dosing_adjustment: 1.0 };
  const adjustedDose = baseDose * adjustment.cyp450_adjustment * adjustment.dosing_adjustment;
  return { adjusted_dose_mg: adjustedDose, factors_applied: adjustment };
}

async function predictSynergy(env: Env, a: string, b: string, tier: string): Promise<any> {
  // Cache check
  const cacheKey = `synergy:${a}:${b}:${tier}`;
  let cached = await env.NEUROBOTANICA_CACHE.get(cacheKey);
  if (cached) return JSON.parse(cached);
  
  // TK enhancement
  let tkBoost = 0.0;
  if (tier === 'tk_enhanced') {
    const tkData = await queryTKData(env, a, b);
    if (tkData && await verifyConsent(env, tkData.consent_id)) {
      tkBoost = 0.2;
    }
  }
  
  // Mock synergy calculation
  const score = Math.min(0.5 + tkBoost, 1.0);
  const result = { synergy_score: score, tk_enhanced: tkBoost > 0 };
  
  // Cache result
  await env.NEUROBOTANICA_CACHE.put(cacheKey, JSON.stringify(result), { expirationTtl: 3600 });
  
  return result;
}

async function analyzePlant(env: Env, plantId: string, tier: string): Promise<any> {
  // Mock plant analysis
  return { profile: 'cannabinoid_ratios', compounds: ['cbd', 'thc'] };
}

async function predictPolysaccharides(env: Env, compoundIds: string[], tier: string): Promise<any> {
  // Mock polysaccharide prediction
  return { effects: 'microbiome_modulation', confidence: 0.85 };
}

async function verifyConsent(env: Env, consentId: string): Promise<boolean> {
  const stmt = env.NEUROBOTANICA_DB.prepare('SELECT consent_status FROM omnipath_consent_artifacts WHERE consent_id = ?');
  const result = await stmt.bind(consentId).all();
  return result.results[0]?.consent_status === 'active';
}

async function queryTKData(env: Env, a: string, b: string): Promise<any> {
  const stmt = env.NEUROBOTANICA_DB.prepare('SELECT consent_id FROM neurobotanica_synergy_predictions WHERE compound_a_id = ? AND compound_b_id = ? AND requires_consent = 1');
  const result = await stmt.bind(a, b).all();
  return result.results[0];
}