/**
 * NeuroBotanica Cloudflare Worker
 * Edge-deployed API for cannabis optimization with TK preservation
 * Performance Target: <200ms global responses
 */

export default {
  async fetch(request, env) {
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
        const reqBody = await request.json();
        const startTime = Date.now();

        // Phase 1: Interactions and Bias
        const interactions = await checkInteractions(env, reqBody.compound_ids, reqBody.customer_tier || 'computational_only');
        const biasCorrected = await applyBiasCorrection(env, 10.0, reqBody.compound_ids[0], reqBody.demographics || {}, reqBody.customer_tier || 'computational_only');

        // Phase 2: Synergy, Plant, Polysaccharides
        const synergy = await predictSynergy(env, reqBody.compound_ids[0], reqBody.compound_ids[1] || reqBody.compound_ids[0], reqBody.customer_tier || 'computational_only');
        const plantProfile = reqBody.plant_id ? await analyzePlant(env, reqBody.plant_id, reqBody.customer_tier || 'computational_only') : {};
        const polysaccharideEffects = await predictPolysaccharides(env, reqBody.compound_ids, reqBody.customer_tier || 'computational_only');

        const result = {
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
        return new Response(JSON.stringify({ error: 'Analysis failed', details: error.message, stack: error.stack }), {
          status: 500,
          headers: { 'Content-Type': 'application/json' }
        });
      }
    }

    return new Response('Not Found', { status: 404 });
  }
};

// Helper functions for engines (simplified for edge deployment)
async function checkInteractions(env, compoundIds, tier) {
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

async function applyBiasCorrection(env, baseDose, compoundId, demographics, tier) {
  // Query demographic factors from D1
  let adjustmentFactor = 1.0;
  
  try {
    if (demographics.age) {
      const ageGroup = demographics.age < 30 ? 'young' : demographics.age > 65 ? 'elderly' : 'adult';
      const stmt = env.NEUROBOTANICA_DB.prepare(`
        SELECT adjustment_factor 
        FROM neurobotanica_demographic_factors 
        WHERE compound_id = ? AND demographic_group = ?
      `);
      const result = await stmt.bind(compoundId, ageGroup).all();
      if (result.results.length > 0) {
        adjustmentFactor *= result.results[0].adjustment_factor || 1.0;
      }
    }
    
    // Weight adjustment
    if (demographics.weight) {
      const weightFactor = demographics.weight < 150 ? 0.9 : demographics.weight > 250 ? 1.1 : 1.0;
      adjustmentFactor *= weightFactor;
    }
    
    // Gender adjustment (if available)
    if (demographics.gender && demographics.gender.toLowerCase() === 'female') {
      adjustmentFactor *= 0.95; // Slightly lower dose for females
    }
  } catch (error) {
    console.error('Bias correction query failed:', error);
    adjustmentFactor = 1.0; // Fallback to no adjustment
  }
  
  const adjustedDose = baseDose * adjustmentFactor;
  return { 
    adjusted_dose_mg: Math.round(adjustedDose * 100) / 100, 
    factors_applied: { adjustment_factor: adjustmentFactor, demographics_considered: Object.keys(demographics) }
  };
}

async function predictSynergy(env, a, b, tier) {
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

async function analyzePlant(env, plantId, tier) {
  // Mock plant analysis
  return { profile: 'cannabinoid_ratios', compounds: ['cbd', 'thc'] };
}

async function predictPolysaccharides(env, compoundIds, tier) {
  // Query D1 for polysaccharide effects
  const stmt = env.NEUROBOTANICA_DB.prepare(`
    SELECT effect_type, confidence_score, modulation_type
    FROM neurobotanica_polysaccharides
    WHERE compound_id = ? AND tier_access <= ?
    ORDER BY confidence_score DESC
    LIMIT 1
  `);
  
  const compoundId = compoundIds[0] || 'cbd';
  const tierLevel = tier === 'tk_enhanced' ? 2 : 1;
  
  try {
    const result = await stmt.bind(compoundId, tierLevel).all();
    if (result.results.length > 0) {
      const data = result.results[0];
      return { 
        effects: data.effect_type || 'microbiome_modulation', 
        confidence: data.confidence_score || 0.85,
        modulation: data.modulation_type || 'beneficial'
      };
    }
  } catch (error) {
    console.error('Polysaccharide query failed:', error);
  }
  
  // Fallback
  return { effects: 'microbiome_modulation', confidence: 0.75 };
}

async function verifyConsent(env, consentId) {
  const stmt = env.NEUROBOTANICA_DB.prepare('SELECT consent_status FROM omnipath_consent_artifacts WHERE consent_id = ?');
  const result = await stmt.bind(consentId).all();
  return result.results[0]?.consent_status === 'active';
}

async function queryTKData(env, a, b) {
  const stmt = env.NEUROBOTANICA_DB.prepare('SELECT consent_id FROM neurobotanica_synergy_predictions WHERE compound_a_id = ? AND compound_b_id = ? AND requires_consent = 1');
  const result = await stmt.bind(a, b).all();
  return result.results[0];
}