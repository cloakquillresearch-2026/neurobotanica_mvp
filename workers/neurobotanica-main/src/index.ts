/**
 * NeuroBotanica Main API Worker (TypeScript)
 * Edge-deployed API for cannabis optimization with TK preservation
 * Performance Target: <200ms global responses, <10s full analysis
 */

export interface Env {
  NEUROBOTANICA_DB: D1Database;
  NEUROBOTANICA_CACHE?: KVNamespace;
  ENVIRONMENT?: string;
}

interface Demographics {
  age?: number;
  gender?: string;
  weight?: number;
  condition?: string;
  ethnicity?: string;
  metabolism_type?: 'fast' | 'normal' | 'slow';
  [key: string]: unknown;
}

interface DemographicFactor {
  adjustment_factor: number;
  evidence_summary: string;
  confidence_level?: number;
  source?: string;
}

interface BiasResult {
  adjusted_dose_mg: number;
  factors_applied: {
    base_dose: number;
    adjustment_factor: number;
    demographics_considered: string[];
    adjustments_breakdown: Record<string, number>;
  };
  evidence: string;
  confidence: number;
  warnings: string[];
}

interface SynergyResult {
  synergy_score: number;
  tk_enhanced: boolean;
  evidence: string;
  confidence_level: number;
  recommended_timing?: string;
}

interface PolysaccharideResult {
  effects: string;
  confidence: number;
  modulation: string;
}

interface AnalysisResult {
  interactions: { warnings: string[]; total_warnings: number };
  bias_correction: BiasResult;
  synergy: SynergyResult;
  plant_profile: Record<string, unknown>;
  polysaccharide_effects: PolysaccharideResult;
  processing_time_ms: number;
}

// CORS headers helper
function getCorsHeaders(request: Request): Record<string, string> {
  const origin = request.headers.get('Origin') || '*';
  return {
    'Access-Control-Allow-Origin': origin,
    'Access-Control-Allow-Methods': 'GET, POST, OPTIONS',
    'Access-Control-Allow-Headers': 'Content-Type, Authorization',
    'Access-Control-Allow-Credentials': 'true',
  };
}

// JSON response helper with error handling
function jsonResponse(data: unknown, status = 200, request?: Request): Response {
  const headers: Record<string, string> = {
    'Content-Type': 'application/json',
    ...(request ? getCorsHeaders(request) : {}),
  };
  return new Response(JSON.stringify(data), { status, headers });
}

// Error response helper
function errorResponse(
  message: string,
  status: number,
  request: Request,
  details?: string
): Response {
  return jsonResponse(
    { error: message, details, timestamp: new Date().toISOString() },
    status,
    request
  );
}

export default {
  async fetch(request: Request, env: Env): Promise<Response> {
    const url = new URL(request.url);
    const corsHeaders = getCorsHeaders(request);

    // Handle CORS preflight
    if (request.method === 'OPTIONS') {
      return new Response(null, { status: 204, headers: corsHeaders });
    }

    // Health check endpoint
    if (url.pathname === '/api/neurobotanica/health') {
      const dbStatus = await checkDatabaseHealth(env);
      return jsonResponse(
        {
          status: dbStatus.healthy ? 'healthy' : 'degraded',
          engines: ['interactions', 'bias', 'synergy', 'plant', 'polysaccharides'],
          database: dbStatus,
          version: '2.0.0',
          timestamp: new Date().toISOString(),
        },
        dbStatus.healthy ? 200 : 503,
        request
      );
    }

    // Analyze endpoint
    if (url.pathname === '/api/neurobotanica/analyze' && request.method === 'POST') {
      return handleAnalyzeRequest(request, env);
    }

    // Bias correction standalone endpoint
    if (url.pathname === '/api/neurobotanica/bias-correction' && request.method === 'POST') {
      return handleBiasCorrectionRequest(request, env);
    }

    return errorResponse('Not Found', 404, request);
  },
};

async function checkDatabaseHealth(env: Env): Promise<{ healthy: boolean; latency_ms: number; error?: string }> {
  const start = Date.now();
  try {
    if (!env.NEUROBOTANICA_DB) {
      return { healthy: false, latency_ms: 0, error: 'Database not configured' };
    }
    await env.NEUROBOTANICA_DB.prepare('SELECT 1').first();
    return { healthy: true, latency_ms: Date.now() - start };
  } catch (error) {
    return {
      healthy: false,
      latency_ms: Date.now() - start,
      error: error instanceof Error ? error.message : 'Unknown database error',
    };
  }
}

async function handleAnalyzeRequest(request: Request, env: Env): Promise<Response> {
  try {
    const reqBody = await request.json() as {
      compound_ids: string[];
      demographics?: Demographics;
      customer_tier?: string;
      plant_id?: string;
    };

    const startTime = Date.now();

    // Validate required fields
    if (!reqBody.compound_ids || !Array.isArray(reqBody.compound_ids) || reqBody.compound_ids.length === 0) {
      return errorResponse('compound_ids is required and must be a non-empty array', 400, request);
    }

    const tier = reqBody.customer_tier || 'computational_only';
    const demographics = reqBody.demographics || {};
    const compoundIds = reqBody.compound_ids;

    // Run Phase 1 and Phase 2 in parallel where possible
    const [interactions, biasCorrected, synergy, polysaccharideEffects] = await Promise.all([
      checkInteractions(env, compoundIds, tier),
      applyBiasCorrection(env, 10.0, compoundIds[0], demographics, tier),
      predictSynergy(env, compoundIds[0], compoundIds[1] || compoundIds[0], tier),
      predictPolysaccharides(env, compoundIds, tier),
    ]);

    // Optional plant profile (sequential as it's optional)
    const plantProfile = reqBody.plant_id
      ? await analyzePlant(env, reqBody.plant_id, tier)
      : {};

    const result: AnalysisResult = {
      interactions,
      bias_correction: biasCorrected,
      synergy,
      plant_profile: plantProfile,
      polysaccharide_effects: polysaccharideEffects,
      processing_time_ms: Date.now() - startTime,
    };

    return jsonResponse(result, 200, request);
  } catch (error) {
    console.error('Analysis error:', error);
    return errorResponse(
      'Analysis failed',
      500,
      request,
      error instanceof Error ? error.message : 'Unknown error'
    );
  }
}

async function handleBiasCorrectionRequest(request: Request, env: Env): Promise<Response> {
  try {
    const reqBody = await request.json() as {
      base_dose: number;
      compound_id: string;
      demographics: Demographics;
      tier?: string;
    };

    if (!reqBody.compound_id || reqBody.base_dose === undefined) {
      return errorResponse('compound_id and base_dose are required', 400, request);
    }

    const result = await applyBiasCorrection(
      env,
      reqBody.base_dose,
      reqBody.compound_id,
      reqBody.demographics || {},
      reqBody.tier || 'computational_only'
    );

    return jsonResponse(result, 200, request);
  } catch (error) {
    console.error('Bias correction error:', error);
    return errorResponse(
      'Bias correction failed',
      500,
      request,
      error instanceof Error ? error.message : 'Unknown error'
    );
  }
}

async function checkInteractions(
  env: Env,
  compoundIds: string[],
  tier: string
): Promise<{ warnings: string[]; total_warnings: number }> {
  // TK consent check for enhanced tier
  if (tier === 'tk_enhanced') {
    const consentValid = await verifyConsent(env, 'sample_consent_id');
    if (!consentValid) {
      throw new Error('TK consent denied');
    }
  }

  try {
    // Query drug interactions from database
    const warnings: string[] = [];

    for (const compoundId of compoundIds) {
      const stmt = env.NEUROBOTANICA_DB.prepare(`
        SELECT interaction_type, severity, description
        FROM neurobotanica_drug_interactions
        WHERE compound_id = ? OR interacting_compound_id = ?
        ORDER BY severity DESC
        LIMIT 5
      `);

      const result = await stmt.bind(compoundId, compoundId).all();

      if (result.results && result.results.length > 0) {
        for (const row of result.results) {
          const r = row as { interaction_type: string; severity: string; description: string };
          if (r.severity === 'high' || r.severity === 'moderate') {
            warnings.push(`${r.interaction_type}: ${r.description}`);
          }
        }
      }
    }

    return { warnings, total_warnings: warnings.length };
  } catch (error) {
    // Log but don't fail - return empty warnings on DB error
    console.error('Interaction check failed:', error);
    return { warnings: [], total_warnings: 0 };
  }
}

/**
 * Enhanced applyBiasCorrection function
 * Queries D1 database for demographic adjustment factors and applies real bias corrections
 * Supports: age, gender, weight, condition, ethnicity, metabolism type
 */
async function applyBiasCorrection(
  env: Env,
  baseDose: number,
  compoundId: string,
  demographics: Demographics,
  tier: string
): Promise<BiasResult> {
  let adjustmentFactor = 1.0;
  const evidenceParts: string[] = [];
  const adjustmentsBreakdown: Record<string, number> = {};
  const warnings: string[] = [];
  let overallConfidence = 0.75; // Default confidence

  // Validate database connection
  if (!env.NEUROBOTANICA_DB) {
    return {
      adjusted_dose_mg: Math.round(baseDose * 100) / 100,
      factors_applied: {
        base_dose: baseDose,
        adjustment_factor: 1.0,
        demographics_considered: [],
        adjustments_breakdown: {},
      },
      evidence: 'Database not configured. Using standard dose.',
      confidence: 0.5,
      warnings: ['D1 database not available'],
    };
  }

  try {
    // 1. Condition-specific adjustments (highest priority)
    if (demographics.condition) {
      const conditionResult = await queryConditionFactor(env, compoundId, demographics.condition);
      if (conditionResult) {
        adjustmentFactor *= conditionResult.adjustment_factor;
        adjustmentsBreakdown['condition'] = conditionResult.adjustment_factor;
        evidenceParts.push(conditionResult.evidence_summary);
        if (conditionResult.confidence_level) {
          overallConfidence = Math.min(overallConfidence, conditionResult.confidence_level);
        }
      }
    }

    // 2. Age-based adjustments
    if (demographics.age !== undefined && demographics.age > 0) {
      const ageResult = await queryAgeFactor(env, compoundId, demographics.age);
      if (ageResult) {
        adjustmentFactor *= ageResult.adjustment_factor;
        adjustmentsBreakdown['age'] = ageResult.adjustment_factor;
        evidenceParts.push(ageResult.evidence_summary);
      } else {
        // Fallback age-based calculation
        const ageFactor = calculateAgeFactor(demographics.age);
        adjustmentFactor *= ageFactor;
        adjustmentsBreakdown['age'] = ageFactor;
        evidenceParts.push(`Age-based adjustment (${demographics.age}y): ${(ageFactor * 100).toFixed(0)}%`);
      }
    }

    // 3. Weight-based adjustments
    if (demographics.weight !== undefined && demographics.weight > 0) {
      const weightResult = await queryWeightFactor(env, compoundId, demographics.weight);
      if (weightResult) {
        adjustmentFactor *= weightResult.adjustment_factor;
        adjustmentsBreakdown['weight'] = weightResult.adjustment_factor;
        evidenceParts.push(weightResult.evidence_summary);
      } else {
        // Fallback weight-based calculation (in lbs)
        const weightFactor = calculateWeightFactor(demographics.weight);
        adjustmentFactor *= weightFactor;
        adjustmentsBreakdown['weight'] = weightFactor;
        evidenceParts.push(`Weight-based adjustment (${demographics.weight}lbs): ${(weightFactor * 100).toFixed(0)}%`);
      }
    }

    // 4. Gender-based adjustments
    if (demographics.gender) {
      const genderResult = await queryGenderFactor(env, compoundId, demographics.gender);
      if (genderResult) {
        adjustmentFactor *= genderResult.adjustment_factor;
        adjustmentsBreakdown['gender'] = genderResult.adjustment_factor;
        evidenceParts.push(genderResult.evidence_summary);
      } else {
        // Fallback gender-based calculation
        const genderFactor = calculateGenderFactor(demographics.gender);
        if (genderFactor !== 1.0) {
          adjustmentFactor *= genderFactor;
          adjustmentsBreakdown['gender'] = genderFactor;
          evidenceParts.push(`Gender-based adjustment: ${(genderFactor * 100).toFixed(0)}%`);
        }
      }
    }

    // 5. Ethnicity-based adjustments (pharmacogenomic considerations)
    if (demographics.ethnicity) {
      const ethnicityResult = await queryEthnicityFactor(env, compoundId, demographics.ethnicity);
      if (ethnicityResult) {
        adjustmentFactor *= ethnicityResult.adjustment_factor;
        adjustmentsBreakdown['ethnicity'] = ethnicityResult.adjustment_factor;
        evidenceParts.push(ethnicityResult.evidence_summary);
      }
    }

    // 6. Metabolism type adjustments
    if (demographics.metabolism_type) {
      const metabolismFactor = calculateMetabolismFactor(demographics.metabolism_type);
      adjustmentFactor *= metabolismFactor;
      adjustmentsBreakdown['metabolism'] = metabolismFactor;
      evidenceParts.push(`Metabolism type (${demographics.metabolism_type}): ${(metabolismFactor * 100).toFixed(0)}%`);
    }

    // Ensure adjustment factor stays within safe bounds (0.5x - 2.0x)
    if (adjustmentFactor < 0.5) {
      warnings.push('Adjustment factor capped at minimum 0.5x for safety');
      adjustmentFactor = 0.5;
    } else if (adjustmentFactor > 2.0) {
      warnings.push('Adjustment factor capped at maximum 2.0x for safety');
      adjustmentFactor = 2.0;
    }

    // TK-enhanced tier can access additional cultural medicine data
    if (tier === 'tk_enhanced') {
      const tkResult = await queryTKDemographicFactors(env, compoundId, demographics);
      if (tkResult) {
        adjustmentFactor *= tkResult.adjustment_factor;
        adjustmentsBreakdown['tk_cultural'] = tkResult.adjustment_factor;
        evidenceParts.push(`TK-enhanced: ${tkResult.evidence_summary}`);
        overallConfidence = Math.max(overallConfidence, 0.85);
      }
    }

  } catch (error) {
    console.error('Bias correction query failed:', error);
    warnings.push('Database query failed, using fallback calculations');

    // Apply fallback calculations on error
    if (demographics.age) {
      adjustmentFactor *= calculateAgeFactor(demographics.age);
    }
    if (demographics.weight) {
      adjustmentFactor *= calculateWeightFactor(demographics.weight);
    }
    if (demographics.gender) {
      adjustmentFactor *= calculateGenderFactor(demographics.gender);
    }

    evidenceParts.push('Fallback calculations applied due to database error');
    overallConfidence = 0.5;
  }

  const adjustedDose = baseDose * adjustmentFactor;
  const evidence = evidenceParts.length > 0
    ? evidenceParts.join('. ')
    : 'Standard adjustment applied - no specific demographic data available.';

  return {
    adjusted_dose_mg: Math.round(adjustedDose * 100) / 100,
    factors_applied: {
      base_dose: baseDose,
      adjustment_factor: Math.round(adjustmentFactor * 1000) / 1000,
      demographics_considered: Object.keys(adjustmentsBreakdown),
      adjustments_breakdown: adjustmentsBreakdown,
    },
    evidence,
    confidence: overallConfidence,
    warnings,
  };
}

// Database query helpers for demographic factors
async function queryConditionFactor(
  env: Env,
  compoundId: string,
  condition: string
): Promise<DemographicFactor | null> {
  // Map common condition names to database values
  const conditionMap: Record<string, string> = {
    weight_management: 'metabolism',
    muscle_spasms: 'spasms',
    chronic_pain: 'pain',
    anxiety: 'anxiety',
    insomnia: 'sleep',
    inflammation: 'inflammation',
  };

  const mappedCondition = conditionMap[condition.toLowerCase()] || condition.toLowerCase();

  const stmt = env.NEUROBOTANICA_DB.prepare(`
    SELECT adjustment_factor, evidence_summary, confidence_level, source
    FROM neurobotanica_demographic_factors
    WHERE compound_id = ? AND condition = ?
    LIMIT 1
  `);

  const result = await stmt.bind(compoundId, mappedCondition).first<DemographicFactor>();
  return result || null;
}

async function queryAgeFactor(
  env: Env,
  compoundId: string,
  age: number
): Promise<DemographicFactor | null> {
  const ageGroup = age < 25 ? 'young' : age > 65 ? 'elderly' : age > 50 ? 'mature' : 'adult';

  const stmt = env.NEUROBOTANICA_DB.prepare(`
    SELECT adjustment_factor, evidence_summary, confidence_level
    FROM neurobotanica_demographic_factors
    WHERE compound_id = ? AND demographic_group = ?
    LIMIT 1
  `);

  const result = await stmt.bind(compoundId, ageGroup).first<DemographicFactor>();
  return result || null;
}

async function queryWeightFactor(
  env: Env,
  compoundId: string,
  weight: number
): Promise<DemographicFactor | null> {
  // Weight ranges in lbs
  const weightGroup = weight < 120 ? 'light' : weight > 220 ? 'heavy' : 'average';

  const stmt = env.NEUROBOTANICA_DB.prepare(`
    SELECT adjustment_factor, evidence_summary, confidence_level
    FROM neurobotanica_demographic_factors
    WHERE compound_id = ? AND weight_category = ?
    LIMIT 1
  `);

  const result = await stmt.bind(compoundId, weightGroup).first<DemographicFactor>();
  return result || null;
}

async function queryGenderFactor(
  env: Env,
  compoundId: string,
  gender: string
): Promise<DemographicFactor | null> {
  const normalizedGender = gender.toLowerCase();

  const stmt = env.NEUROBOTANICA_DB.prepare(`
    SELECT adjustment_factor, evidence_summary, confidence_level
    FROM neurobotanica_demographic_factors
    WHERE compound_id = ? AND gender = ?
    LIMIT 1
  `);

  const result = await stmt.bind(compoundId, normalizedGender).first<DemographicFactor>();
  return result || null;
}

async function queryEthnicityFactor(
  env: Env,
  compoundId: string,
  ethnicity: string
): Promise<DemographicFactor | null> {
  const stmt = env.NEUROBOTANICA_DB.prepare(`
    SELECT adjustment_factor, evidence_summary, confidence_level
    FROM neurobotanica_demographic_factors
    WHERE compound_id = ? AND ethnicity = ?
    LIMIT 1
  `);

  const result = await stmt.bind(compoundId, ethnicity.toLowerCase()).first<DemographicFactor>();
  return result || null;
}

async function queryTKDemographicFactors(
  env: Env,
  compoundId: string,
  demographics: Demographics
): Promise<DemographicFactor | null> {
  // Query traditional knowledge database for cultural medicine considerations
  const stmt = env.NEUROBOTANICA_DB.prepare(`
    SELECT adjustment_factor, evidence_summary, confidence_level
    FROM neurobotanica_tk_demographic_factors
    WHERE compound_id = ?
    LIMIT 1
  `);

  const result = await stmt.bind(compoundId).first<DemographicFactor>();
  return result || null;
}

// Fallback calculation functions
function calculateAgeFactor(age: number): number {
  if (age < 25) return 0.85; // Younger metabolism, lower dose
  if (age > 75) return 0.8;  // Elderly, reduced clearance
  if (age > 65) return 0.85; // Senior, slightly reduced
  if (age > 50) return 0.95; // Mature adult
  return 1.0; // Adult (25-50)
}

function calculateWeightFactor(weightLbs: number): number {
  if (weightLbs < 120) return 0.85;
  if (weightLbs < 150) return 0.9;
  if (weightLbs > 280) return 1.15;
  if (weightLbs > 250) return 1.1;
  if (weightLbs > 200) return 1.05;
  return 1.0; // 150-200 lbs
}

function calculateGenderFactor(gender: string): number {
  const normalized = gender.toLowerCase();
  if (normalized === 'female' || normalized === 'f') return 0.92;
  if (normalized === 'male' || normalized === 'm') return 1.0;
  return 1.0; // Non-binary or unknown
}

function calculateMetabolismFactor(type: 'fast' | 'normal' | 'slow'): number {
  switch (type) {
    case 'fast': return 1.15;   // Fast metabolizers need higher dose
    case 'slow': return 0.85;   // Slow metabolizers need lower dose
    default: return 1.0;
  }
}

async function predictSynergy(
  env: Env,
  a: string,
  b: string,
  tier: string
): Promise<SynergyResult> {
  // Cache check
  const cacheKey = `synergy:${a}:${b}:${tier}`;
  let cached: string | null = null;

  try {
    if (env.NEUROBOTANICA_CACHE) {
      cached = await env.NEUROBOTANICA_CACHE.get(cacheKey);
    }
  } catch (cacheError) {
    console.warn('Cache not available:', cacheError instanceof Error ? cacheError.message : 'Unknown');
  }

  if (cached) {
    return JSON.parse(cached);
  }

  let synergyScore = 0.5;
  let tkEnhanced = false;
  let evidence = 'Computational prediction';
  let confidenceLevel = 0.6;
  let recommendedTiming: string | undefined;

  try {
    const stmt = env.NEUROBOTANICA_DB.prepare(`
      SELECT synergy_score, confidence_level, evidence_summary, tk_enhanced, recommended_timing
      FROM neurobotanica_synergy_predictions
      WHERE (compound_a_id = ? AND compound_b_id = ?) OR (compound_a_id = ? AND compound_b_id = ?)
      ORDER BY confidence_level DESC
      LIMIT 1
    `);

    const result = await stmt.bind(a, b, b, a).first<{
      synergy_score: number;
      confidence_level: number;
      evidence_summary: string;
      tk_enhanced: boolean;
      recommended_timing?: string;
    }>();

    if (result) {
      synergyScore = result.synergy_score || 0.5;
      tkEnhanced = result.tk_enhanced || false;
      evidence = result.evidence_summary || 'Database prediction';
      confidenceLevel = result.confidence_level || 0.6;
      recommendedTiming = result.recommended_timing;
    } else {
      // Fallback synergy calculation based on known compound relationships
      synergyScore = calculateFallbackSynergy(a, b);
      evidence = 'Calculated synergy - no database data available';
    }
  } catch (error) {
    console.error('Synergy query failed:', error);
    synergyScore = calculateFallbackSynergy(a, b);
    evidence = 'Query failed, using calculated synergy';
  }

  // TK enhancement
  if (tier === 'tk_enhanced' && !tkEnhanced) {
    try {
      const tkData = await queryTKData(env, a, b);
      if (tkData && await verifyConsent(env, tkData.consent_id)) {
        synergyScore = Math.min(synergyScore + 0.15, 1.0);
        tkEnhanced = true;
        evidence += ' (TK-enhanced)';
        confidenceLevel = Math.min(confidenceLevel + 0.1, 0.95);
      }
    } catch (error) {
      console.error('TK enhancement failed:', error);
    }
  }

  const result: SynergyResult = {
    synergy_score: Math.round(synergyScore * 1000) / 1000,
    tk_enhanced: tkEnhanced,
    evidence,
    confidence_level: confidenceLevel,
    recommended_timing: recommendedTiming,
  };

  // Cache result
  try {
    if (env.NEUROBOTANICA_CACHE) {
      await env.NEUROBOTANICA_CACHE.put(cacheKey, JSON.stringify(result), { expirationTtl: 3600 });
    }
  } catch (cacheError) {
    console.warn('Cache storage failed:', cacheError instanceof Error ? cacheError.message : 'Unknown');
  }

  return result;
}

function calculateFallbackSynergy(a: string, b: string): number {
  // Known synergistic pairs
  const synergyPairs: Record<string, number> = {
    'cbd:thc': 0.75,
    'thc:cbd': 0.75,
    'cbd:cbg': 0.7,
    'cbg:cbd': 0.7,
    'thc:myrcene': 0.8,
    'myrcene:thc': 0.8,
    'cbd:linalool': 0.72,
    'linalool:cbd': 0.72,
  };

  const key = `${a.toLowerCase()}:${b.toLowerCase()}`;
  return synergyPairs[key] || 0.5 + (Math.random() * 0.2 - 0.1); // 0.4-0.6 for unknown pairs
}

async function analyzePlant(
  env: Env,
  plantId: string,
  tier: string
): Promise<Record<string, unknown>> {
  try {
    const stmt = env.NEUROBOTANICA_DB.prepare(`
      SELECT strain_name, cannabinoid_profile, terpene_profile, effects, lineage
      FROM neurobotanica_plant_profiles
      WHERE plant_id = ?
      LIMIT 1
    `);

    const result = await stmt.bind(plantId).first<{
      strain_name: string;
      cannabinoid_profile: string;
      terpene_profile: string;
      effects: string;
      lineage: string;
    }>();

    if (result) {
      return {
        profile: 'full_spectrum',
        strain: result.strain_name,
        cannabinoids: JSON.parse(result.cannabinoid_profile || '{}'),
        terpenes: JSON.parse(result.terpene_profile || '[]'),
        effects: JSON.parse(result.effects || '[]'),
        lineage: result.lineage,
      };
    }
  } catch (error) {
    console.error('Plant analysis failed:', error);
  }

  // Fallback
  return { profile: 'cannabinoid_ratios', compounds: ['cbd', 'thc'] };
}

async function predictPolysaccharides(
  env: Env,
  compoundIds: string[],
  tier: string
): Promise<PolysaccharideResult> {
  const compoundId = compoundIds[0] || 'cbd';
  const tierLevel = tier === 'tk_enhanced' ? 2 : 1;

  try {
    const stmt = env.NEUROBOTANICA_DB.prepare(`
      SELECT effect_type, confidence_score, modulation_type, mechanism
      FROM neurobotanica_polysaccharides
      WHERE compound_id = ? AND tier_access <= ?
      ORDER BY confidence_score DESC
      LIMIT 1
    `);

    const result = await stmt.bind(compoundId, tierLevel).first<{
      effect_type: string;
      confidence_score: number;
      modulation_type: string;
      mechanism?: string;
    }>();

    if (result) {
      return {
        effects: result.effect_type || 'microbiome_modulation',
        confidence: result.confidence_score || 0.85,
        modulation: result.modulation_type || 'beneficial',
      };
    }
  } catch (error) {
    console.error('Polysaccharide query failed:', error);
  }

  // Fallback
  return { effects: 'microbiome_modulation', confidence: 0.75, modulation: 'beneficial' };
}

async function verifyConsent(env: Env, consentId: string): Promise<boolean> {
  try {
    const stmt = env.NEUROBOTANICA_DB.prepare(
      'SELECT consent_status, expiry_date FROM omnipath_consent_artifacts WHERE consent_id = ?'
    );
    const result = await stmt.bind(consentId).first<{ consent_status: string; expiry_date?: string }>();

    if (!result || result.consent_status !== 'active') {
      return false;
    }

    // Check expiry if present
    if (result.expiry_date) {
      const expiry = new Date(result.expiry_date);
      if (expiry < new Date()) {
        return false;
      }
    }

    return true;
  } catch (error) {
    console.error('Consent verification failed:', error);
    return false;
  }
}

async function queryTKData(
  env: Env,
  a: string,
  b: string
): Promise<{ consent_id: string } | null> {
  try {
    const stmt = env.NEUROBOTANICA_DB.prepare(
      'SELECT consent_id FROM neurobotanica_synergy_predictions WHERE compound_a_id = ? AND compound_b_id = ? AND requires_consent = 1'
    );
    const result = await stmt.bind(a, b).first<{ consent_id: string }>();
    return result || null;
  } catch {
    return null;
  }
}
