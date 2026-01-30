// =============================================================================
// DISPENSARY ENDPOINTS - Stub implementations for Budtender App
// =============================================================================

// Helper: Generate a mock profile ID
function generateProfileId() {
  return `prof_${Math.random().toString(36).substring(2, 8).toUpperCase()}`;
}

// Helper: Generate a profile code
function generateProfileCode() {
  return `NB-${Math.random().toString(36).substring(2, 7).toUpperCase()}`;
}

// In-memory storage (will reset on worker restart)
const profiles = new Map();
const recommendations = new Map();
const feedback = new Map();

// --- PROFILE ENDPOINTS ---

// POST /api/dispensary/profile - Create customer profile
async function handleCreateProfile(request) {
  const body = await request.json();
  const profileId = generateProfileId();
  const profileCode = generateProfileCode();

  const profile = {
    profile_id: profileId,
    profile_code: profileCode,
    created_at: new Date().toISOString(),
    completeness_score: 0.75,
    primary_condition: body.conditions?.[0]?.name || null,
    ...body
  };

  profiles.set(profileId, profile);

  return new Response(JSON.stringify({
    profile_id: profileId,
    profile_code: profileCode,
    created_at: profile.created_at,
    completeness_score: profile.completeness_score,
    primary_condition: profile.primary_condition
  }), {
    status: 200,
    headers: { 'Content-Type': 'application/json' }
  });
}

// GET /api/dispensary/profile/:id - Get customer profile
async function handleGetProfile(profileId) {
  const profile = profiles.get(profileId);

  if (!profile) {
    return new Response(JSON.stringify({ detail: 'Profile not found' }), {
      status: 404,
      headers: { 'Content-Type': 'application/json' }
    });
  }

  return new Response(JSON.stringify(profile), {
    status: 200,
    headers: { 'Content-Type': 'application/json' }
  });
}

// PUT /api/dispensary/profile/:id - Update customer profile
async function handleUpdateProfile(request, profileId) {
  const body = await request.json();
  let profile = profiles.get(profileId);

  if (!profile) {
    // Create if doesn't exist
    profile = {
      profile_id: profileId,
      profile_code: generateProfileCode(),
      created_at: new Date().toISOString()
    };
  }

  const updatedProfile = {
    ...profile,
    ...body,
    updated_at: new Date().toISOString(),
    completeness_score: 0.85,
    primary_condition: body.conditions?.[0]?.name || profile.primary_condition
  };

  profiles.set(profileId, updatedProfile);

  return new Response(JSON.stringify({
    profile_id: profileId,
    created_at: updatedProfile.created_at,
    completeness_score: updatedProfile.completeness_score,
    primary_condition: updatedProfile.primary_condition
  }), {
    status: 200,
    headers: { 'Content-Type': 'application/json' }
  });
}

// --- TRANSACTION ENDPOINT ---

// POST /api/dispensary/transaction - Create transaction
async function handleCreateTransaction(request) {
  const body = await request.json();
  const transactionId = `txn_${Math.random().toString(36).substring(2, 10).toUpperCase()}`;

  return new Response(JSON.stringify({
    transaction_id: transactionId,
    status: 'completed',
    created_at: new Date().toISOString(),
    profile_id: body.profile_id || null,
    products: body.products || []
  }), {
    status: 200,
    headers: { 'Content-Type': 'application/json' }
  });
}

// --- RECOMMENDATION ENDPOINT ---

// POST /api/dispensary/recommend - Product recommendations
async function handleRecommend(request) {
  const body = await request.json();
  const recId = `rec_${Math.random().toString(36).substring(2, 8).toUpperCase()}`;

  // Mock recommendation based on primary condition
  const primaryCondition = body.customer_profile?.conditions?.[0]?.name || 'general_wellness';

  const mockRecommendations = [
    {
      product_id: 'prod_001',
      product_name: 'Blue Dream - Hybrid',
      product_type: 'flower',
      thc_percent: 18.5,
      cbd_percent: 2.0,
      confidence_score: 0.87,
      rationale: `Recommended for ${primaryCondition} based on 505+ clinical studies`,
      terpene_profile: { myrcene: 0.4, limonene: 0.3, pinene: 0.2 }
    },
    {
      product_id: 'prod_002',
      product_name: 'CBD Relief Tincture',
      product_type: 'tincture',
      thc_percent: 1.0,
      cbd_percent: 25.0,
      confidence_score: 0.82,
      rationale: 'High CBD for therapeutic benefits without psychoactivity',
      terpene_profile: { linalool: 0.5, beta_caryophyllene: 0.3 }
    }
  ];

  recommendations.set(recId, { id: recId, recommendations: mockRecommendations });

  return new Response(JSON.stringify({
    recommendation_id: recId,
    recommendations: mockRecommendations,
    clinical_studies_referenced: 505,
    generated_at: new Date().toISOString()
  }), {
    status: 200,
    headers: { 'Content-Type': 'application/json' }
  });
}

// --- FEEDBACK ENDPOINT ---

// POST /api/dispensary/feedback - Submit feedback
async function handleFeedback(request) {
  const body = await request.json();

  if (!body.recommendation_id) {
    return new Response(JSON.stringify({ detail: 'recommendation_id required' }), {
      status: 400,
      headers: { 'Content-Type': 'application/json' }
    });
  }

  const feedbackList = feedback.get(body.recommendation_id) || [];
  feedbackList.push({ ...body, submitted_at: new Date().toISOString() });
  feedback.set(body.recommendation_id, feedbackList);

  return new Response(JSON.stringify({
    status: 'Feedback recorded',
    recommendation_id: body.recommendation_id
  }), {
    status: 200,
    headers: { 'Content-Type': 'application/json' }
  });
}

// --- STATISTICS ENDPOINT ---

// GET /api/dispensary/statistics - Get analytics
async function handleStatistics() {
  return new Response(JSON.stringify({
    total_profiles: profiles.size,
    total_recommendations: recommendations.size,
    total_feedback_entries: [...feedback.values()].reduce((sum, f) => sum + f.length, 0),
    clinical_studies_used: 505,
    conditions_covered: 22,
    trade_secret_engines: 6,
    pricing: {
      single_location: '$1,200/month',
      multi_location: '$900/month per location',
      enterprise: '$700/month per location'
    }
  }), {
    status: 200,
    headers: { 'Content-Type': 'application/json' }
  });
}

// --- INFLAMMATORY SYNERGY ENDPOINT (TS-PS-001) ---

// POST /api/dispensary/inflammatory-synergy
async function handleInflammatorySynergy(request) {
  const body = await request.json();
  const biomarkers = body.biomarkers || {};
  const availableKingdoms = body.available_kingdoms || ['cannabis', 'fungal', 'marine', 'plant'];

  // Simple heuristic matching the Python stub
  const tnf = biomarkers.tnf_alpha || 0;
  const crp = biomarkers.crp || 0;
  const il6 = biomarkers.il6 || 0;

  const primary = (tnf > 5 || crp > 3 || il6 > 2) && availableKingdoms.includes('cannabis')
    ? 'cannabis'
    : availableKingdoms[0] || 'cannabis';

  const secondary = availableKingdoms.filter(k => k !== primary).slice(0, 2);
  const rawScore = Math.min(1.0, (tnf * 0.02) + (crp * 0.03) + (il6 * 0.01));
  const confidence = 0.6 + Math.min(0.4, rawScore * 0.4);

  let recommended = ['Standard Botanical Extract'];
  if (primary === 'cannabis') recommended = ['CBD', 'CBG'];
  else if (primary === 'fungal') recommended = ['Beta-Glucan Extract'];

  const warning = (tnf === 0 && crp === 0 && il6 === 0)
    ? 'Biomarkers appear empty — results are a heuristic fallback.'
    : null;

  return new Response(JSON.stringify({
    primary_kingdom: primary,
    secondary_kingdoms: secondary,
    synergy_score: rawScore,
    confidence_level: confidence,
    recommended_compounds: recommended,
    dosing_guidance: { primary: 'start_low_titrate', frequency: 'twice_daily' },
    expected_reduction: { tnf_alpha: rawScore * 10, crp: rawScore * 5 },
    warning
  }), {
    status: 200,
    headers: { 'Content-Type': 'application/json' }
  });
}

// --- INFLAMMATORY PROFILE ENDPOINT ---

// POST /api/dispensary/profile/inflammatory
async function handleCreateInflammatoryProfile(request) {
  const body = await request.json();
  const profileId = generateProfileId();
  const profileCode = `NB-IF-${Math.random().toString(36).substring(2, 8).toUpperCase()}`;

  const profile = {
    profile_id: profileId,
    profile_code: profileCode,
    profile_type: 'inflammatory',
    created_at: new Date().toISOString(),
    completeness_score: 0.80,
    primary_condition: body.conditions?.[0]?.name || null,
    biomarkers: body.biomarkers || {},
    ...body
  };

  profiles.set(profileId, profile);

  return new Response(JSON.stringify({
    profile_id: profileId,
    created_at: profile.created_at,
    completeness_score: profile.completeness_score,
    primary_condition: profile.primary_condition
  }), {
    status: 200,
    headers: { 'Content-Type': 'application/json' }
  });
}

// --- ADJUVANT OPTIMIZATION ENDPOINT ---

// POST /api/dispensary/adjuvants/optimize
async function handleAdjuvantOptimize(request) {
  const body = await request.json();
  const primaryCompound = body.primary_compound || 'CBD';
  const therapeuticTarget = body.therapeutic_target || 'general_wellness';

  // Mock adjuvant recommendations
  const adjuvants = [
    {
      name: 'Magnesium Glycinate',
      timing: '-30 minutes',
      rationale: 'NMDA antagonist priming enhances cannabinoid receptor sensitivity',
      confidence: 0.85
    },
    {
      name: 'Omega-3 Fatty Acids',
      timing: 'with dose',
      rationale: 'Enhances bioavailability and anti-inflammatory synergy',
      confidence: 0.78
    }
  ];

  return new Response(JSON.stringify({
    primary_compound: primaryCompound,
    therapeutic_target: therapeuticTarget,
    recommended_adjuvants: adjuvants,
    optimization_score: 0.82,
    generated_at: new Date().toISOString()
  }), {
    status: 200,
    headers: { 'Content-Type': 'application/json' }
  });
}

export default {
  async fetch(request) {
    const url = new URL(request.url);
    const path = url.pathname;
    const method = request.method;

    const json = (obj) => new Response(JSON.stringify(obj), { headers: { 'Content-Type': 'application/json' } });

    // Root
    if (path === '/') {
      return json({
        message: '🌿 NeuroBotanica API - Worker shim',
        version: '0.4.0',
        status: 'operational',
        endpoints: { health: '/health', api_docs: '/docs' },
        trade_secret_engines: ['ChemPath', 'ToxPath', 'RegPath', 'BioPath', 'ClinPath', 'GenomePath']
      });
    }

    // Health
    if (path === '/health') {
      return json({ status: 'healthy', results: { ml_models: [], features: { trade_secret_engines: ['ChemPath','ToxPath','RegPath','BioPath','ClinPath','GenomePath'] } } });
    }

    // Stats
    if (path === '/api/v1/stats') {
      return json({ db: { studies: 368, compounds: 63 } });
    }

    // Dispensary routes
    if (path === '/api/dispensary/profile' && method === 'POST') {
      return handleCreateProfile(request);
    }
    if (path.startsWith('/api/dispensary/profile/') && path !== '/api/dispensary/profile/inflammatory') {
      const profileId = path.split('/').pop();
      if (method === 'GET') return handleGetProfile(profileId);
      if (method === 'PUT') return handleUpdateProfile(request, profileId);
    }
    if (path === '/api/dispensary/profile/inflammatory' && method === 'POST') {
      return handleCreateInflammatoryProfile(request);
    }
    if (path === '/api/dispensary/transaction' && method === 'POST') {
      return handleCreateTransaction(request);
    }
    if (path === '/api/dispensary/recommend' && method === 'POST') {
      return handleRecommend(request);
    }
    if (path === '/api/dispensary/feedback' && method === 'POST') {
      return handleFeedback(request);
    }
    if (path === '/api/dispensary/statistics' && method === 'GET') {
      return handleStatistics();
    }
    if (path === '/api/dispensary/inflammatory-synergy' && method === 'POST') {
      return handleInflammatorySynergy(request);
    }
    if (path === '/api/dispensary/adjuvants/optimize' && method === 'POST') {
      return handleAdjuvantOptimize(request);
    }

    // Trade-secret stubs (basic existence checks used by integration tests)
    const mappings = [
      ['/api/v1/chempath/analyze', 200],
      ['/api/v1/chempath/validate-smiles', 200],
      ['/api/v1/toxpath/assess', 200],
      ['/api/v1/toxpath/statistics', 200],
      ['/api/v1/regpath/strategy', 200],
      ['/api/v1/regpath/pathways', 200],
      ['/api/genomepath/tk-to-genomic', 200],
      ['/api/genomepath/statistics', 200],
      ['/api/biopath/validate', 200],
      ['/api/biopath/validate-from-studies', 200],
      ['/api/biopath/statistics', 200],
      ['/api/clinpath/optimize', 200],
      ['/api/clinpath/predict-approval', 200],
      ['/api/clinpath/jurisdiction-sequence', 200],
      ['/api/clinpath/statistics', 200],
      ['/api/dimers/predict/homodimer', 200],
      ['/api/dimers/predict/heterodimer', 200],
    ];

    for (const [p, code] of mappings) {
      if (path === p) {
        return json({ status: 'stub', path: p });
      }
    }

    return new Response('Not Found', { status: 404 });
  }
};
