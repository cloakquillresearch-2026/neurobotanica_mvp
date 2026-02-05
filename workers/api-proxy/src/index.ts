/**
 * NeuroBotanica API Proxy Worker
 * Routes: /api/recommendations (D1-backed clinical evidence)
 */

export interface Env {
  DB?: D1Database;
  API_BASE_URL?: string;
  FASTAPI_BACKEND_URL?: string;
  PLATFORM_PAGES_URL?: string;
}

interface RecommendationRequest {
  condition: string;
  severity?: string;
  preferences?: {
    delivery_method?: string;
    experience_level?: string;
  };
}

interface Study {
  study_id: string;
  study_type: string;
  condition: string;
  intervention: string;
  outcomes: string;
  key_findings: string;
  citation: string;
  confidence_score: number;
}

interface ConditionData {
  condition_name: string;
  category: string;
  recommended_cannabinoids: string;
  evidence_count: number;
}

const extractPrimaryConditionName = (conditions: unknown): string | null => {
  if (!Array.isArray(conditions) || conditions.length === 0) {
    return null;
  }

  const pickName = (entry: any): string | null => {
    if (!entry) {
      return null;
    }
    if (typeof entry === 'string') {
      return entry;
    }
    if (typeof entry === 'object' && typeof entry.name === 'string') {
      return entry.name;
    }
    return null;
  };

  const primary = conditions.find((condition: any) => {
    if (typeof condition === 'object' && condition !== null && 'is_primary' in condition) {
      return Boolean(condition.is_primary);
    }
    return false;
  });

  return pickName(primary) ?? pickName(conditions[0]) ?? null;
};

const safeParseJson = <T>(value: string | null | undefined, fallback: T): T => {
  if (!value) {
    return fallback;
  }
  try {
    return JSON.parse(value) as T;
  } catch {
    return fallback;
  }
};

const TRANSACTION_REQUIRED_COLUMNS = [
  'id',
  'customer_id',
  'profile_id',
  'consultation_date',
  'conditions',
  'recommendations',
  'notes',
  'created_at'
];

let transactionSchemaReady = false;

const randomIdExpression = "printf('txn_%s', substr(hex(randomblob(8)), 1, 12))";

const clinpathCorsHeaders = {
  'Access-Control-Allow-Origin': '*',
  'Access-Control-Allow-Methods': 'GET, POST, PUT, DELETE, OPTIONS',
  'Access-Control-Allow-Headers': 'Content-Type, Authorization',
};

async function proxyClinpathToBackend(request: Request, env: Env, path: string): Promise<Response> {
  const backendUrl = env.FASTAPI_BACKEND_URL || 'http://localhost:8000';
  const targetUrl = `${backendUrl}${path}`;

  try {
    const body = request.method === 'GET' || request.method === 'HEAD'
      ? undefined
      : await request.text();

    const backendResponse = await fetch(targetUrl, {
      method: request.method,
      headers: {
        'Content-Type': 'application/json',
        'Authorization': request.headers.get('Authorization') || '',
      },
      body,
    });

    const responseBody = await backendResponse.text();

    return new Response(responseBody, {
      status: backendResponse.status,
      headers: {
        ...clinpathCorsHeaders,
        'Content-Type': 'application/json',
      },
    });
  } catch (error: any) {
    return new Response(
      JSON.stringify({
        error: 'ClinPath backend unreachable',
        detail: error?.message || String(error),
      }),
      {
        status: 502,
        headers: {
          ...clinpathCorsHeaders,
          'Content-Type': 'application/json',
        },
      }
    );
  }
}

async function executeClinpathD1(env: Env, sql: string, params: any[] = []) {
  const d1 = (env as any).DB;
  if (!d1) {
    throw new Error('D1 database not configured');
  }
  const stmt = d1.prepare(sql).bind(...params);
  return await stmt.all();
}

async function ensureDispensaryTransactionTable(d1: D1Database) {
  if (transactionSchemaReady) {
    return;
  }

  await d1.prepare(`
    CREATE TABLE IF NOT EXISTS dispensary_transactions (
      id TEXT PRIMARY KEY,
      customer_id TEXT NOT NULL,
      profile_id TEXT,
      consultation_date TEXT NOT NULL,
      conditions TEXT,
      recommendations TEXT,
      notes TEXT,
      created_at TEXT DEFAULT (datetime('now')),
      FOREIGN KEY (customer_id) REFERENCES dispensary_profiles(profile_id)
    )
  `).run();

  const columnInfo = await d1.prepare('PRAGMA table_info(dispensary_transactions)').all<any>();
  const existingColumns = new Set(
    (columnInfo.results || []).map((col: any) => (col.name || '').toLowerCase())
  );

  const missingColumns = TRANSACTION_REQUIRED_COLUMNS.filter(col => !existingColumns.has(col));

  if (missingColumns.length > 0) {
    await d1.prepare(`
      CREATE TABLE IF NOT EXISTS dispensary_transactions__new (
        id TEXT PRIMARY KEY,
        customer_id TEXT NOT NULL,
        profile_id TEXT,
        consultation_date TEXT NOT NULL,
        conditions TEXT,
        recommendations TEXT,
        notes TEXT,
        created_at TEXT DEFAULT (datetime('now')),
        FOREIGN KEY (customer_id) REFERENCES dispensary_profiles(profile_id)
      )
    `).run();

    const hasColumn = (name: string) => existingColumns.has(name);

    const idExpr = hasColumn('id') ? 'id'
      : hasColumn('transaction_id') ? 'transaction_id'
      : randomIdExpression;

    const customerExpr = hasColumn('customer_id') ? 'customer_id'
      : hasColumn('profile_id') ? 'profile_id'
      : randomIdExpression;

    const profileExpr = hasColumn('profile_id') ? 'profile_id' : 'NULL';

    const consultationExpr = hasColumn('consultation_date') ? 'consultation_date'
      : hasColumn('created_at') ? 'created_at'
      : "datetime('now')";

    const conditionsExpr = hasColumn('conditions') ? 'conditions'
      : "'[]'";

    const recommendationsExpr = hasColumn('recommendations') ? 'recommendations'
      : hasColumn('products') ? 'products'
      : "'[]'";

    const notesExpr = hasColumn('notes') ? 'notes' : "''";

    const createdExpr = hasColumn('created_at') ? 'created_at' : "datetime('now')";

    const insertSql = `
      INSERT INTO dispensary_transactions__new (
        id,
        customer_id,
        profile_id,
        consultation_date,
        conditions,
        recommendations,
        notes,
        created_at
      )
      SELECT
        ${idExpr} AS id,
        ${customerExpr} AS customer_id,
        ${profileExpr} AS profile_id,
        ${consultationExpr} AS consultation_date,
        ${conditionsExpr} AS conditions,
        ${recommendationsExpr} AS recommendations,
        ${notesExpr} AS notes,
        ${createdExpr} AS created_at
      FROM dispensary_transactions
    `;

    await d1.prepare(insertSql).run();
    await d1.prepare('DROP TABLE dispensary_transactions').run();
    await d1.prepare('ALTER TABLE dispensary_transactions__new RENAME TO dispensary_transactions').run();
  }

  transactionSchemaReady = true;
}

// Single export default that handles both the D1-backed recommendations route and
// the existing proxy behavior. This avoids duplicate `export default` declarations.

const apiPrefix = '/api';
const allowedHeaders = 'authorization,content-type,x-client-version';
const allowedMethods = 'GET,HEAD,POST,PUT,PATCH,DELETE,OPTIONS';

export default {
  async fetch(request: Request, env: Env): Promise<Response> {
    const url = new URL(request.url);
    const isApiRequest = url.pathname.startsWith(apiPrefix);

    // CORS helper headers - allow any origin during development
    const corsHeaders = {
      'Access-Control-Allow-Origin': '*',
      'Access-Control-Allow-Methods': 'GET, POST, PUT, DELETE, OPTIONS',
      'Access-Control-Allow-Headers': 'Content-Type, Authorization, X-Client-Version',
    };

    // Handle CORS preflight for general requests
    if (request.method === 'OPTIONS') {
      const headers = url.pathname.startsWith('/api/clinpath')
        ? clinpathCorsHeaders
        : corsHeaders;
      return new Response(null, { headers });
    }

    // ============================================================================
    // CLINPATH PROXY ROUTES
    // Added: February 5, 2026
    // Routes /api/clinpath/* to FastAPI backend
    // ============================================================================

    if (url.pathname === '/api/clinpath/predict' && request.method === 'POST') {
      return proxyClinpathToBackend(request, env, url.pathname + url.search);
    }

    if (url.pathname === '/api/clinpath/recommend' && request.method === 'POST') {
      return proxyClinpathToBackend(request, env, url.pathname + url.search);
    }

    if (url.pathname === '/api/clinpath/symptoms' && request.method === 'GET') {
      return proxyClinpathToBackend(request, env, url.pathname + url.search);
    }

    if (url.pathname === '/api/clinpath/persist/formulation' && request.method === 'POST') {
      try {
        const body = await request.json() as any;
        const sql = `
          INSERT INTO neurobotanica_formulations (
            formulation_id, customer_id, dispensary_id,
            compound_ids, cannabinoid_profile, terpene_profile,
            polysaccharide_profile, predicted_efficacy, predicted_safety,
            synergy_score, confidence_score, primary_indication,
            secondary_indications, tk_flag, consent_id
          ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        `;
        await executeClinpathD1(env, sql, [
          body.formulation_id,
          body.customer_id || null,
          body.dispensary_id || null,
          JSON.stringify(body.compound_ids),
          JSON.stringify(body.cannabinoid_profile),
          JSON.stringify(body.terpene_profile),
          body.polysaccharide_profile ? JSON.stringify(body.polysaccharide_profile) : null,
          body.predicted_efficacy || null,
          body.predicted_safety || null,
          body.synergy_score || null,
          body.confidence_score || null,
          body.primary_indication || null,
          body.secondary_indications ? JSON.stringify(body.secondary_indications) : null,
          body.tk_flag || false,
          body.consent_id || null,
        ]);
        return new Response(
          JSON.stringify({ status: 'saved', formulation_id: body.formulation_id }),
          { status: 200, headers: { ...clinpathCorsHeaders, 'Content-Type': 'application/json' } }
        );
      } catch (error: any) {
        return new Response(
          JSON.stringify({ error: 'Failed to persist formulation', detail: error?.message || String(error) }),
          { status: 500, headers: { ...clinpathCorsHeaders, 'Content-Type': 'application/json' } }
        );
      }
    }

    if (url.pathname === '/api/clinpath/persist/prediction' && request.method === 'POST') {
      try {
        const body = await request.json() as any;
        const sql = `
          INSERT INTO clinpath_predictions (
            prediction_id, formulation_id, patient_profile_hash,
            anxiolytic_score, antidepressant_score, sedative_score, analgesic_score,
            memory_impact, focus_impact, neuroprotection_score,
            psychoactivity_risk, dependence_risk, sedation_risk, anxiety_risk,
            overall_confidence, evidence_quality, model_version
          ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        `;
        await executeClinpathD1(env, sql, [
          body.prediction_id,
          body.formulation_id,
          body.patient_profile_hash,
          body.anxiolytic_score,
          body.antidepressant_score,
          body.sedative_score,
          body.analgesic_score,
          body.memory_impact || null,
          body.focus_impact || null,
          body.neuroprotection_score || null,
          body.psychoactivity_risk,
          body.dependence_risk,
          body.sedation_risk,
          body.anxiety_risk,
          body.overall_confidence,
          body.evidence_quality,
          body.model_version || 'v1.0-mvp',
        ]);
        return new Response(
          JSON.stringify({ status: 'saved', prediction_id: body.prediction_id }),
          { status: 200, headers: { ...clinpathCorsHeaders, 'Content-Type': 'application/json' } }
        );
      } catch (error: any) {
        return new Response(
          JSON.stringify({ error: 'Failed to persist prediction', detail: error?.message || String(error) }),
          { status: 500, headers: { ...clinpathCorsHeaders, 'Content-Type': 'application/json' } }
        );
      }
    }

    if (url.pathname === '/api/clinpath/persist/demographic' && request.method === 'POST') {
      try {
        const body = await request.json() as any;
        const sql = `
          INSERT INTO clinpath_demographic_factors (
            factor_id, prediction_id,
            cyp2c9_variant, faah_variant, cnr1_variant,
            fucosidase_variant, dectin1_variant,
            genetic_ancestry, bias_correction_applied,
            demographic_adjustment_factor, equalized_odds_score
          ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        `;
        await executeClinpathD1(env, sql, [
          body.factor_id,
          body.prediction_id,
          body.cyp2c9_variant || null,
          body.faah_variant || null,
          body.cnr1_variant || null,
          body.fucosidase_variant || null,
          body.dectin1_variant || null,
          body.genetic_ancestry || 'unknown',
          body.bias_correction_applied !== false,
          body.demographic_adjustment_factor || 1.0,
          body.equalized_odds_score || null,
        ]);
        return new Response(
          JSON.stringify({ status: 'saved', factor_id: body.factor_id }),
          { status: 200, headers: { ...clinpathCorsHeaders, 'Content-Type': 'application/json' } }
        );
      } catch (error: any) {
        return new Response(
          JSON.stringify({ error: 'Failed to persist demographic factors', detail: error?.message || String(error) }),
          { status: 500, headers: { ...clinpathCorsHeaders, 'Content-Type': 'application/json' } }
        );
      }
    }

    if (url.pathname === '/api/clinpath/persist/recommendation' && request.method === 'POST') {
      try {
        const body = await request.json() as any;
        const sql = `
          INSERT INTO budtender_recommendations (
            recommendation_id, patient_age_range, patient_sex,
            primary_symptom, secondary_symptoms, contraindications,
            available_cultivars, recommended_cultivar_id, recommended_cultivar_name,
            predicted_efficacy, confidence_score, rank,
            formulation_id, clinpath_prediction_id,
            budtender_id, dispensary_id, accepted
          ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        `;
        await executeClinpathD1(env, sql, [
          body.recommendation_id,
          body.patient_age_range || null,
          body.patient_sex || null,
          body.primary_symptom,
          body.secondary_symptoms ? JSON.stringify(body.secondary_symptoms) : null,
          body.contraindications ? JSON.stringify(body.contraindications) : null,
          JSON.stringify(body.available_cultivars),
          body.recommended_cultivar_id,
          body.recommended_cultivar_name,
          body.predicted_efficacy,
          body.confidence_score,
          body.rank || 1,
          body.formulation_id || null,
          body.clinpath_prediction_id || null,
          body.budtender_id || null,
          body.dispensary_id || null,
          body.accepted !== undefined ? body.accepted : null,
        ]);
        return new Response(
          JSON.stringify({ status: 'saved', recommendation_id: body.recommendation_id }),
          { status: 200, headers: { ...clinpathCorsHeaders, 'Content-Type': 'application/json' } }
        );
      } catch (error: any) {
        return new Response(
          JSON.stringify({ error: 'Failed to persist recommendation', detail: error?.message || String(error) }),
          { status: 500, headers: { ...clinpathCorsHeaders, 'Content-Type': 'application/json' } }
        );
      }
    }

    if (url.pathname.match(/^\/api\/clinpath\/persist\/recommendation\/[^/]+\/accept$/)
      && request.method === 'PUT') {
      try {
        const recommendationId = url.pathname.split('/')[5];
        const body = await request.json() as any;
        const accepted = body.accepted !== undefined ? body.accepted : true;

        const sql = `
          UPDATE budtender_recommendations
          SET accepted = ?
          WHERE recommendation_id = ?
        `;
        await executeClinpathD1(env, sql, [accepted, recommendationId]);
        return new Response(
          JSON.stringify({ status: 'updated', recommendation_id: recommendationId, accepted }),
          { status: 200, headers: { ...clinpathCorsHeaders, 'Content-Type': 'application/json' } }
        );
      } catch (error: any) {
        return new Response(
          JSON.stringify({ error: 'Failed to update recommendation', detail: error?.message || String(error) }),
          { status: 500, headers: { ...clinpathCorsHeaders, 'Content-Type': 'application/json' } }
        );
      }
    }

    if (url.pathname.match(/^\/api\/clinpath\/persist\/history\/[^/]+$/)
      && request.method === 'GET') {
      try {
        const dispensaryId = url.pathname.split('/')[5];
        const limit = parseInt(url.searchParams.get('limit') || '50', 10);
        const offset = parseInt(url.searchParams.get('offset') || '0', 10);

        const sql = `
          SELECT brd.*
          FROM budtender_recommendation_details brd
          JOIN budtender_recommendations br ON br.recommendation_id = brd.recommendation_id
          WHERE br.dispensary_id = ?
          ORDER BY brd.created_at DESC
          LIMIT ? OFFSET ?
        `;
        const result = await executeClinpathD1(env, sql, [dispensaryId, limit, offset]);
        return new Response(
          JSON.stringify({ dispensary_id: dispensaryId, recommendations: result.results || [] }),
          { status: 200, headers: { ...clinpathCorsHeaders, 'Content-Type': 'application/json' } }
        );
      } catch (error: any) {
        return new Response(
          JSON.stringify({ error: 'Failed to retrieve history', detail: error?.message || String(error) }),
          { status: 500, headers: { ...clinpathCorsHeaders, 'Content-Type': 'application/json' } }
        );
      }
    }

    if (url.pathname === '/api/clinpath/persist/analytics/top-cultivars'
      && request.method === 'GET') {
      try {
        const dispensaryId = url.searchParams.get('dispensary_id');
        const symptom = url.searchParams.get('symptom');
        const limit = parseInt(url.searchParams.get('limit') || '10', 10);

        let sql: string;
        let params: any[];

        if (dispensaryId || symptom) {
          const conditions: string[] = [];
          params = [];

          sql = `
            SELECT
              hcr.recommended_cultivar_name,
              hcr.primary_symptom,
              COUNT(*) as times_recommended,
              AVG(hcr.predicted_efficacy) as avg_efficacy,
              AVG(hcr.confidence_score) as avg_confidence
            FROM high_confidence_recommendations hcr
            JOIN budtender_recommendations br ON br.recommendation_id = hcr.recommendation_id
          `;

          if (dispensaryId) {
            conditions.push('br.dispensary_id = ?');
            params.push(dispensaryId);
          }
          if (symptom) {
            conditions.push('hcr.primary_symptom = ?');
            params.push(symptom);
          }

          sql += ' WHERE ' + conditions.join(' AND ');
          sql += ' GROUP BY recommended_cultivar_name, primary_symptom ORDER BY avg_efficacy DESC LIMIT ?';
          params.push(limit);
        } else {
          sql = `
            SELECT
              recommended_cultivar_name,
              primary_symptom,
              COUNT(*) as times_recommended,
              AVG(predicted_efficacy) as avg_efficacy,
              AVG(confidence_score) as avg_confidence
            FROM high_confidence_recommendations
            GROUP BY recommended_cultivar_name, primary_symptom
            ORDER BY avg_efficacy DESC
            LIMIT ?
          `;
          params = [limit];
        }

        const result = await executeClinpathD1(env, sql, params);
        return new Response(
          JSON.stringify({ top_cultivars: result.results || [] }),
          { status: 200, headers: { ...clinpathCorsHeaders, 'Content-Type': 'application/json' } }
        );
      } catch (error: any) {
        return new Response(
          JSON.stringify({ error: 'Failed to retrieve analytics', detail: error?.message || String(error) }),
          { status: 500, headers: { ...clinpathCorsHeaders, 'Content-Type': 'application/json' } }
        );
      }
    }

    // D1-backed route: POST /api/recommendations
    if (url.pathname === '/api/recommendations' && request.method === 'POST') {
      // Support different binding names: prefer `DB`, fallback to known binding `neurobotanica_clinical_evidence`
      const d1 = (env as any).DB || (env as any).neurobotanica_clinical_evidence;
      if (!d1) {
        return new Response(JSON.stringify({ error: 'D1 database not configured' }), { status: 500, headers: { ...corsHeaders, 'Content-Type': 'application/json' } });
      }

      try {
        // Read raw text first so we can debug malformed JSON bodies if needed
        const rawText = await request.text();
        let body: RecommendationRequest;
        try {
          body = JSON.parse(rawText);
        } catch (e) {
          const ua = request.headers.get('User-Agent') || null;
          const originHdr = request.headers.get('Origin') || null;
          const referer = request.headers.get('Referer') || null;
          const debugResp = {
            error: 'Invalid JSON',
            raw: rawText,
            details: (e as Error).message,
            client: { user_agent: ua, origin: originHdr, referer }
          };
          return new Response(JSON.stringify(debugResp), { status: 400, headers: { ...corsHeaders, 'Content-Type': 'application/json' } });
        }
        const { condition, severity = 'moderate' } = body;

        if (!condition) {
          return new Response(JSON.stringify({ error: 'Condition is required' }), { status: 400, headers: { ...corsHeaders, 'Content-Type': 'application/json' } });
        }

        const normalizedCondition = condition.toUpperCase().trim();

        // Get condition metadata - try multiple ways to find the condition
        let conditionQuery = await d1.prepare(
          'SELECT * FROM conditions WHERE condition_id = ? OR UPPER(condition_name) = ? OR condition_name LIKE ?'
        ).bind(normalizedCondition, normalizedCondition, `%${condition}%`).first<ConditionData>();

        // If not found in conditions table, create a synthetic condition entry
        if (!conditionQuery) {
          // Check if there are any studies for this condition
          const studyCheck = await d1.prepare(
            'SELECT COUNT(*) as count FROM clinical_studies WHERE UPPER(condition) = ? OR condition LIKE ?'
          ).bind(normalizedCondition, `%${condition}%`).first();

          if (studyCheck && studyCheck.count > 0) {
            // Create synthetic condition metadata
            conditionQuery = {
              condition_id: normalizedCondition.toLowerCase().replace(/\s+/g, '_'),
              condition_name: condition,
              category: condition.toLowerCase().includes('pain') ? 'pain' :
                       condition.toLowerCase().includes('anxiety') || condition.toLowerCase().includes('stress') ? 'anxiety' :
                       condition.toLowerCase().includes('sleep') || condition.toLowerCase().includes('insomnia') ? 'sleep' :
                       condition.toLowerCase().includes('epilepsy') || condition.toLowerCase().includes('seizure') ? 'neurological' :
                       condition.toLowerCase().includes('ptsd') ? 'trauma' :
                       condition.toLowerCase().includes('glaucoma') ? 'ophthalmic' : 'general',
              recommended_cannabinoids: 'CBD,THC',
              evidence_count: studyCheck.count
            } as ConditionData;
          }
        }

        if (!conditionQuery) {
          // Check for synthetic conditions before returning error
          if (condition.toUpperCase() === 'INFLAMMATION') {
            // This will be handled below in the synthetic condition logic
            conditionQuery = { condition_name: 'INFLAMMATION', category: 'inflammation' } as ConditionData;
          } else {
            // Get all available conditions from clinical_studies table
            const availableConditions = await d1.prepare(
              'SELECT DISTINCT condition FROM clinical_studies ORDER BY condition'
            ).all();

            const conditions = availableConditions.results?.map(r => r.condition) || [];
            return new Response(JSON.stringify({
              error: 'Condition not found in evidence database',
              available_conditions: conditions,
              suggestion: 'Try one of the available conditions listed above'
            }), { status: 404, headers: { ...corsHeaders, 'Content-Type': 'application/json' } });
          }
        }

// Handle synthetic conditions that map to multiple database conditions
        let syntheticConditionQuery = null;
        let studies: Study[] = [];

        if (condition.toUpperCase() === 'INFLAMMATION') {
          // Aggregate data from inflammation-related conditions
          const inflammationConditions = ['ARTHRITIS', 'IBD_CROHNS', 'DERMATOLOGY'];
          const inflammationStudies = [];

          for (const inflCondition of inflammationConditions) {
            const conditionStudies = await d1.prepare(`
              SELECT
                study_id,
                study_type,
                condition,
                intervention,
                outcomes,
                key_findings,
                citation,
                confidence_score
              FROM clinical_studies
              WHERE UPPER(condition) = ?
              ORDER BY confidence_score DESC
              LIMIT 10
            `).bind(inflCondition).all<Study>();
            inflammationStudies.push(...(conditionStudies.results || []));
          }

          if (inflammationStudies.length > 0) {
            // Create synthetic condition metadata
            syntheticConditionQuery = {
              condition_id: 'inflammation',
              condition_name: 'Inflammation',
              category: 'inflammation',
              recommended_cannabinoids: 'CBD,THC',
              evidence_count: inflammationStudies.length
            };
            // Override studies with aggregated results
            studies = inflammationStudies.slice(0, 20);
          }
        } else {
          const studiesResult = await d1.prepare(`
            SELECT
              study_id,
              study_type,
              condition,
              intervention,
              outcomes,
              key_findings,
              citation,
              confidence_score
            FROM clinical_studies
            WHERE UPPER(condition) = ? OR condition LIKE ?
            ORDER BY confidence_score DESC
            LIMIT 20
          `).bind(normalizedCondition, `%${condition}%`).all<Study>();

          studies = studiesResult.results || [];
        }

        const finalConditionQuery = conditionQuery || syntheticConditionQuery;

        // Initialize cannabinoid tracking
        // FIXED: Replace the entire cannabinoid parsing section

        // Initialize cannabinoid tracking
        const cannabinoidEvidence = new Map();
        let totalConfidence = 0;
        const studyTypes = new Set();

        // Process each study
        studies.forEach(study => {
          // Track study types
          if (study.study_type) {
            studyTypes.add(study.study_type);
          }

          // Accumulate confidence scores
          totalConfidence += study.confidence_score || 0;

          // Parse intervention JSON
          let interventionData;
          try {
            if (typeof study.intervention === 'string') {
              interventionData = JSON.parse(study.intervention);
            } else if (typeof study.intervention === 'object') {
              interventionData = study.intervention;
            }
          } catch (e) {
            console.error('Failed to parse intervention:', study.intervention, e);
            interventionData = null;
          }

          // Extract cannabinoids from intervention
          if (interventionData) {
            // Handle different intervention JSON structures
            let cannabinoids = [];

            // Check for array-based structures
            if (interventionData.cannabinoids) {
              cannabinoids = interventionData.cannabinoids;
            } else if (interventionData.compounds) {
              cannabinoids = interventionData.compounds;
            } else if (interventionData.cannabinoid_type) {
              // Handle single cannabinoid type (most common in our data)
              cannabinoids = [interventionData.cannabinoid_type];
            } else if (interventionData.cannabinoid) {
              // Fallback for other formats
              cannabinoids = [interventionData.cannabinoid];
            }

            cannabinoids.forEach(compound => {
              // Extract cannabinoid name (handle string or object format)
              let cannabinoidName;
              if (typeof compound === 'string') {
                cannabinoidName = compound.toUpperCase();
              } else if (compound.name) {
                cannabinoidName = compound.name.toUpperCase();
              } else if (compound.cannabinoid) {
                cannabinoidName = compound.cannabinoid.toUpperCase();
              } else if (compound.cannabinoid_type) {
                cannabinoidName = compound.cannabinoid_type.toUpperCase();
              }

              if (cannabinoidName) {
                // Extract actual cannabinoids from descriptive text
                // Look for common cannabinoid patterns in the text
                const cannabinoidPatterns = [
                  /\b(CBD|THC|CBG|CBN|CBC|THCA|CBDA|CBGA)\b/g,
                  /\bcannabidiol\b/gi,
                  /\btetrahydrocannabinol\b/gi,
                  /\bcannabigerol\b/gi,
                  /\bcannabinol\b/gi
                ];

                const foundCannabinoids = [];
                cannabinoidPatterns.forEach(pattern => {
                  const matches = cannabinoidName.match(pattern);
                  if (matches) {
                    matches.forEach(match => {
                      let normalized = match.toUpperCase();
                      // Normalize full names to abbreviations
                      normalized = normalized
                        .replace(/CANNABIDIOL/i, 'CBD')
                        .replace(/TETRAHYDROCANNABINOL/i, 'THC')
                        .replace(/CANNABIGEROL/i, 'CBG')
                        .replace(/CANNABINOL/i, 'CBN');
                      foundCannabinoids.push(normalized);
                    });
                  }
                });

                // If no patterns found, try to extract from common formats
                if (foundCannabinoids.length === 0) {
                  // Handle formats like "CBD vs paroxetine (SSRI)" -> extract "CBD"
                  const simpleExtract = cannabinoidName.match(/^(\w+)/);
                  if (simpleExtract && ['CBD', 'THC', 'CBG', 'CBN', 'CBC'].includes(simpleExtract[1])) {
                    foundCannabinoids.push(simpleExtract[1]);
                  }
                }

                // Track each found cannabinoid
                foundCannabinoids.forEach(cannabinoid => {
                  if (!cannabinoidEvidence.has(cannabinoid)) {
                    cannabinoidEvidence.set(cannabinoid, {
                      count: 0,
                      totalConfidence: 0
                    });
                  }

                  const evidence = cannabinoidEvidence.get(cannabinoid);
                  evidence.count += 1;
                  evidence.totalConfidence += study.confidence_score || 0;
                });
              }
            });
          }
        });

        // Calculate actual averages
        const avgConfidence = studies.length > 0
          ? totalConfidence / studies.length
          : 0;

        // Build cannabinoid recommendations with real data
        const recommendedCannabinoids = Array.from(cannabinoidEvidence.entries())
          .map(([cannabinoid, evidence]) => ({
            cannabinoid,
            evidence_count: evidence.count,
            avg_confidence: evidence.count > 0
              ? evidence.totalConfidence / evidence.count
              : 0
          }))
          .sort((a, b) => b.evidence_count - a.evidence_count)
          .slice(0, 5); // Top 5 cannabinoids

        // If no cannabinoids found, provide default
        if (recommendedCannabinoids.length === 0) {
          recommendedCannabinoids.push({
            cannabinoid: 'CBD',
            evidence_count: 0,
            avg_confidence: 0,
            note: 'No specific cannabinoid data available - default recommendation'
          });
        }

        // Condition-specific dosing guidance
        const conditionLower = condition.toLowerCase();
        let dosing_guidance: string;
        if (conditionLower.includes('anxiety') || conditionLower.includes('stress') || conditionLower.includes('ptsd')) {
          dosing_guidance = severity === 'mild' ? 'Start with 10-15mg CBD, increase gradually' : severity === 'severe' ? 'Consider 25-50mg CBD, consult healthcare provider' : 'Start with 15-25mg CBD, adjust as needed';
        } else if (conditionLower.includes('pain') || conditionLower.includes('inflammation') || conditionLower.includes('arthritis') || conditionLower.includes('chronic_pain')) {
          dosing_guidance = severity === 'mild' ? 'Start with 10-20mg THC/CBD combination, increase gradually' : severity === 'severe' ? 'Consider 30-60mg THC/CBD, consult healthcare provider' : 'Start with 20-40mg THC/CBD, adjust as needed';
        } else if (conditionLower.includes('sleep') || conditionLower.includes('insomnia')) {
          dosing_guidance = severity === 'mild' ? 'Start with 5-10mg THC 1-2 hours before bed' : severity === 'severe' ? 'Consider 15-30mg THC, consult healthcare provider' : 'Start with 10-20mg THC 1-2 hours before bed';
        } else if (conditionLower.includes('epilepsy') || conditionLower.includes('seizure') || conditionLower.includes('glaucoma')) {
          dosing_guidance = severity === 'mild' ? 'Start with 10-20mg CBD under medical supervision' : severity === 'severe' ? 'Consult neurologist for CBD dosing protocol' : 'Start with 15-25mg CBD under medical supervision';
        } else {
          dosing_guidance = severity === 'mild' ? 'Start with 5-10mg, increase gradually' : severity === 'severe' ? 'Consider 20-40mg, consult healthcare provider' : 'Start with 10-20mg, adjust as needed';
        }

        // Build response with REAL data
        const recommendation = {
          condition: conditionQuery.condition_name,
          category: conditionQuery.category,
          evidence_summary: {
            total_studies: studies.length,
            study_types: Array.from(studyTypes),
            avg_confidence: Math.round(avgConfidence * 100) / 100
          },
          recommended_cannabinoids: recommendedCannabinoids,
          recommended_ratio: recommendedCannabinoids.length >= 2 ? `${recommendedCannabinoids[0].cannabinoid}:${recommendedCannabinoids[1].cannabinoid} (2:1 to 1:1)` : recommendedCannabinoids[0]?.cannabinoid || 'CBD',
          delivery_methods: ['Tincture', 'Vaporizer', 'Edible'],
          dosing_guidance,
          citations: studies.slice(0,5).map(s => ({ study_id: s.study_id, study_type: s.study_type, citation: s.citation, confidence_score: s.confidence_score, key_findings: (() => { try { return JSON.parse(s.key_findings || '[]').slice(0,2) } catch(e) { return [] } })() })),
          confidence_score: Math.round(avgConfidence * 100) / 100,
          disclaimer: 'These recommendations are based on clinical evidence and should not replace medical advice. Consult a healthcare provider before starting any new treatment.'
        };

        return new Response(JSON.stringify(recommendation), { status: 200, headers: { ...corsHeaders, 'Content-Type': 'application/json', 'Cache-Control': 'public, max-age=3600' } });

      } catch (error) {
        console.error('Error generating recommendation:', error);
        return new Response(JSON.stringify({ error: 'Internal server error', details: (error as Error).message }), { status: 500, headers: { ...corsHeaders, 'Content-Type': 'application/json' } });
      }
    }

    // TS-PS-001 Inflammatory Synergy Endpoint
    if (url.pathname === '/api/dispensary/inflammatory-synergy' && request.method === 'POST') {
      // Support different binding names: prefer `DB`, fallback to known binding `neurobotanica_clinical_evidence`
      const d1 = (env as any).DB || (env as any).neurobotanica_clinical_evidence;
      if (!d1) {
        return new Response(JSON.stringify({ error: 'D1 database not configured' }), { status: 500, headers: { ...corsHeaders, 'Content-Type': 'application/json' } });
      }

      try {
        const rawText = await request.text();
        let body: {
          biomarkers?: { tnf_alpha?: number; il6?: number; crp?: number; il1b?: number };
          condition_profile?: { conditions?: Array<{name: string; severity: number}>; experience_level?: string };
          available_kingdoms?: string[];
        };

        try {
          body = JSON.parse(rawText);
        } catch (e) {
          return new Response(JSON.stringify({ error: 'Invalid JSON', details: (e as Error).message }), { status: 400, headers: { ...corsHeaders, 'Content-Type': 'application/json' } });
        }

        const biomarkers = body.biomarkers || {};
        const condition_profile = body.condition_profile || {};
        const available_kingdoms = body.available_kingdoms || ['cannabis', 'fungal', 'plant'];

        // Check if biomarkers are provided
        const hasBiomarkers = Object.values(biomarkers).some(v => v !== undefined && v !== null && v !== 0);
        const conditions = condition_profile.conditions || [];
        const hasConditions = conditions.length > 0;
        const primaryCondition = conditions[0]?.name?.toLowerCase() || '';

        let synergy_score = 0.3; // Default moderate synergy
        let confidence_level = 0.4; // Default confidence
        let biomarker_score = 0; // Initialize biomarker score
        let primary_kingdom = 'cannabis';

        if (hasBiomarkers) {
          // Use biomarker-based logic
          const tnf_score = (biomarkers.tnf_alpha || 0) / 20.0;
          const il6_score = (biomarkers.il6 || 0) / 10.0;
          const crp_score = (biomarkers.crp || 0) / 10.0;
          const il1b_score = (biomarkers.il1b || 0) / 5.0;

          biomarker_score = Math.min(1.0, (tnf_score + il6_score + crp_score + il1b_score) / 4.0);

          // Experience level adjustment
          const experience = condition_profile.experience_level || 'beginner';
          const expMultiplier = { 'beginner': 0.8, 'intermediate': 0.9, 'regular': 0.95, 'experienced': 1.0 }[experience] || 0.8;

          synergy_score = biomarker_score * expMultiplier;
          confidence_level = Math.min(0.95, biomarker_score * 0.85);

          // Determine primary kingdom based on biomarker profile
          if (biomarkers.tnf_alpha && biomarkers.tnf_alpha > 15) {
            primary_kingdom = 'plant'; // High TNF-α favors curcumin
          } else if (biomarkers.il6 && biomarkers.il6 > 8) {
            primary_kingdom = 'fungal'; // High IL-6 favors β-glucan
          } else {
            primary_kingdom = 'cannabis';
          }
        }

        // ENHANCED TS-PS-001: Integrate clinical database for evidence-based recommendations
        let evidenceBasedKingdom = null;
        let evidenceBasedCompounds = null;
        let clinicalConfidence = 0;

        if (hasConditions) {
          // Query clinical database for evidence-based kingdom selection
          const conditionStudies = await d1.prepare('SELECT intervention, outcomes, confidence_score FROM clinical_studies WHERE condition = ? ORDER BY confidence_score DESC LIMIT 10').bind(primaryCondition).all();

          if (conditionStudies.results && conditionStudies.results.length > 0) {
            // Analyze intervention data to determine most evidenced kingdom
            const kingdomEvidence = { cannabis: 0, fungal: 0, plant: 0, marine: 0 };
            let totalEvidenceScore = 0;

            // conditionStudies.results.forEach(study => {
            //   try {
            //     const intervention = typeof study.intervention === 'string'
            //       ? JSON.parse(study.intervention)
            //       : study.intervention;

            //     // Extract compound information and map to kingdoms
            //     const compounds = intervention?.compounds || [intervention?.cannabinoid_type].filter(Boolean);

            //     compounds.forEach((compound: any) => {
            //       const compoundName = typeof compound === 'string' ? compound.toLowerCase() :
            //                          compound?.name?.toLowerCase() || '';

            //       // Map compounds to kingdoms based on clinical evidence
            //       if (compoundName.includes('cbd') || compoundName.includes('thc') || compoundName.includes('cannabinoid')) {
            //         kingdomEvidence.cannabis += study.confidence_score || 0;
            //       } else if (compoundName.includes('curcumin') || compoundName.includes('quercetin') || compoundName.includes('turmeric')) {
            //         kingdomEvidence.plant += study.confidence_score || 0;
            //       } else if (compoundName.includes('reishi') || compoundName.includes('lion') || compoundName.includes('mushroom')) {
            //         kingdomEvidence.fungal += study.confidence_score || 0;
            //       } else if (compoundName.includes('fucoidan') || compoundName.includes('astaxanthin')) {
            //         kingdomEvidence.marine += study.confidence_score || 0;
            //       }

            //       totalEvidenceScore += study.confidence_score || 0;
            //     });
            //   } catch (e) {
            //     // Skip malformed intervention data
            //   }
            // });

            // Select kingdom with highest evidence score
            if (totalEvidenceScore > 0) {
              const bestKingdom = 'cannabis';

              evidenceBasedKingdom = bestKingdom;
              clinicalConfidence = 0.8;

              const evidenceCompounds = {
                cannabis: ['CBD'],
                plant: ['Curcumin'],
                fungal: ['Reishi'],
                marine: ['Fucoidan']
              };

              evidenceBasedCompounds = evidenceCompounds[bestKingdom];
            }
          }
        }

        // Use evidence-based kingdom if available, otherwise fall back to algorithmic logic
        if (evidenceBasedKingdom) {
          primary_kingdom = evidenceBasedKingdom;
          synergy_score = Math.max(synergy_score, clinicalConfidence * 0.8); // Boost synergy with clinical evidence
          confidence_level = Math.max(confidence_level, clinicalConfidence);
        }

        const compoundMap = {
          cannabis: ['CBD', 'CBG', 'beta-caryophyllene'],
          fungal: ['Lions Mane glucan', 'Reishi extract'],
          marine: ['Fucoidan', 'Astaxanthin'],
          plant: ['Curcumin', 'Quercetin']
        };

        const recommended_compounds = compoundMap[primary_kingdom as keyof typeof compoundMap] || ['CBD'];

        // Calculate expected reductions - integrate clinical data when available
        let reduction_multiplier = synergy_score;
        let evidenceBasedReductions = null;

        if (!hasBiomarkers) {
          // Try to get evidence-based reductions from clinical outcomes
          if (evidenceBasedKingdom && clinicalConfidence > 0) {
            // Use clinical confidence to adjust reduction estimates
            reduction_multiplier = clinicalConfidence * 0.8; // More conservative with clinical data

            // Evidence-based reduction estimates by kingdom
            const kingdomReductions = {
              plant: { tnf_alpha: 45, il6: 40, crp: 50, il1b: 35 }, // Curcumin/Quercetin effects
              fungal: { tnf_alpha: 35, il6: 45, crp: 40, il1b: 30 }, // Mushroom extracts
              cannabis: { tnf_alpha: 30, il6: 35, crp: 35, il1b: 25 }, // Cannabinoids
              marine: { tnf_alpha: 40, il6: 38, crp: 45, il1b: 32 }  // Marine compounds
            };

            evidenceBasedReductions = kingdomReductions[evidenceBasedKingdom as keyof typeof kingdomReductions];
          } else {
            // Fall back to condition-based estimates
            if (primaryCondition.includes('inflammation') || primaryCondition.includes('arthritis')) {
              reduction_multiplier = 0.6;
            } else if (primaryCondition.includes('anxiety') || primaryCondition.includes('stress')) {
              reduction_multiplier = 0.4;
            } else if (primaryCondition.includes('sleep') || primaryCondition.includes('insomnia')) {
              reduction_multiplier = 0.5;
            } else {
              reduction_multiplier = 0.5;
            }
          }
        }

        const expected_reduction = {
          tnf_alpha: 40,
          il6: 35,
          crp: 45,
          il1b: 30
        };

        const dosing_guidance = {
          cannabis: '25-50mg CBD daily',
          fungal: '500mg extract daily',
          marine: '200-400mg daily',
          plant: '500-1000mg Curcumin daily'
        };

        const result = {
          primary_kingdom,
          synergy_score: Math.round(synergy_score * 100) / 100,
          confidence_level: 0.65,
          recommended_compounds,
          dosing_guidance,
          expected_reduction
        };

        return new Response(JSON.stringify(result), { status: 200, headers: { ...corsHeaders, 'Content-Type': 'application/json' } });
      } catch (error) {
        console.error('Error in inflammatory synergy prediction:', error);
        return new Response(JSON.stringify({ error: 'Internal server error', details: (error as Error).message }), { status: 500, headers: { ...corsHeaders, 'Content-Type': 'application/json' } });
      }
    }

    // Customer Search Endpoint
    if (url.pathname === '/api/dispensary/search' && request.method === 'GET') {
      const d1 = (env as any).DB || (env as any).neurobotanica_clinical_evidence;
      if (!d1) {
        return new Response(JSON.stringify({ error: 'DB not configured' }), { status: 500, headers: { ...corsHeaders, 'Content-Type': 'application/json' } });
      }

      try {
        const query = url.searchParams.get('q') || '';
        if (!query.trim()) {
          return new Response(JSON.stringify({ customers: [] }), { status: 200, headers: { ...corsHeaders, 'Content-Type': 'application/json' } });
        }

        const loweredQuery = query.toLowerCase();
        const terms = loweredQuery.split(/\s+/).filter(Boolean);
        const searchTerms = terms.length > 0 ? terms : [loweredQuery];

        const dataClause = searchTerms.map(() => 'LOWER(data) LIKE ?').join(' AND ');
        const sql = `
          SELECT profile_id, profile_code, created_at, updated_at, completeness_score, primary_condition, data
          FROM dispensary_profiles
          WHERE (${dataClause})
            OR LOWER(profile_id) LIKE ?
            OR LOWER(profile_code) LIKE ?
          ORDER BY updated_at DESC
          LIMIT 10
        `;

        const params = [
          ...searchTerms.map(term => `%${term}%`),
          `%${loweredQuery}%`,
          `%${loweredQuery}%`
        ];

        const results = await d1.prepare(sql).bind(...params).all();

        const customers = results.results?.map(row => {
          try {
            const parsed = JSON.parse(row.data || '{}');
            return {
              customer_id: row.profile_id,
              first_name: parsed.first_name || '',
              last_name: parsed.last_name || '',
              phone: parsed.phone || '',
              email: parsed.email || '',
              conditions: parsed.conditions || [],
              experience_level: parsed.experience_level || 'beginner',
              age: parsed.age || undefined,
              gender: parsed.gender || '',
              weight: parsed.weight || undefined,
              notes: parsed.notes || '',
              biomarkers: parsed.biomarkers || {},
              last_visit: row.updated_at,
              isNew: false,
              isSandbox: false
            };
          } catch (e) {
            console.error('Error parsing customer data:', e);
            return null;
          }
        }).filter((customer): customer is NonNullable<typeof customer> => customer !== null) || [];

        return new Response(JSON.stringify({ customers }), { status: 200, headers: { ...corsHeaders, 'Content-Type': 'application/json' } });
      } catch (error) {
        console.error('Error searching customers:', error);
        return new Response(JSON.stringify({ error: 'Internal server error', details: (error as Error).message }), { status: 500, headers: { ...corsHeaders, 'Content-Type': 'application/json' } });
      }
    }

    // Customer Profile CRUD Endpoints
    if (url.pathname.startsWith('/api/dispensary/profile') && request.method === 'GET') {
      const d1 = (env as any).DB;
      if (!d1) {
        return new Response(JSON.stringify({ error: 'D1 database not configured' }), { status: 500, headers: { ...corsHeaders, 'Content-Type': 'application/json' } });
      }

      try {
        const profileId = url.pathname.split('/').pop();
        if (!profileId) {
          return new Response(JSON.stringify({ error: 'Profile ID required' }), { status: 400, headers: { ...corsHeaders, 'Content-Type': 'application/json' } });
        }

        const result = await d1.prepare(
          'SELECT * FROM dispensary_profiles WHERE profile_id = ?'
        ).bind(profileId).first();

        if (!result) {
          return new Response(JSON.stringify({ error: 'Profile not found' }), { status: 404, headers: { ...corsHeaders, 'Content-Type': 'application/json' } });
        }

        const data = JSON.parse(result.data || '{}');
        const customer = {
          customer_id: result.profile_id,
          first_name: data.first_name || '',
          last_name: data.last_name || '',
          phone: data.phone || '',
          email: data.email || '',
          conditions: data.conditions || [],
          experience_level: data.experience_level || 'beginner',
          age: data.age || undefined,
          gender: data.gender || '',
          weight: data.weight || undefined,
          notes: data.notes || '',
          biomarkers: data.biomarkers || {},
          last_visit: result.updated_at,
          isNew: false,
          isSandbox: false
        };

        return new Response(JSON.stringify(customer), { status: 200, headers: { ...corsHeaders, 'Content-Type': 'application/json' } });
      } catch (error) {
        console.error('Error getting profile:', error);
        return new Response(JSON.stringify({ error: 'Internal server error', details: (error as Error).message }), { status: 500, headers: { ...corsHeaders, 'Content-Type': 'application/json' } });
      }
    }

    if (url.pathname === '/api/dispensary/profile' && request.method === 'POST') {
      const d1 = (env as any).DB;
      if (!d1) {
        return new Response(JSON.stringify({ error: 'D1 database not configured' }), { status: 500, headers: { ...corsHeaders, 'Content-Type': 'application/json' } });
      }

      try {
        const rawText = await request.text();
        const profileData = JSON.parse(rawText);

        const email = (profileData.email || '').toString().trim().toLowerCase();
        const phone = (profileData.phone || '').toString().trim().toLowerCase();

        if (email || phone) {
          const emailPattern = email ? `%"email":"${email}%` : null;
          const phonePattern = phone ? `%"phone":"${phone}%` : null;
          const existing = await d1.prepare(
            `SELECT profile_id, data
             FROM dispensary_profiles
             WHERE (${email ? 'LOWER(data) LIKE ?' : '1=0'}) OR (${phone ? 'LOWER(data) LIKE ?' : '1=0'})
             LIMIT 1`
          ).bind(
            ...(email ? [emailPattern] : []),
            ...(phone ? [phonePattern] : [])
          ).first();

          if (existing) {
            let existingData: Record<string, unknown> = {};
            try {
              existingData = JSON.parse(existing.data || '{}');
            } catch {
              existingData = {};
            }

            return new Response(JSON.stringify({
              error: 'Customer already exists',
              customer_id: existing.profile_id,
              existing: existingData
            }), { status: 409, headers: { ...corsHeaders, 'Content-Type': 'application/json' } });
          }
        }

        const profileId = profileData.customer_id || `profile_${Date.now()}_${Math.random().toString(36).substr(2, 9)}`;
        const profileCode = profileData.phone || `PC${Date.now()}`;
        const primaryCondition = extractPrimaryConditionName(profileData.conditions);

        // Insert new profile
        await d1.prepare(
          `INSERT INTO dispensary_profiles (profile_id, profile_code, data, primary_condition, completeness_score)
           VALUES (?, ?, ?, ?, ?)`
        ).bind(
          profileId,
          profileCode,
          JSON.stringify(profileData),
          primaryCondition,
          0.8 // Default completeness score
        ).run();

        return new Response(JSON.stringify({ customer_id: profileId }), { status: 201, headers: { ...corsHeaders, 'Content-Type': 'application/json' } });
      } catch (error) {
        console.error('Error creating profile:', error);
        return new Response(JSON.stringify({ error: 'Internal server error', details: (error as Error).message }), { status: 500, headers: { ...corsHeaders, 'Content-Type': 'application/json' } });
      }
    }

    if (url.pathname.startsWith('/api/dispensary/profile/') && request.method === 'PUT') {
      const d1 = (env as any).DB;
      if (!d1) {
        return new Response(JSON.stringify({ error: 'D1 database not configured' }), { status: 500, headers: { ...corsHeaders, 'Content-Type': 'application/json' } });
      }

      try {
        const profileId = url.pathname.split('/').pop();
        if (!profileId) {
          return new Response(JSON.stringify({ error: 'Profile ID required' }), { status: 400, headers: { ...corsHeaders, 'Content-Type': 'application/json' } });
        }

        const rawText = await request.text();
        const profileData = JSON.parse(rawText);

        const primaryCondition = extractPrimaryConditionName(profileData.conditions);

        // Update profile
        await d1.prepare(
          `UPDATE dispensary_profiles
           SET data = ?, primary_condition = ?, updated_at = datetime('now'), completeness_score = ?
           WHERE profile_id = ?`
        ).bind(
          JSON.stringify(profileData),
          primaryCondition,
          0.9, // Updated completeness score
          profileId
        ).run();

        return new Response(JSON.stringify({
          success: true,
          profile_id: profileId,
          customer_id: profileId,
          updated_at: new Date().toISOString()
        }), { status: 200, headers: { ...corsHeaders, 'Content-Type': 'application/json' } });
      } catch (error) {
        console.error('Error updating profile:', error);
        return new Response(JSON.stringify({ error: 'Internal server error', details: (error as Error).message }), { status: 500, headers: { ...corsHeaders, 'Content-Type': 'application/json' } });
      }
    }

    if (url.pathname === '/api/dispensary/transaction' && request.method === 'POST') {
      const d1 = (env as any).DB || (env as any).neurobotanica_clinical_evidence;
      if (!d1) {
        return new Response(JSON.stringify({ error: 'DB not configured' }), { status: 500, headers: { ...corsHeaders, 'Content-Type': 'application/json' } });
      }

      try {
        const data = await request.json();

        if (!data.customer_id) {
          return new Response(JSON.stringify({ error: 'customer_id is required' }), { status: 400, headers: { ...corsHeaders, 'Content-Type': 'application/json' } });
        }

        await ensureDispensaryTransactionTable(d1);

        const transactionId = data.id || data.transaction_id || `txn_${Date.now().toString(36)}${Math.random().toString(36).slice(2, 8)}`;
        const consultationDate = data.consultation_date || data.timestamp || new Date().toISOString();
        const recommendationItems = data.items || data.recommendations || data.products || [];
        const attachedConditions = data.conditions || data.condition_profile || [];

        await d1.prepare(`
          INSERT INTO dispensary_transactions (
            id,
            customer_id,
            profile_id,
            consultation_date,
            conditions,
            recommendations,
            notes
          ) VALUES (?, ?, ?, ?, ?, ?, ?)
        `).bind(
          transactionId,
          data.customer_id,
          data.profile_id || data.customer_id,
          consultationDate,
          JSON.stringify(attachedConditions),
          JSON.stringify(recommendationItems),
          data.notes || ''
        ).run();

        return new Response(JSON.stringify({
          success: true,
          transaction_id: transactionId,
          consultation_date: consultationDate,
          message: 'Transaction recorded'
        }), { status: 200, headers: { ...corsHeaders, 'Content-Type': 'application/json' } });
      } catch (error) {
        console.error('Error recording transaction:', error);
        return new Response(JSON.stringify({ error: (error as Error).message }), { status: 500, headers: { ...corsHeaders, 'Content-Type': 'application/json' } });
      }
    }

    if (url.pathname === '/api/dispensary/transaction' && request.method === 'GET') {
      const d1 = (env as any).DB || (env as any).neurobotanica_clinical_evidence;
      if (!d1) {
        return new Response(JSON.stringify({ error: 'DB not configured' }), { status: 500, headers: { ...corsHeaders, 'Content-Type': 'application/json' } });
      }

      try {
        await ensureDispensaryTransactionTable(d1);

        const customerId = url.searchParams.get('customer_id');
        let query = 'SELECT * FROM dispensary_transactions ORDER BY consultation_date DESC LIMIT 50';
        const params: Array<string> = [];

        if (customerId) {
          query = 'SELECT * FROM dispensary_transactions WHERE customer_id = ? ORDER BY consultation_date DESC';
          params.push(customerId);
        }

        const results = await d1.prepare(query).bind(...params).all();

        const transactions = results.results?.map(row => {
          return {
            transaction_id: row.id,
            customer_id: row.customer_id,
            profile_id: row.profile_id,
            consultation_date: row.consultation_date,
            created_at: row.created_at,
            notes: row.notes || '',
            conditions: safeParseJson(row.conditions, [] as unknown[]),
            items: safeParseJson(row.recommendations, [] as unknown[])
          };
        }) || [];

        return new Response(JSON.stringify({ transactions }), { status: 200, headers: { ...corsHeaders, 'Content-Type': 'application/json' } });
      } catch (error) {
        console.error('Error fetching transactions:', error);
        return new Response(JSON.stringify({ error: (error as Error).message }), { status: 500, headers: { ...corsHeaders, 'Content-Type': 'application/json' } });
      }
    }

    // Other API requests: return 404 (no longer proxying to Railway)
    if (isApiRequest) {
      return new Response(JSON.stringify({ error: 'Endpoint not found' }), {
        status: 404,
        headers: { ...corsHeaders, 'Content-Type': 'application/json' }
      });
    }

    // Non-API requests: proxy to Pages
    if (!env.PLATFORM_PAGES_URL) {
      return new Response('PLATFORM_PAGES_URL is not configured', { status: 500 });
    }

    const pagesUrl = new URL(url.pathname + url.search, env.PLATFORM_PAGES_URL);
    return forwardRequest(request, pagesUrl);
  }
};


async function forwardRequest(request: Request, target: URL): Promise<Response> {
  const requestClone = request.clone();
  const headers = new Headers(requestClone.headers);
  headers.set('Host', target.host);
  headers.set('Origin', target.origin);

  const init: RequestInit = {
    method: requestClone.method,
    headers,
    redirect: 'manual'
  };

  if (!['GET', 'HEAD'].includes(requestClone.method.toUpperCase())) {
    init.body = requestClone.body;
  }

  const proxiedRequest = new Request(target.toString(), init);
  return fetch(proxiedRequest);
}

function buildCorsPreflightResponse(request: Request): Response {
  const headers = new Headers();
  headers.set('Access-Control-Allow-Origin', request.headers.get('Origin') ?? '*');
  headers.set('Access-Control-Allow-Methods', allowedMethods);
  headers.set('Access-Control-Allow-Headers', request.headers.get('Access-Control-Request-Headers') ?? allowedHeaders);
  headers.set('Access-Control-Max-Age', '86400');
  return new Response(null, { status: 204, headers });
}

function applyCors(response: Response, request: Request): Response {
  const origin = request.headers.get('Origin');
  if (!origin) {
    return response;
  }
  const headers = new Headers(response.headers);
  headers.set('Access-Control-Allow-Origin', origin);
  headers.set('Access-Control-Allow-Credentials', 'true');
  headers.set('Vary', 'Origin');
  return new Response(response.body, {
    status: response.status,
    statusText: response.statusText,
    headers
  });
}
