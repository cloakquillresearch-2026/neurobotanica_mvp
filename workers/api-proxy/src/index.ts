/**
 * NeuroBotanica API Proxy Worker
 * Routes: /api/recommendations (D1-backed clinical evidence)
 */

export interface Env {
  DB?: D1Database;
  API_BASE_URL?: string;
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

// Single export default that handles both the D1-backed recommendations route and
// the existing proxy behavior. This avoids duplicate `export default` declarations.

const apiPrefix = '/api';
const allowedHeaders = 'authorization,content-type,x-client-version';
const allowedMethods = 'GET,HEAD,POST,PUT,PATCH,DELETE,OPTIONS';

export default {
  async fetch(request: Request, env: Env): Promise<Response> {
    const url = new URL(request.url);
    const isApiRequest = url.pathname.startsWith(apiPrefix);

    // CORS helper headers
    const corsHeaders = {
      'Access-Control-Allow-Origin': '*',
      'Access-Control-Allow-Methods': 'GET, POST, OPTIONS',
      'Access-Control-Allow-Headers': 'Content-Type',
    };

    // Handle CORS preflight for general requests
    if (request.method === 'OPTIONS') {
      return new Response(null, { headers: corsHeaders });
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

        // Get condition metadata
        const conditionQuery = await d1.prepare(
          'SELECT * FROM conditions WHERE condition_id = ? OR condition_name LIKE ?'
        ).bind(normalizedCondition, `%${normalizedCondition}%`).first<ConditionData>();

        if (!conditionQuery) {
          return new Response(JSON.stringify({ error: 'Condition not found in evidence database', suggestion: 'Try: ANXIETY, CHRONIC PAIN, PTSD, EPILEPSY, INSOMNIA, GLAUCOMA' }), { status: 404, headers: { ...corsHeaders, 'Content-Type': 'application/json' } });
        }

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
          WHERE condition = ?
          ORDER BY confidence_score DESC
          LIMIT 20
        `).bind(normalizedCondition).all<Study>();

        const studies = studiesResult.results || [];

        const cannabinoidCounts: Record<string, { count: number; totalConfidence: number }> = {};

        studies.forEach(study => {
          try {
            const intervention = JSON.parse(study.intervention || '{}');
            const profile = (intervention.cannabinoid_profile || intervention.cannabinoid || '').toString();
            if (profile.toUpperCase().includes('CBD')) {
              cannabinoidCounts['CBD'] = cannabinoidCounts['CBD'] || { count: 0, totalConfidence: 0 };
              cannabinoidCounts['CBD'].count++;
              cannabinoidCounts['CBD'].totalConfidence += study.confidence_score || 0;
            }
            if (profile.toUpperCase().includes('THC')) {
              cannabinoidCounts['THC'] = cannabinoidCounts['THC'] || { count: 0, totalConfidence: 0 };
              cannabinoidCounts['THC'].count++;
              cannabinoidCounts['THC'].totalConfidence += study.confidence_score || 0;
            }
          } catch (e) {}
        });

        const topCannabinoids = Object.entries(cannabinoidCounts)
          .sort((a, b) => b[1].totalConfidence - a[1].totalConfidence)
          .slice(0, 3)
          .map(([name, data]) => ({
            cannabinoid: name,
            evidence_count: data.count,
            avg_confidence: data.totalConfidence / data.count
          }));

        const recommendation = {
          condition: conditionQuery.condition_name,
          category: conditionQuery.category,
          evidence_summary: {
            total_studies: studies.length,
            study_types: [...new Set(studies.map(s => s.study_type))],
            avg_confidence: studies.length ? studies.reduce((sum, s) => sum + (s.confidence_score || 0), 0) / studies.length : 0
          },
          recommended_cannabinoids: topCannabinoids,
          recommended_ratio: topCannabinoids.length >= 2 ? `${topCannabinoids[0].cannabinoid}:${topCannabinoids[1].cannabinoid} (2:1 to 1:1)` : topCannabinoids[0]?.cannabinoid || 'CBD',
          delivery_methods: ['Tincture', 'Vaporizer', 'Edible'],
          dosing_guidance: severity === 'mild' ? 'Start with 5-10mg, increase gradually' : severity === 'severe' ? 'Consider 20-40mg, consult healthcare provider' : 'Start with 10-20mg, adjust as needed',
          citations: studies.slice(0,5).map(s => ({ study_id: s.study_id, study_type: s.study_type, citation: s.citation, confidence_score: s.confidence_score, key_findings: (() => { try { return JSON.parse(s.key_findings || '[]').slice(0,2) } catch(e) { return [] } })() })),
          confidence_score: studies.length ? studies.reduce((sum, s) => sum + (s.confidence_score || 0), 0) / studies.length : 0,
          disclaimer: 'These recommendations are based on clinical evidence and should not replace medical advice. Consult a healthcare provider before starting any new treatment.'
        };

        return new Response(JSON.stringify(recommendation), { status: 200, headers: { ...corsHeaders, 'Content-Type': 'application/json', 'Cache-Control': 'public, max-age=3600' } });

      } catch (error) {
        console.error('Error generating recommendation:', error);
        return new Response(JSON.stringify({ error: 'Internal server error', details: (error as Error).message }), { status: 500, headers: { ...corsHeaders, 'Content-Type': 'application/json' } });
      }
    }

    // Other API requests: proxy to upstream API_BASE_URL
    if (isApiRequest) {
      if (!env.API_BASE_URL) {
        return new Response('API_BASE_URL is not configured', { status: 500 });
      }

      if (request.method === 'OPTIONS') {
        return buildCorsPreflightResponse(request);
      }

      const targetPath = url.pathname.replace(apiPrefix, '') || '/';
      const upstreamUrl = new URL(targetPath + url.search, env.API_BASE_URL);
      const proxyResponse = await forwardRequest(request, upstreamUrl);

      return applyCors(proxyResponse, request);
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
