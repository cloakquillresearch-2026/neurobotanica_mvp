export async function onRequest(context) {
  // Enable CORS
  const corsHeaders = {
    'Access-Control-Allow-Origin': '*',
    'Access-Control-Allow-Methods': 'GET, POST, OPTIONS',
    'Access-Control-Allow-Headers': 'Content-Type, Authorization',
  };

  // Handle preflight requests
  if (context.request.method === 'OPTIONS') {
    return new Response(null, { headers: corsHeaders });
  }

  try {
    const { request } = context;
    const url = new URL(request.url);
    const compound = url.searchParams.get('compound') || 'THC';

    // Mock ToxPath analysis data
    const toxPathData = {
      compound: compound,
      toxicity_profile: {
        acute_toxicity: {
          ld50_oral: compound === 'THC' ? '482 mg/kg' : '1270 mg/kg',
          ld50_dermal: compound === 'THC' ? '>2000 mg/kg' : '>5000 mg/kg',
          risk_level: compound === 'THC' ? 'Low' : 'Very Low'
        },
        chronic_toxicity: {
          carcinogenicity: 'Non-carcinogenic',
          mutagenicity: 'Non-mutagenic',
          reproductive_toxicity: 'No evidence',
          developmental_toxicity: 'No evidence'
        },
        organ_systems: {
          cardiovascular: 'Minimal risk',
          neurological: compound === 'THC' ? 'Psychoactive effects' : 'Minimal',
          hepatic: 'Low risk',
          renal: 'Low risk',
          respiratory: 'Low risk'
        }
      },
      safety_assessment: {
        therapeutic_index: compound === 'THC' ? 'High' : 'Very High',
        maximum_safe_dose: compound === 'THC' ? '10 mg/kg/day' : '50 mg/kg/day',
        contraindications: compound === 'THC' ? ['Pregnancy', 'Cardiac conditions'] : ['None identified'],
        drug_interactions: compound === 'THC' ? ['CNS depressants', 'MAOIs'] : ['Minimal']
      },
      regulatory_status: {
        fda_schedule: compound === 'THC' ? 'Schedule I' : 'Not scheduled',
        ema_classification: compound === 'THC' ? 'Not approved' : 'Investigational',
        nevada_status: 'Medical cannabis approved'
      },
      timestamp: new Date().toISOString(),
      analysis_version: 'ToxPath-v1.0'
    };

    return new Response(JSON.stringify(toxPathData), {
      headers: {
        ...corsHeaders,
        'Content-Type': 'application/json',
      },
    });

  } catch (error) {
    console.error('ToxPath analysis error:', error);
    return new Response(JSON.stringify({
      error: 'ToxPath analysis failed',
      message: error.message
    }), {
      status: 500,
      headers: {
        ...corsHeaders,
        'Content-Type': 'application/json',
      },
    });
  }
}