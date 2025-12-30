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
    const compound = url.searchParams.get('compound') || 'CBD';
    const jurisdiction = url.searchParams.get('jurisdiction') || 'nevada';

    // Mock RegPath analysis data
    const regPathData = {
      compound: compound,
      jurisdiction: jurisdiction,
      regulatory_status: {
        federal_level: {
          dea_schedule: compound === 'CBD' ? 'Not scheduled' : 'Schedule I',
          fda_status: compound === 'CBD' ? 'Investigational new drug' : 'Not approved',
          controlled_substances_act: compound === 'CBD' ? 'Exempt' : 'Controlled'
        },
        state_level: {
          [jurisdiction]: {
            medical_use: 'Approved',
            recreational_use: jurisdiction === 'nevada' ? 'Approved' : 'Not approved',
            cultivation: 'Licensed dispensaries only',
            possession_limits: jurisdiction === 'nevada' ? '1 oz flower, 8g concentrate' : 'Varies by state'
          }
        },
        international_status: {
          who_classification: compound === 'CBD' ? 'No public health risk' : 'High abuse potential',
          un_conventions: compound === 'CBD' ? 'Not controlled' : 'Single Convention 1961',
          eu_status: compound === 'CBD' ? 'Novel food regulation' : 'Not approved'
        }
      },
      compliance_requirements: {
        labeling: {
          required_info: ['THC content', 'CBD content', 'Batch number', 'Expiration date'],
          testing_requirements: ['Potency', 'Contaminants', 'Microbiological'],
          child_resistant: 'Required for all products'
        },
        distribution: {
          licensed_facilities: 'Required',
          tracking_system: 'State seed-to-sale tracking',
          age_verification: '21+ for recreational, varies for medical'
        },
        quality_standards: {
          gmp_compliance: 'Required for manufacturers',
          iso_certification: 'Recommended',
          third_party_testing: 'Mandatory'
        }
      },
      legal_considerations: {
        employment: {
          drug_testing: 'Varies by employer',
          workers_compensation: 'Covered in most states',
          unemployment_benefits: 'Protected in medical states'
        },
        driving: {
          dui_thresholds: '5 ng/mL THC in blood',
          commercial_driving: 'Zero tolerance in most cases',
          medical_exceptions: 'May apply with documentation'
        },
        insurance: {
          health_insurance: 'Covered in some medical plans',
          auto_insurance: 'Rate impacts vary',
          life_insurance: 'Typically not covered'
        }
      },
      business_considerations: {
        taxation: {
          excise_tax: jurisdiction === 'nevada' ? '15% on retail sales' : 'Varies by state',
          sales_tax: 'Standard state sales tax applies',
          cultivation_tax: 'Additional licensing fees'
        },
        licensing: {
          facility_types: ['Cultivation', 'Manufacturing', 'Dispensing', 'Distribution'],
          application_process: 'State regulatory agency review',
          renewal_requirements: 'Annual compliance audits'
        }
      },
      timestamp: new Date().toISOString(),
      analysis_version: 'RegPath-v1.0'
    };

    return new Response(JSON.stringify(regPathData), {
      headers: {
        ...corsHeaders,
        'Content-Type': 'application/json',
      },
    });

  } catch (error) {
    console.error('RegPath analysis error:', error);
    return new Response(JSON.stringify({
      error: 'RegPath analysis failed',
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