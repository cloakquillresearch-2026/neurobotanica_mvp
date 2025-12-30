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
    const condition = url.searchParams.get('condition') || 'anxiety';

    // Mock BioPath analysis data
    const bioPathData = {
      compound: compound,
      target_condition: condition,
      pharmacokinetic_profile: {
        absorption: {
          bioavailability: compound === 'CBD' ? '13-19%' : '6-20%',
          tmax: compound === 'CBD' ? '2-4 hours' : '1-3 hours',
          route: 'Oral administration'
        },
        distribution: {
          volume_distribution: compound === 'CBD' ? '32 L/kg' : '2.5 L/kg',
          protein_binding: compound === 'CBD' ? '94%' : '97%',
          blood_brain_barrier: compound === 'CBD' ? 'Crosses readily' : 'Limited'
        },
        metabolism: {
          primary_enzyme: compound === 'CBD' ? 'CYP2C19, CYP3A4' : 'CYP2C9, CYP3A4',
          metabolites: compound === 'CBD' ? ['7-OH-CBD', '6-OH-CBD'] : ['11-OH-THC', 'THC-COOH'],
          half_life: compound === 'CBD' ? '18-32 hours' : '25-36 hours'
        },
        excretion: {
          primary_route: 'Hepatic metabolism',
          renal_clearance: '<1%',
          fecal_clearance: '70-80%'
        }
      },
      pharmacodynamic_profile: {
        receptor_affinity: {
          cb1_receptor: compound === 'CBD' ? 'Low affinity' : 'High affinity',
          cb2_receptor: compound === 'CBD' ? 'Moderate affinity' : 'Moderate affinity',
          other_targets: compound === 'CBD' ? ['TRPV1', 'GPR55', '5-HT1A'] : ['CB1', 'CB2']
        },
        mechanism_of_action: compound === 'CBD' ? 'Anxiolytic, anti-inflammatory' : 'Analgesic, antiemetic',
        therapeutic_window: compound === 'CBD' ? 'Wide' : 'Moderate'
      },
      clinical_efficacy: {
        condition_specific: {
          [condition]: {
            evidence_level: 'Strong clinical evidence',
            response_rate: compound === 'CBD' ? '65-75%' : '70-80%',
            onset_time: '2-4 weeks',
            maintenance_dose: compound === 'CBD' ? '25-50 mg/day' : '5-15 mg/day'
          }
        },
        biomarkers: ['Cytokine levels', 'Inflammation markers', 'Neurotransmitter balance']
      },
      timestamp: new Date().toISOString(),
      analysis_version: 'BioPath-v1.0'
    };

    return new Response(JSON.stringify(bioPathData), {
      headers: {
        ...corsHeaders,
        'Content-Type': 'application/json',
      },
    });

  } catch (error) {
    console.error('BioPath analysis error:', error);
    return new Response(JSON.stringify({
      error: 'BioPath analysis failed',
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