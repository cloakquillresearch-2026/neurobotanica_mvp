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
    const genetic_markers = url.searchParams.get('markers')?.split(',') || ['CYP2C19', 'COMT'];

    // Mock GenomePath analysis data
    const genomePathData = {
      compound: compound,
      genetic_analysis: {
        pharmacogenomic_markers: genetic_markers.map(marker => ({
          gene: marker,
          variant: marker === 'CYP2C19' ? '*2/*2' : 'Val/Val',
          phenotype: marker === 'CYP2C19' ? 'Poor metabolizer' : 'Normal metabolizer',
          impact_on_metabolism: marker === 'CYP2C19' ? 'Reduced clearance' : 'Normal clearance',
          dosing_recommendation: marker === 'CYP2C19' ? 'Reduce dose by 50%' : 'Standard dose'
        })),
        metabolic_pathways: {
          phase1_metabolism: {
            enzymes: ['CYP2C19', 'CYP3A4', 'CYP2C9'],
            genetic_influence: 'High',
            variability_factor: '3-5x'
          },
          phase2_metabolism: {
            enzymes: ['UGT1A9', 'UGT2B7'],
            genetic_influence: 'Moderate',
            variability_factor: '2-3x'
          }
        },
        transporter_proteins: {
          abc_transporters: ['ABCB1', 'ABCC1'],
          slc_transporters: ['SLC21A6'],
          genetic_impact: 'Moderate influence on bioavailability'
        }
      },
      personalized_recommendations: {
        optimal_dose_range: compound === 'CBD' ? '15-30 mg/day' : '2.5-10 mg/day',
        titration_schedule: 'Start low, go slow - increase by 25% weekly',
        monitoring_parameters: ['Drug levels', 'Clinical response', 'Adverse effects'],
        alternative_compounds: compound === 'CBD' ? ['CBG', 'CBN'] : ['THCA', 'CBDV']
      },
      risk_assessment: {
        drug_interactions: {
          high_risk: ['Strong CYP3A4 inhibitors'],
          moderate_risk: ['Moderate CYP2C19 inhibitors'],
          low_risk: ['Most other medications']
        },
        adverse_effect_risk: 'Low to moderate',
        therapeutic_efficacy: 'High likelihood of response'
      },
      research_evidence: {
        clinical_studies: 45,
        genetic_associations: 23,
        confidence_score: '85%',
        last_updated: '2024-12-01'
      },
      timestamp: new Date().toISOString(),
      analysis_version: 'GenomePath-v1.0'
    };

    return new Response(JSON.stringify(genomePathData), {
      headers: {
        ...corsHeaders,
        'Content-Type': 'application/json',
      },
    });

  } catch (error) {
    console.error('GenomePath analysis error:', error);
    return new Response(JSON.stringify({
      error: 'GenomePath analysis failed',
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