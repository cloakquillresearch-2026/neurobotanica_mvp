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

    // Mock MetaPath analysis data
    const metaPathData = {
      compound: compound,
      target_condition: condition,
      meta_analysis: {
        systematic_reviews: {
          total_reviews: 23,
          high_quality_reviews: 18,
          pooled_effect_size: 0.78,
          heterogeneity_i2: '35%',
          confidence_interval: '0.65-0.91'
        },
        clinical_trials: {
          randomized_controlled: 45,
          total_participants: 3247,
          average_jadad_score: 4.2,
          publication_bias: 'Low risk',
          quality_assessment: 'High quality evidence'
        },
        observational_studies: {
          cohort_studies: 67,
          case_control_studies: 34,
          cross_sectional_studies: 89,
          total_participants: 15643,
          effect_consistency: 'High'
        }
      },
      evidence_synthesis: {
        strength_of_evidence: {
          overall_grade: 'A (High)',
          consistency: 'Consistent findings',
          directness: 'Direct evidence',
          precision: 'Precise estimates',
          reporting_bias: 'Low risk'
        },
        clinical_recommendations: {
          primary_recommendation: `Strong recommendation for ${compound} in ${condition}`,
          strength: 'Grade 1A',
          conditions_apply: ['Adult patients', 'Medical supervision', 'Regular monitoring']
        },
        comparative_effectiveness: {
          vs_placebo: 'SMD 0.82 (95% CI: 0.65-0.99)',
          vs_standard_care: 'SMD 0.45 (95% CI: 0.28-0.62)',
          vs_other_cannabinoids: 'Similar efficacy profile',
          network_meta_analysis: 'Ranked 2nd of 8 interventions'
        }
      },
      safety_profile: {
        adverse_events: {
          common: ['Drowsiness', 'Dry mouth', 'Dizziness'],
          serious: ['Psychosis (rare)', 'Cardiovascular events (rare)'],
          withdrawal: 'Minimal withdrawal symptoms',
          dependency: 'Low dependency potential'
        },
        risk_benefit_ratio: {
          overall: 'Favorable',
          number_needed_to_treat: 4,
          number_needed_to_harm: 28,
          quality_adjusted_life_years: '+0.15 QALYs'
        }
      },
      implementation_guidelines: {
        dosing_strategy: {
          initiation: 'Start low (10-25 mg/day)',
          titration: 'Increase by 25% weekly',
          maintenance: '25-100 mg/day based on response',
          maximum_dose: '200 mg/day'
        },
        monitoring_protocol: {
          clinical_assessment: 'Weekly for first month',
          therapeutic_drug_monitoring: 'Optional',
          adverse_event_monitoring: 'Ongoing',
          efficacy_assessment: 'Patient reported outcomes'
        },
        patient_selection: {
          inclusion_criteria: ['Diagnosis confirmed', 'Failed conventional treatments', 'No contraindications'],
          exclusion_criteria: ['Pregnancy', 'Severe psychiatric illness', 'Cardiac disease'],
          special_populations: ['Elderly: reduce dose', 'Children: limited data', 'Hepatic impairment: caution']
        }
      },
      research_gaps: {
        identified_gaps: [
          'Long-term safety data (>2 years)',
          'Pediatric population studies',
          'Drug interaction studies',
          'Optimal dosing regimens'
        ],
        ongoing_research: 12,
        priority_areas: ['Comparative effectiveness', 'Cost-effectiveness', 'Implementation science']
      },
      timestamp: new Date().toISOString(),
      analysis_version: 'MetaPath-v1.0'
    };

    return new Response(JSON.stringify(metaPathData), {
      headers: {
        ...corsHeaders,
        'Content-Type': 'application/json',
      },
    });

  } catch (error) {
    console.error('MetaPath analysis error:', error);
    return new Response(JSON.stringify({
      error: 'MetaPath analysis failed',
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