export async function onRequest(context) {
  const { request } = context;

  // Handle CORS
  if (request.method === 'OPTIONS') {
    return new Response(null, {
      status: 200,
      headers: {
        'Access-Control-Allow-Origin': '*',
        'Access-Control-Allow-Methods': 'GET, POST, PUT, DELETE, OPTIONS',
        'Access-Control-Allow-Headers': 'Content-Type, Authorization, X-API-Key',
      },
    });
  }

  if (request.method !== 'POST') {
    return new Response(JSON.stringify({ error: 'Method not allowed' }), {
      status: 405,
      headers: { 'Content-Type': 'application/json' },
    });
  }

  try {
    const requestData = await request.json();

    // Mock recommendation logic based on condition
    const condition = requestData.condition || 'chronic_pain';
    const experienceLevel = requestData.experience_level || 'intermediate';

    // Mock recommendations based on condition
    const recommendations = generateMockRecommendations(condition, experienceLevel);

    const response = {
      recommendation_id: `rec_${Date.now()}`,
      recommendations: recommendations,
      products_to_avoid: [
        {
          product_name: "High THC Strain",
          reason: "May be too potent for current tolerance level"
        }
      ],
      education_notes: [
        "Start low and go slow with dosage",
        "Monitor effects for 2 hours after consumption",
        "Consider functional requirements for daily activities"
      ],
      generated_at: new Date().toISOString(),
      disclaimer: "This is educational information, not medical advice. Consult healthcare provider."
    };

    return new Response(JSON.stringify(response), {
      status: 200,
      headers: {
        'Content-Type': 'application/json',
        'Access-Control-Allow-Origin': '*',
        'Access-Control-Allow-Headers': 'Content-Type, Authorization, X-API-Key',
      },
    });

  } catch (error) {
    return new Response(JSON.stringify({
      error: 'Invalid request format',
      message: error.message
    }), {
      status: 400,
      headers: { 'Content-Type': 'application/json' },
    });
  }
}

function generateMockRecommendations(condition, experienceLevel) {
  const conditionProfiles = {
    chronic_pain: {
      primary: { thc: 0.8, cbd: 0.7, cbg: 0.6 },
      strains: ["Blue Dream", "OG Kush", "Northern Lights"],
      terpenes: ["myrcene", "limonene", "beta-caryophyllene"]
    },
    anxiety: {
      primary: { thc: 0.3, cbd: 0.85, cbg: 0.6 },
      strains: ["ACDC", "Harlequin", "Cannatonic"],
      terpenes: ["linalool", "limonene", "beta-caryophyllene"]
    },
    insomnia: {
      primary: { thc: 0.7, cbd: 0.5, cbn: 0.85 },
      strains: ["Granddaddy Purple", "Afghani", "Northern Lights"],
      terpenes: ["myrcene", "terpinolene", "beta-caryophyllene"]
    },
    inflammation: {
      primary: { thc: 0.6, cbd: 0.8, cbg: 0.75 },
      strains: ["Charlotte's Web", "Ringo's Gift", "Critical Mass"],
      terpenes: ["beta-caryophyllene", "limonene", "pinene"]
    }
  };

  const profile = conditionProfiles[condition] || conditionProfiles.chronic_pain;

  // Adjust for experience level
  const dosageMultiplier = {
    beginner: 0.5,
    intermediate: 0.75,
    advanced: 1.0
  }[experienceLevel] || 0.75;

  return [
    {
      rank: 1,
      time_of_day: "evening",
      product_name: profile.strains[0],
      match_score: 0.87,
      cannabinoid_profile: {
        THC: Math.round(profile.primary.thc * 20 * dosageMultiplier * 10) / 10,
        CBD: Math.round(profile.primary.cbd * 20 * dosageMultiplier * 10) / 10,
        CBG: Math.round((profile.primary.cbg || 0) * 20 * dosageMultiplier * 10) / 10
      },
      key_terpenes: profile.terpenes.slice(0, 2).map(terpene => ({
        name: terpene,
        percent: Math.random() * 0.5 + 0.3,
        effects: getTerpeneEffects(terpene)
      })),
      why_recommended: `High match for ${condition.replace('_', ' ')} based on clinical evidence`,
      expected_benefits: getExpectedBenefits(condition),
      dosage_guidance: {
        starting_dose: `${Math.round(2.5 * dosageMultiplier * 10) / 10}-${Math.round(5 * dosageMultiplier * 10) / 10}mg THC`,
        wait_time_minutes: 30,
        max_dose_per_session: `${Math.round(10 * dosageMultiplier * 10) / 10}mg THC`,
        frequency: "2-3 times daily as needed"
      },
      adjuvant_optimization: {
        note: "Consider magnesium glycinate for enhanced relaxation",
        timing: "30 minutes before cannabis",
        expected_synergy: "GABA receptor priming",
        recommended_addition: "300mg magnesium glycinate"
      }
    }
  ];
}

function getTerpeneEffects(terpene) {
  const effects = {
    myrcene: ["sedative", "muscle_relaxant"],
    limonene: ["mood_elevation", "stress_relief"],
    linalool: ["anxiolytic", "sedative"],
    "beta-caryophyllene": ["anti_inflammatory", "pain_relief"],
    terpinolene: ["sedative", "antioxidant"],
    pinene: ["alertness", "anti_inflammatory"]
  };
  return effects[terpene] || ["relaxing"];
}

function getExpectedBenefits(condition) {
  const benefits = {
    chronic_pain: ["pain reduction", "muscle relaxation", "better sleep"],
    anxiety: ["reduced anxiety", "calmness", "stress relief"],
    insomnia: ["improved sleep", "relaxation", "reduced sleep latency"],
    inflammation: ["reduced inflammation", "pain relief", "immune support"]
  };
  return benefits[condition] || ["therapeutic benefits"];
}