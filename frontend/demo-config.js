// frontend/demo-config.js
// Demo environment configuration for Cloudflare Pages demos

export const DEMO_CONFIG = {
  // Demo mode settings
  demoMode: true,
  showDebugInfo: false, // Set to true for development demos

  // Use mock data instead of live API calls
  useMockData: true,

  // Demo API endpoint (can point to a demo backend or use client-side mocks)
  demoApiUrl: 'https://demo-api.neurobotanica.com',

  // Mock response data for demonstrations
  mockResponses: {
    // Sample recommendation response
    recommendations: [
      {
        rank: 1,
        time_of_day: "evening",
        product_name: "Blue Dream (Demo)",
        match_score: 0.87,
        cannabinoid_profile: {
          THC: 18,
          CBD: 1,
          CBG: 1.2,
          CBN: 0.1
        },
        key_terpenes: [
          {
            name: "myrcene",
            percent: 0.8,
            effects: ["sedative", "muscle_relaxant", "anti_inflammatory"]
          },
          {
            name: "limonene",
            percent: 0.3,
            effects: ["mood_elevation", "stress_relief", "energizing"]
          }
        ],
        why_recommended: "High THC for pain relief, myrcene for muscle relaxation, limonene for mood support",
        expected_benefits: ["pain reduction", "muscle relaxation", "better sleep", "mood enhancement"],
        dosage_guidance: {
          starting_dose: "5-10mg THC",
          wait_time_minutes: 30,
          max_dose_per_session: "20mg THC",
          frequency: "2-3 times daily as needed"
        },
        adjuvant_optimization: {
          note: "Consider magnesium glycinate for enhanced relaxation",
          timing: "30 minutes before cannabis",
          expected_synergy: "GABA receptor priming",
          recommended_addition: "300mg magnesium glycinate",
          evidence_level: "Strong clinical evidence"
        }
      },
      {
        rank: 2,
        time_of_day: "morning",
        product_name: "Sour Diesel (Demo)",
        match_score: 0.82,
        cannabinoid_profile: {
          THC: 22,
          CBD: 0.5,
          CBG: 0.8,
          CBC: 0.3
        },
        key_terpenes: [
          {
            name: "limonene",
            percent: 0.9,
            effects: ["energizing", "mood_elevation", "focus"]
          },
          {
            name: "beta_caryophyllene",
            percent: 0.4,
            effects: ["anti_inflammatory", "pain_relief", "CB2_activation"]
          }
        ],
        why_recommended: "High limonene content for morning energy and focus",
        expected_benefits: ["energy boost", "mood elevation", "mental clarity"],
        dosage_guidance: {
          starting_dose: "2.5-5mg THC",
          wait_time_minutes: 20,
          max_dose_per_session: "15mg THC",
          frequency: "1-2 times daily as needed"
        }
      }
    ],

    // Sample condition profiles
    conditions: [
      {
        id: "chronic_pain",
        name: "Chronic Pain",
        description: "Long-term pain management",
        cannabinoid_recommendations: {
          primary: "THC",
          secondary: "CBD",
          ratio: "2:1 THC:CBD"
        },
        terpene_focus: ["myrcene", "beta_caryophyllene"],
        clinical_studies: 89
      },
      {
        id: "anxiety",
        name: "Anxiety",
        description: "Generalized anxiety and stress",
        cannabinoid_recommendations: {
          primary: "CBD",
          secondary: "CBG",
          ratio: "1:1 CBD:THC or CBD dominant"
        },
        terpene_focus: ["linalool", "limonene"],
        clinical_studies: 67
      },
      {
        id: "insomnia",
        name: "Insomnia",
        description: "Sleep difficulties and disorders",
        cannabinoid_recommendations: {
          primary: "CBN",
          secondary: "THC",
          ratio: "Low dose THC with CBN"
        },
        terpene_focus: ["myrcene", "terpinolene"],
        clinical_studies: 43
      }
    ],

    // Sample adjuvant recommendations
    adjuvants: [
      {
        name: "Magnesium Glycinate",
        dosage: "300-400mg",
        timing: "30 minutes before cannabis",
        mechanism: "GABA receptor modulation",
        evidence_level: "Strong",
        pmids: ["PMC123456", "PMC789012"]
      },
      {
        name: "L-Theanine",
        dosage: "200mg",
        timing: "With cannabis",
        mechanism: "Glutamate receptor modulation",
        evidence_level: "Moderate",
        pmids: ["PMC345678"]
      },
      {
        name: "Curcumin",
        dosage: "500mg",
        timing: "With food",
        mechanism: "Inflammation reduction",
        evidence_level: "Strong",
        pmids: ["PMC901234", "PMC567890"]
      }
    ]
  },

  // Demo UI customizations
  ui: {
    showDemoBadge: true,
    demoBadgeText: "DEMO ENVIRONMENT",
    demoBadgeColor: "bg-blue-500",

    // Show additional info in demo mode
    showTradeSecretInfo: true,
    showClinicalEvidenceCount: true,
    showProcessingTime: true,

    // Mock loading states
    mockApiDelay: 800, // milliseconds
  },

  // Demo analytics (optional)
  analytics: {
    trackDemoUsage: true,
    demoEvents: [
      'demo_recommendation_generated',
      'demo_condition_selected',
      'demo_adjuvant_viewed'
    ]
  }
};

export default DEMO_CONFIG;