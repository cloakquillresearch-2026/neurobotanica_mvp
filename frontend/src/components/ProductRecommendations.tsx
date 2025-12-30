import { useState, useEffect } from 'react'
import { dispensaryAPI, adjuvantAPI } from '@/utils/api'

interface ProductRecommendationsProps {
  customer: any
  recommendations: any[]
  onRecommendationsUpdate: (recommendations: any[]) => void
}

interface Recommendation {
  product_id: string
  product_name: string
  product_type: string
  thc_percent: number
  cbd_percent: number
  confidence_score?: number
  therapeutic_rationale?: string
  terpene_profile?: string[]
  clinical_studies?: number
  adjuvants?: AdjuvantInfo[]
}

interface AdjuvantInfo {
  name: string
  mechanism: string
  dosage: string
  timing: string
  evidence: string
}

const TERPENE_INFO: Record<string, { emoji: string; color: string; effect: string; description: string }> = {
  'Myrcene': { emoji: '🥭', color: 'terpene-myrcene', effect: 'Relaxation', description: 'Enhances THC absorption, promotes sedation' },
  'Limonene': { emoji: '🍋', color: 'terpene-limonene', effect: 'Mood', description: 'Elevates mood, relieves stress and anxiety' },
  'Pinene': { emoji: '🌲', color: 'terpene-pinene', effect: 'Focus', description: 'Improves alertness, counteracts THC memory effects' },
  'Linalool': { emoji: '💜', color: 'terpene-linalool', effect: 'Calm', description: 'Reduces anxiety, promotes relaxation' },
  'Caryophyllene': { emoji: '🌶️', color: 'terpene-caryophyllene', effect: 'Pain', description: 'Anti-inflammatory, activates CB2 receptors' },
  'Beta-Caryophyllene': { emoji: '🌶️', color: 'terpene-caryophyllene', effect: 'Pain', description: 'Anti-inflammatory, activates CB2 receptors' },
}

const ADJUVANT_DETAILS: Record<string, AdjuvantInfo> = {
  'Magnesium Glycinate': {
    name: 'Magnesium Glycinate',
    mechanism: 'Enhances GABA receptor activity, relaxes muscles, reduces nerve excitability',
    dosage: '200-400mg daily',
    timing: 'Take 30-60 min before cannabis for enhanced relaxation',
    evidence: 'Supported by 47 clinical studies on pain and sleep'
  },
  'Curcumin': {
    name: 'Curcumin (Turmeric Extract)',
    mechanism: 'Potent anti-inflammatory, inhibits COX-2, enhances cannabinoid receptor sensitivity',
    dosage: '500-1000mg with black pepper (piperine)',
    timing: 'Take with meals, 1-2 hours before cannabis',
    evidence: 'Over 200 studies support anti-inflammatory synergy'
  },
  'L-Theanine': {
    name: 'L-Theanine',
    mechanism: 'Promotes alpha brain waves, reduces anxiety without sedation, balances THC stimulation',
    dosage: '100-200mg',
    timing: 'Take 15-30 min before cannabis',
    evidence: '34 studies show anxiety reduction and focus enhancement'
  },
  'Omega-3': {
    name: 'Omega-3 Fatty Acids',
    mechanism: 'Supports endocannabinoid system function, reduces neuroinflammation',
    dosage: '1000-2000mg EPA/DHA daily',
    timing: 'Take consistently with meals',
    evidence: 'Enhances CB1 receptor density and sensitivity'
  },
  'Black Pepper Extract': {
    name: 'Black Pepper Extract (Piperine)',
    mechanism: 'Activates CB2 receptors, enhances bioavailability of cannabinoids by 2000%',
    dosage: '5-20mg',
    timing: 'Take with cannabis products',
    evidence: 'Clinically proven to enhance absorption of multiple compounds'
  },
  'Vitamin D3': {
    name: 'Vitamin D3',
    mechanism: 'Regulates endocannabinoid tone, supports immune function and mood',
    dosage: '2000-5000 IU daily',
    timing: 'Take with fatty meal for absorption',
    evidence: 'Deficiency linked to reduced cannabinoid receptor expression'
  },
  'NAC': {
    name: 'N-Acetyl Cysteine (NAC)',
    mechanism: 'Glutathione precursor, modulates glutamate, reduces oxidative stress and may reduce cannabis tolerance',
    dosage: '600-1800mg daily',
    timing: 'Take 1 hour before cannabis',
    evidence: 'Clinical studies show benefits for PTSD, depression, and neuroprotection'
  },
  'Taurine': {
    name: 'Taurine',
    mechanism: 'GABA-A receptor modulation, glycine receptor agonism, neuroprotective',
    dosage: '500-3000mg daily',
    timing: 'Take 30 min before cannabis',
    evidence: 'Supports anxiety reduction and seizure threshold'
  },
  'Ashwagandha': {
    name: 'Ashwagandha (KSM-66)',
    mechanism: 'Reduces cortisol, GABA mimetic effects, adaptogenic stress relief',
    dosage: '300-600mg daily',
    timing: 'Take with meals, morning or evening',
    evidence: 'Multiple RCTs show 50% anxiety reduction and improved sleep'
  },
  'CoQ10': {
    name: 'Coenzyme Q10 (Ubiquinol)',
    mechanism: 'Mitochondrial energy production, powerful antioxidant, supports cellular health',
    dosage: '100-400mg daily (ubiquinol form)',
    timing: 'Take with fatty meal',
    evidence: 'Strong evidence for migraine prevention and neuroprotection'
  },
  'Alpha-Lipoic Acid': {
    name: 'Alpha-Lipoic Acid (ALA)',
    mechanism: 'Universal antioxidant, regenerates other antioxidants, supports nerve health',
    dosage: '300-1200mg daily (R-ALA preferred)',
    timing: 'Take 30 min before cannabis on empty stomach',
    evidence: 'Clinically proven for neuropathic pain and neuroprotection'
  },
  'Phosphatidylserine': {
    name: 'Phosphatidylserine (PS)',
    mechanism: 'Blunts cortisol response, supports cell membranes, enhances neurotransmitter release',
    dosage: '100-300mg daily',
    timing: 'Take with meals',
    evidence: 'Studies show reduced stress response and improved cognitive function'
  },
  'Berberine': {
    name: 'Berberine',
    mechanism: 'Activates AMPK (master metabolic switch), improves insulin sensitivity, modulates GLP-1 pathway similar to Ozempic',
    dosage: '500-1500mg daily in divided doses',
    timing: 'Take 30 min before meals',
    evidence: 'Meta-analyses show comparable efficacy to metformin for blood sugar control'
  },
  'Inositol': {
    name: 'Myo-Inositol',
    mechanism: 'Enhances insulin signaling, supports healthy ovarian function, reduces cravings',
    dosage: '2000-4000mg daily',
    timing: 'Take with or without food',
    evidence: 'Tier 1 evidence for PCOS and metabolic syndrome'
  },
  'Chromium': {
    name: 'Chromium Picolinate',
    mechanism: 'Enhances insulin receptor sensitivity, improves glucose transport, reduces carb cravings',
    dosage: '200-1000mcg daily',
    timing: 'Take with meals',
    evidence: 'Clinical studies show reduced food cravings and improved glucose metabolism'
  },
  'EGCG': {
    name: 'Green Tea Extract (EGCG)',
    mechanism: 'Increases thermogenesis, enhances fat oxidation, activates AMPK for metabolism boost',
    dosage: '250-500mg daily',
    timing: 'Take 30 min before exercise or morning',
    evidence: 'Multiple RCTs show increased fat burning and metabolic rate'
  },
  // MARINE KINGDOM ADJUVANTS
  'Krill Oil': {
    name: 'Krill Oil (Omega-3 + Astaxanthin)',
    mechanism: 'Phospholipid-bound omega-3 + astaxanthin antioxidant + enhanced absorption',
    dosage: '500-1000mg daily',
    timing: 'Take with meals',
    evidence: 'Superior bioavailability compared to fish oil, enhanced joint and cognitive support'
  },
  'Astaxanthin': {
    name: 'Astaxanthin',
    mechanism: 'Powerful antioxidant + mitochondrial protection + anti-inflammatory + crosses blood-brain barrier',
    dosage: '4-12mg daily',
    timing: 'Take with fatty meal',
    evidence: 'Most powerful carotenoid antioxidant, supports eye, skin, and brain health'
  },
  'Marine Omega-3': {
    name: 'Marine Omega-3 Concentrate',
    mechanism: 'High-EPA/DHA ratio + phospholipid form + enhanced endocannabinoid system support',
    dosage: '750-1500mg daily',
    timing: 'Take with meals',
    evidence: 'Clinical studies show superior cardiovascular and anti-inflammatory benefits'
  },
  // FUNGAL KINGDOM ADJUVANTS
  'Reishi': {
    name: 'Reishi Mushroom (Ganoderma lucidum)',
    mechanism: 'Cortisol modulation + immune enhancement + triterpenes + polysaccharides',
    dosage: '500-1000mg daily',
    timing: 'Take with meals',
    evidence: 'Traditional medicine validated by modern studies for stress and immune support'
  },
  'Lion Mane': {
    name: 'Lion\'s Mane Mushroom (Hericium erinaceus)',
    mechanism: 'NGF stimulation + myelin repair + cognitive enhancement + erinacines',
    dosage: '500-1000mg daily',
    timing: 'Take with meals',
    evidence: 'Clinical studies show nerve regeneration and cognitive enhancement'
  },
  'Cordyceps': {
    name: 'Cordyceps Mushroom (Cordyceps sinensis)',
    mechanism: 'AMPK activation + mitochondrial enhancement + adenosine + cordycepin',
    dosage: '500-1000mg daily',
    timing: 'Take 30 min before exercise',
    evidence: 'Olympic athletes use for endurance, clinical studies support energy production'
  },
  'Turkey Tail': {
    name: 'Turkey Tail Mushroom (Trametes versicolor)',
    mechanism: 'PSP/PSG polysaccharides + immune modulation + beta-glucans',
    dosage: '500-1000mg daily',
    timing: 'Take with meals',
    evidence: 'Strong evidence for immune support and cancer adjunct therapy'
  },
  // PLANT KINGDOM ADJUVANTS
  'Bacopa': {
    name: 'Bacopa Monnieri',
    mechanism: 'BDNF enhancement + synaptic plasticity + bacosides + antioxidant',
    dosage: '150-300mg daily',
    timing: 'Take with meals',
    evidence: 'Ayurvedic herb with clinical validation for cognitive enhancement'
  },
  'Rhodiola': {
    name: 'Rhodiola Rosea',
    mechanism: 'Cortisol modulation + serotonin enhancement + salidrosides + rosavins',
    dosage: '100-200mg daily',
    timing: 'Take morning',
    evidence: 'Adaptogenic herb with clinical studies for fatigue and mental performance'
  },
  'Holy Basil': {
    name: 'Holy Basil (Tulsi)',
    mechanism: 'Cortisol reduction + GABA enhancement + eugenol + ursolic acid',
    dosage: '150-300mg daily',
    timing: 'Take with meals',
    evidence: 'Sacred herb with modern validation for stress and metabolic support'
  },
  'Gotu Kola': {
    name: 'Gotu Kola (Centella asiatica)',
    mechanism: 'Triterpenes + asiaticoside + cognitive enhancement + collagen synthesis',
    dosage: '250-500mg daily',
    timing: 'Take with meals',
    evidence: 'Traditional wound healer with clinical studies for circulation and cognition'
  },
  'Milk Thistle': {
    name: 'Milk Thistle (Silymarin)',
    mechanism: 'Hepatoprotection + glutathione enhancement + antioxidant + silibinin',
    dosage: '150-300mg daily',
    timing: 'Take with meals',
    evidence: 'Gold standard liver support with extensive clinical validation'
  },
  'Boswellia': {
    name: 'Boswellia (Frankincense)',
    mechanism: '5-LOX inhibition + AKBA + boswellic acids + COX-2 modulation',
    dosage: '250-500mg daily',
    timing: 'Take with meals',
    evidence: 'Clinical studies show joint health and anti-inflammatory benefits'
  },
  'Quercetin': {
    name: 'Quercetin',
    mechanism: 'Flavonoid + antioxidant + mast cell stabilization + CYP modulation',
    dosage: '250-500mg daily',
    timing: 'Take with meals',
    evidence: 'Bioflavonoid with clinical studies for allergies and cardiovascular health'
  },
  'Resveratrol': {
    name: 'Resveratrol',
    mechanism: 'SIRT1 activation + AMPK enhancement + antioxidant + phytoalexin',
    dosage: '50-100mg daily',
    timing: 'Take with meals',
    evidence: 'Red wine compound with clinical studies for longevity and metabolic health'
  },
  'Ginger': {
    name: 'Ginger',
    mechanism: 'Anti-nausea + anti-inflammatory + gingerols + shogaols',
    dosage: '500-1000mg daily',
    timing: 'Take 30 min before meals',
    evidence: 'Clinical studies show effectiveness for nausea and digestive support'
  },
  'Turmeric': {
    name: 'Turmeric',
    mechanism: 'Curcumin + anti-inflammatory + antioxidant + COX-2 inhibition',
    dosage: '500-1000mg daily with black pepper',
    timing: 'Take with meals',
    evidence: 'Extensive studies support anti-inflammatory and joint health benefits'
  },
}

export function ProductRecommendations({
  customer,
  recommendations,
  onRecommendationsUpdate,
}: ProductRecommendationsProps) {
  const [loading, setLoading] = useState(false)
  const [expandedProduct, setExpandedProduct] = useState<string | null>(null)
  const [expandedAdjuvant, setExpandedAdjuvant] = useState<string | null>(null)

  useEffect(() => {
    if (customer && !customer.isNew) {
      fetchRecommendations()
    }
  }, [customer])

  const fetchRecommendations = async () => {
    setLoading(true)
    try {
      const conditions = customer.conditions || []
      
      // Build recommendations based on conditions
      const mockRecommendations: Recommendation[] = []
      
      if (conditions.includes('chronic_pain') || conditions.includes('inflammation')) {
        mockRecommendations.push({
          product_id: 'PAIN_001',
          product_name: 'High THC Indica Formulation',
          product_type: 'Recommended Profile',
          thc_percent: 18,
          cbd_percent: 2,
          confidence_score: 0.92,
          therapeutic_rationale: 'Indica-dominant strains with high myrcene content provide effective pain relief through CB1 receptor activation and anti-inflammatory terpene synergy. The sedative properties also help with pain-related sleep disturbances.',
          terpene_profile: ['Myrcene', 'Beta-Caryophyllene', 'Linalool'],
          clinical_studies: 47,
          adjuvants: [
            ADJUVANT_DETAILS['Magnesium Glycinate'],
            ADJUVANT_DETAILS['Curcumin'],
            ADJUVANT_DETAILS['Alpha-Lipoic Acid'],
            ADJUVANT_DETAILS['Krill Oil'], // Marine kingdom
            ADJUVANT_DETAILS['Boswellia'], // Plant kingdom
          ],
        })
      }
      
      if (conditions.includes('anxiety') || conditions.includes('ptsd')) {
        mockRecommendations.push({
          product_id: 'ANXIETY_001',
          product_name: 'CBD-Dominant Calming Blend',
          product_type: 'Recommended Profile',
          thc_percent: 2,
          cbd_percent: 20,
          confidence_score: 0.89,
          therapeutic_rationale: 'CBD-dominant formulations reduce anxiety by modulating serotonin receptors (5-HT1A) without psychoactive effects. Linalool and limonene terpenes provide additional anxiolytic support.',
          terpene_profile: ['Linalool', 'Limonene', 'Pinene'],
          clinical_studies: 63,
          adjuvants: [
            ADJUVANT_DETAILS['L-Theanine'],
            ADJUVANT_DETAILS['Ashwagandha'],
            ADJUVANT_DETAILS['NAC'],
            ADJUVANT_DETAILS['Reishi'], // Fungal kingdom
            ADJUVANT_DETAILS['Holy Basil'], // Plant kingdom
          ],
        })
      }
      
      if (conditions.includes('insomnia')) {
        mockRecommendations.push({
          product_id: 'SLEEP_001',
          product_name: 'Nighttime Rest Formula',
          product_type: 'Recommended Profile',
          thc_percent: 15,
          cbd_percent: 5,
          confidence_score: 0.87,
          therapeutic_rationale: 'CBN-rich indica strains with high myrcene promote deep, restorative sleep. The sedative terpene profile helps reduce sleep latency and improves sleep quality without morning grogginess.',
          terpene_profile: ['Myrcene', 'Linalool', 'Caryophyllene'],
          clinical_studies: 31,
          adjuvants: [
            ADJUVANT_DETAILS['Magnesium Glycinate'],
            ADJUVANT_DETAILS['L-Theanine'],
            ADJUVANT_DETAILS['Reishi'], // Fungal kingdom
            ADJUVANT_DETAILS['Holy Basil'], // Plant kingdom
          ],
        })
      }

      if (conditions.includes('depression')) {
        mockRecommendations.push({
          product_id: 'MOOD_001',
          product_name: 'Mood Elevation Sativa',
          product_type: 'Recommended Profile',
          thc_percent: 12,
          cbd_percent: 8,
          confidence_score: 0.85,
          therapeutic_rationale: 'Balanced sativa strains with limonene provide mood elevation and motivation. The moderate THC with CBD prevents anxiety while promoting positive mental state.',
          terpene_profile: ['Limonene', 'Pinene', 'Beta-Caryophyllene'],
          clinical_studies: 28,
          adjuvants: [
            ADJUVANT_DETAILS['Omega-3'],
            ADJUVANT_DETAILS['NAC'],
            ADJUVANT_DETAILS['Phosphatidylserine'],
            ADJUVANT_DETAILS['Lion Mane'], // Fungal kingdom
            ADJUVANT_DETAILS['Rhodiola'], // Plant kingdom
          ],
        })
      }

      if (conditions.includes('nausea')) {
        mockRecommendations.push({
          product_id: 'NAUSEA_001',
          product_name: 'Digestive Support Blend',
          product_type: 'Recommended Profile',
          thc_percent: 15,
          cbd_percent: 3,
          confidence_score: 0.91,
          therapeutic_rationale: 'THC activates CB1 receptors in the gut-brain axis, reducing nausea and stimulating appetite. Limonene provides additional anti-nausea effects.',
          terpene_profile: ['Limonene', 'Myrcene', 'Caryophyllene'],
          clinical_studies: 42,
          adjuvants: [
            ADJUVANT_DETAILS['Black Pepper Extract'],
            ADJUVANT_DETAILS['Ginger'], // Plant kingdom
            ADJUVANT_DETAILS['Reishi'], // Fungal kingdom
            ADJUVANT_DETAILS['Turmeric'], // Plant kingdom
          ],
        })
      }

      if (conditions.includes('weight_management')) {
        mockRecommendations.push({
          product_id: 'METABOLIC_001',
          product_name: 'THCV Metabolic Formula',
          product_type: 'Recommended Profile',
          thc_percent: 5,
          cbd_percent: 10,
          confidence_score: 0.88,
          therapeutic_rationale: 'THCV (tetrahydrocannabivarin) is known as "diet weed" - it suppresses appetite unlike THC, improves insulin sensitivity, and provides energizing effects. Combined with CBD for metabolic support without munchies.',
          terpene_profile: ['Limonene', 'Pinene', 'Beta-Caryophyllene'],
          clinical_studies: 23,
          adjuvants: [
            ADJUVANT_DETAILS['Berberine'],
            ADJUVANT_DETAILS['EGCG'],
            ADJUVANT_DETAILS['Inositol'],
            ADJUVANT_DETAILS['Cordyceps'], // Fungal kingdom
            ADJUVANT_DETAILS['EGCG'], // Plant kingdom
          ],
        })
      }

      // Default if no conditions match
      if (mockRecommendations.length === 0) {
        mockRecommendations.push({
          product_id: 'BALANCED_001',
          product_name: 'Balanced Wellness Profile',
          product_type: 'Recommended Profile',
          thc_percent: 10,
          cbd_percent: 10,
          confidence_score: 0.80,
          therapeutic_rationale: 'A balanced 1:1 THC:CBD ratio provides therapeutic benefits with minimal psychoactive intensity. Ideal for new users or those seeking general wellness support.',
          terpene_profile: ['Limonene', 'Pinene', 'Myrcene'],
          clinical_studies: 55,
          adjuvants: [
            ADJUVANT_DETAILS['L-Theanine'],
            ADJUVANT_DETAILS['Omega-3'],
            ADJUVANT_DETAILS['Reishi'], // Fungal kingdom
            ADJUVANT_DETAILS['Turmeric'], // Plant kingdom
          ],
        })
      }

      onRecommendationsUpdate(mockRecommendations)
    } catch (error) {
      console.error('Failed to fetch recommendations:', error)
    } finally {
      setLoading(false)
    }
  }

  const getMatchBadgeClass = (score: number) => {
    if (score >= 0.9) return 'match-badge-high'
    if (score >= 0.75) return 'match-badge-medium'
    return 'match-badge-low'
  }

  if (loading) {
    return (
      <div className="text-center py-12">
        <div className="spinner mx-auto mb-4"></div>
        <p className="text-white/80 font-medium">Analyzing therapeutic profile...</p>
        <p className="text-white/50 text-sm mt-1">Consulting 505+ clinical studies</p>
      </div>
    )
  }

  if (recommendations.length === 0) {
    return (
      <div className="text-center py-8">
        <div className="text-4xl mb-3">🌿</div>
        <p className="text-white/60">Select conditions to see therapeutic insights</p>
      </div>
    )
  }

  return (
    <div className="space-y-4">
      {/* Therapeutic Insights */}
      {recommendations.map((product: Recommendation, index: number) => (
        <div
          key={product.product_id || index}
          className="product-card overflow-hidden"
        >
          {/* Product Header */}
          <div className="flex justify-between items-start mb-3">
            <div className="flex-1">
              <div className="flex items-center gap-2 flex-wrap">
                <h3 className="font-bold text-gray-800 text-lg">{product.product_name}</h3>
                {product.confidence_score && (
                  <span className={getMatchBadgeClass(product.confidence_score)}>
                    {Math.round(product.confidence_score * 100)}% Match
                  </span>
                )}
              </div>
              <p className="text-emerald-600 text-sm font-medium">{product.product_type}</p>
            </div>
          </div>

          {/* Cannabinoid Profile */}
          <div className="mb-4 bg-gradient-to-r from-gray-50 to-gray-100 rounded-xl p-4">
            <div className="text-sm font-medium text-gray-600 mb-2">Cannabinoid Profile</div>
            <div className="flex gap-4">
              <div className="flex-1">
                <div className="flex justify-between text-xs text-gray-500 mb-1">
                  <span className="font-semibold text-orange-600">THC</span>
                  <span>{product.thc_percent}%</span>
                </div>
                <div className="h-3 bg-gray-200 rounded-full overflow-hidden">
                  <div 
                    className="bg-gradient-to-r from-orange-400 to-red-500 h-full rounded-full transition-all duration-500"
                    style={{ width: `${Math.min(product.thc_percent * 4, 100)}%` }}
                  />
                </div>
              </div>
              <div className="flex-1">
                <div className="flex justify-between text-xs text-gray-500 mb-1">
                  <span className="font-semibold text-blue-600">CBD</span>
                  <span>{product.cbd_percent}%</span>
                </div>
                <div className="h-3 bg-gray-200 rounded-full overflow-hidden">
                  <div 
                    className="bg-gradient-to-r from-blue-400 to-purple-500 h-full rounded-full transition-all duration-500"
                    style={{ width: `${Math.min(product.cbd_percent * 4, 100)}%` }}
                  />
                </div>
              </div>
            </div>
          </div>

          {/* Recommended Terpenes */}
          {product.terpene_profile && (
            <div className="mb-4">
              <div className="text-sm font-medium text-gray-600 mb-2">Recommended Terpenes</div>
              <div className="space-y-2">
                {product.terpene_profile.map((terpene, i) => {
                  const info = TERPENE_INFO[terpene] || { emoji: '🌿', color: 'terpene-default', effect: '', description: '' }
                  return (
                    <div key={i} className="flex items-center gap-3 bg-white rounded-lg p-3 border border-gray-100">
                      <span className="text-2xl">{info.emoji}</span>
                      <div className="flex-1">
                        <div className="font-semibold text-gray-800">{terpene}</div>
                        <div className="text-gray-500 text-xs">{info.description}</div>
                      </div>
                      <span className={`${info.color} text-xs`}>{info.effect}</span>
                    </div>
                  )
                })}
              </div>
            </div>
          )}

          {/* Therapeutic Rationale */}
          {product.therapeutic_rationale && (
            <div className="info-box-clinical mb-4">
              <div className="flex items-center gap-2 font-bold text-purple-800 mb-2">
                <span>🧬</span> Clinical Rationale
              </div>
              <p className="text-purple-700 text-sm">{product.therapeutic_rationale}</p>
              {product.clinical_studies && (
                <div className="mt-3 flex items-center gap-2 text-purple-600 text-xs">
                  <span>📊</span>
                  <span>Based on <strong>{product.clinical_studies}</strong> peer-reviewed studies</span>
                </div>
              )}
            </div>
          )}

          {/* Synergistic Adjuvants - Expanded Education */}
          {product.adjuvants && product.adjuvants.length > 0 && (
            <div className="bg-gradient-to-r from-blue-50 to-indigo-50 rounded-xl p-4 border border-blue-200">
              <div className="font-bold text-blue-800 mb-3 flex items-center gap-2">
                <span className="text-xl">💊</span> 
                Synergistic Adjuvants
                <span className="text-xs font-normal text-blue-600 ml-auto">Tap for details</span>
              </div>
              <p className="text-blue-700 text-sm mb-3">
                These supplements can enhance therapeutic outcomes when used alongside cannabis:
              </p>
              <div className="space-y-2">
                {product.adjuvants.map((adj, i) => (
                  <div key={i}>
                    <button
                      onClick={() => setExpandedAdjuvant(
                        expandedAdjuvant === adj.name ? null : adj.name
                      )}
                      className="w-full text-left bg-white rounded-lg p-3 border border-blue-100 hover:border-blue-300 transition-all"
                    >
                      <div className="flex items-center justify-between">
                        <div className="flex items-center gap-2">
                          <span className="text-xl">💎</span>
                          <span className="font-semibold text-blue-800">{adj.name}</span>
                        </div>
                        <svg 
                          className={`w-5 h-5 text-blue-400 transition-transform ${expandedAdjuvant === adj.name ? 'rotate-180' : ''}`} 
                          fill="none" 
                          stroke="currentColor" 
                          viewBox="0 0 24 24"
                        >
                          <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M19 9l-7 7-7-7" />
                        </svg>
                      </div>
                    </button>
                    
                    {expandedAdjuvant === adj.name && (
                      <div className="mt-2 bg-white rounded-lg p-4 border border-blue-100 space-y-3">
                        <div>
                          <div className="text-xs font-semibold text-blue-600 uppercase mb-1">How It Works</div>
                          <p className="text-gray-700 text-sm">{adj.mechanism}</p>
                        </div>
                        <div className="grid grid-cols-2 gap-3">
                          <div>
                            <div className="text-xs font-semibold text-emerald-600 uppercase mb-1">Suggested Dosage</div>
                            <p className="text-gray-700 text-sm">{adj.dosage}</p>
                          </div>
                          <div>
                            <div className="text-xs font-semibold text-purple-600 uppercase mb-1">Timing</div>
                            <p className="text-gray-700 text-sm">{adj.timing}</p>
                          </div>
                        </div>
                        <div>
                          <div className="text-xs font-semibold text-amber-600 uppercase mb-1">Evidence</div>
                          <p className="text-gray-700 text-sm">{adj.evidence}</p>
                        </div>
                      </div>
                    )}
                  </div>
                ))}
              </div>
            </div>
          )}
        </div>
      ))}
    </div>
  )
}
