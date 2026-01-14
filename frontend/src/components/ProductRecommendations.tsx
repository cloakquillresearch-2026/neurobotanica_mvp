import { useState, useEffect, useCallback } from 'react'
import { dispensaryAPI } from '@/utils/api'
import type {
  CustomerProfileData,
  ConditionPayload,
  InflammatorySynergyPayload,
} from '@/types/customer'

interface ProductRecommendationsProps {
  customer: CustomerProfileData | null
  recommendations: Recommendation[]
  onRecommendationsUpdate: (recommendations: Recommendation[]) => void
}

export interface Recommendation {
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
  // D1-backed fields
  condition?: string
  category?: string
  evidence_summary?: {
    total_studies: number
    study_types: string[]
    avg_confidence: number
  }
  recommended_cannabinoids?: Array<{
    cannabinoid: string
    evidence_count: number
    avg_confidence: number
  }>
  recommended_ratio?: string
  delivery_methods?: string[]
  dosing_guidance?: string | DosingGuidance
  citations?: Array<{
    study_id: string
    study_type: string
    citation: string
    confidence_score: number
    key_findings: string[]
  }>
  disclaimer?: string
}

interface AdjuvantInfo {
  name: string
  mechanism: string
  dosage: string
  timing: string
  evidence: string
}

interface DosingGuidance {
  primary?: string
  preferred_route?: string
  frequency?: string
  notes?: string[]
  history_summary?: string
}

interface SynergyResponse {
  primary_kingdom: string
  secondary_kingdoms: string[]
  synergy_score: number
  confidence_level: number
  recommended_compounds: string[]
  dosing_guidance?: DosingGuidance
  expected_reduction?: Record<string, string | number>
  warning?: string
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
  'Nutmeg': {
    name: 'Nutmeg',
    mechanism: 'Myristicin modulates serotonin receptors, provides anxiogenic and antidepressant effects',
    dosage: '100-300mg daily (low dose)',
    timing: 'Take with meals, avoid high doses',
    evidence: 'Multiple clinical studies show psychotropic and anti-inflammatory effects'
  },
  'Cinnamon': {
    name: 'Cinnamon',
    mechanism: 'Cinnamaldehyde inhibits inflammatory pathways, improves insulin sensitivity',
    dosage: '500-2000mg daily',
    timing: 'Take with meals',
    evidence: 'Meta-analyses confirm glycemic control and cardiovascular benefits'
  },
  'Mace': {
    name: 'Mace',
    mechanism: 'Similar to nutmeg with anti-inflammatory and gastroprotective effects',
    dosage: '100-300mg daily (low dose)',
    timing: 'Take with meals',
    evidence: 'Clinical studies support digestive and anti-inflammatory benefits'
  },
}

const SANDBOX_SYNERGY: SynergyResponse = {
  primary_kingdom: 'Cannabis',
  secondary_kingdoms: ['Fungal', 'Plant'],
  synergy_score: 0.82,
  confidence_level: 0.78,
  recommended_compounds: ['CBD', 'CBG', 'Linalool', 'Myrcene'],
  dosing_guidance: {
    primary: 'gentle_titration',
    preferred_route: 'sublingual_tincture',
    frequency: 'twice_daily',
    notes: [
      'Use tincture micro-doses to demonstrate slow onboarding',
      'Pair with slow diaphragmatic breathing for anxiety spikes',
    ],
    history_summary: 'Sandbox persona with inflammatory biomarkers + anxiety'
  },
  expected_reduction: {
    TNF_alpha: '-28% (6 weeks)',
    IL6: '-22% (8 weeks)',
    CRP: '-18% (6 weeks)',
  },
  warning: 'Sandbox data for training only. Do not dispense based on this panel.'
}

export function ProductRecommendations({
  customer,
  recommendations,
  onRecommendationsUpdate,
}: ProductRecommendationsProps) {
  const [loading, setLoading] = useState(false)
  const [expandedAdjuvant, setExpandedAdjuvant] = useState<string | null>(null)
  const [synergyData, setSynergyData] = useState<SynergyResponse | null>(null)
  const [synergyLoading, setSynergyLoading] = useState(false)
  const [synergyError, setSynergyError] = useState<string | null>(null)

  const fetchRecommendations = useCallback(async () => {
    if (!customer || customer.isNew) {
      return
    }
    setLoading(true)
    try {
      const conditions = customer.conditions || []

      // Production: call D1-backed recommendations endpoint
      if (conditions.length > 0) {
        try {
          const allRecommendations: Recommendation[] = []
          
          // Fetch recommendation for EACH condition
          for (const condition of conditions) {
            try {
              const resp = await dispensaryAPI.getRecommendations({ 
                condition: condition.toUpperCase(), 
                severity: 'moderate' 
              })
              const data = resp.data
              
              const rec: Recommendation = {
                product_id: `${condition.toUpperCase()}_REC`,
                product_name: `${data.condition} Evidence-based Recommendation`,
                product_type: data.category || 'Evidence-backed recommendation',
                thc_percent: 0,
                cbd_percent: 0,
                confidence_score: data.confidence_score,
                therapeutic_rationale: data.disclaimer || undefined,
                terpene_profile: undefined,
                clinical_studies: data.evidence_summary?.total_studies,
                adjuvants: undefined,
                condition: data.condition,
                category: data.category,
                evidence_summary: data.evidence_summary,
                recommended_cannabinoids: data.recommended_cannabinoids,
                recommended_ratio: data.recommended_ratio,
                delivery_methods: data.delivery_methods,
                dosing_guidance: data.dosing_guidance,
                citations: data.citations,
                disclaimer: data.disclaimer,
              }
              
              allRecommendations.push(rec)
            } catch (err) {
              console.error(`Failed to fetch recommendation for ${condition}:`, err)
              // Continue with other conditions even if one fails
            }
          }
          
          if (allRecommendations.length > 0) {
            onRecommendationsUpdate(allRecommendations)
            setLoading(false)
            return
          }
        } catch (err) {
          console.error('D1 recommendations request failed:', err)
        }
      }

      // Fallback to mock if API not available or failed
      const fallback: Recommendation[] = [{
        product_id: 'FALLBACK_001',
        product_name: 'Evidence lookup unavailable',
        product_type: 'Fallback',
        thc_percent: 0,
        cbd_percent: 0,
        confidence_score: 0,
        therapeutic_rationale: 'Unable to retrieve evidence-backed recommendations at this time.',
      }]

      onRecommendationsUpdate(fallback)
    } catch (error) {
      console.error('Failed to fetch recommendations:', error)
    } finally {
      setLoading(false)
    }
  }, [customer, onRecommendationsUpdate])

  const fetchSynergyInsights = useCallback(async () => {
    if (!customer) {
      setSynergyData(null)
      return
    }

    setSynergyLoading(true)
    setSynergyError(null)

    const biomarkerPayload = Object.entries(customer.biomarkers || {})
      .reduce<Record<string, number>>((acc, [key, value]) => {
        const parsed = typeof value === 'string' ? parseFloat(value) : value
        if (!Number.isNaN(parsed ?? NaN)) {
          acc[key] = Number(parsed)
        }
        return acc
      }, {})

    const normalizedConditions: ConditionPayload[] = (customer.conditions || []).map((condition, index) => ({
      name: condition,
      severity: 7,
      is_primary: index === 0,
    }))

    if (normalizedConditions.length === 0 && Object.keys(biomarkerPayload).length === 0) {
      setSynergyData(null)
      return
    }

    setSynergyLoading(true)
    setSynergyError(null)
    try {
      const payload: InflammatorySynergyPayload = {
        biomarkers: biomarkerPayload,
        condition_profile: {
          conditions: normalizedConditions,
          experience_level: customer.experience_level || 'beginner',
          administration_preferences: ['inhalation'],
          primary_goal: normalizedConditions[0]?.name || 'inflammation',
        },
        available_kingdoms: ['cannabis', 'fungal', 'marine', 'plant'],
      }
      const response = await dispensaryAPI.predictInflammatorySynergy(payload)
      setSynergyData(response.data)
    } catch (error) {
      const message = error instanceof Error ? error.message : 'Unable to fetch TS-PS-001 insights right now.'
      setSynergyError(message)
    } finally {
      setSynergyLoading(false)
    }
  }, [customer])

  useEffect(() => {
    fetchRecommendations()
  }, [fetchRecommendations])

  useEffect(() => {
    fetchSynergyInsights()
  }, [fetchSynergyInsights])

  const getMatchBadgeClass = (score: number) => {
    if (score >= 0.9) return 'match-badge-high'
    if (score >= 0.75) return 'match-badge-medium'
    return 'match-badge-low'
  }

  const renderSynergyPanel = () => {
    if (!customer) return null

    return (
      <section
        className="bg-white/5 border border-white/10 rounded-2xl p-4"
        aria-live="polite"
        aria-label="TS-PS-001 insights"
        data-testid="synergy-panel"
      >
        <div className="flex items-center gap-2 mb-3">
          <span className="text-xl" aria-hidden="true">🧬</span>
          <div>
            <p className="font-semibold text-white">TS-PS-001 Cross-Kingdom Insights</p>
            <p className="text-white/60 text-xs">Powered by inflammatory biomarkers + patient history</p>
          </div>
          {synergyLoading && (
            <span className="ml-auto flex items-center gap-2 text-emerald-300 text-xs" role="status">
              <span className="spinner w-4 h-4 border-2 border-emerald-200 border-t-transparent" aria-hidden="true" />
              Computing…
            </span>
          )}
        </div>

        {synergyError && (
          <div className="info-box-warning text-sm">
            <strong className="block">Unable to load TS-PS-001 insights.</strong>
            <span>{synergyError}</span>
          </div>
        )}

        {!synergyError && !synergyLoading && !synergyData && (
          <p className="text-white/70 text-sm">
            Add biomarkers or conditions to unlock personalized cross-kingdom guidance.
          </p>
        )}

        {synergyData && !synergyError && (
          <div className="space-y-4 text-white/90">
            <div className="grid gap-3 md:grid-cols-3">
              <div className="bg-slate-900/60 rounded-xl p-3 border border-white/10">
                <p className="text-xs text-white/50 uppercase font-semibold">Primary Kingdom</p>
                <p className="text-lg font-bold text-white text-elevated">{synergyData.primary_kingdom}</p>
                {synergyData.secondary_kingdoms?.length > 0 && (
                  <p className="text-xs text-white/60 mt-1">
                    Adjacent: {synergyData.secondary_kingdoms.join(', ')}
                  </p>
                )}
              </div>
              <div className="bg-slate-900/60 rounded-xl p-3 border border-white/10">
                <p className="text-xs text-white/50 uppercase font-semibold">Synergy Score</p>
                <p className="text-lg font-bold text-emerald-300" data-testid="synergy-score">
                  {(synergyData.synergy_score * 100).toFixed(1)}%
                </p>
                <p className="text-xs text-white/60 mt-1">Confidence {Math.round((synergyData.confidence_level || 0) * 100)}%</p>
              </div>
              <div className="bg-slate-900/60 rounded-xl p-3 border border-white/10">
                <p className="text-xs text-white/50 uppercase font-semibold">Dosing Strategy</p>
                <p className="text-sm font-semibold text-white">
                  {synergyData.dosing_guidance?.primary?.replace('_', ' ') || 'Standard titration'}
                </p>
                <p className="text-xs text-white/60 mt-1 capitalize">
                  {synergyData.dosing_guidance?.preferred_route?.replace('_', ' ') || 'oral_tincture'} • {synergyData.dosing_guidance?.frequency || 'twice_daily'}
                </p>
              </div>
            </div>

            {synergyData.recommended_compounds?.length > 0 && (
              <div>
                <p className="text-xs text-white/50 uppercase font-semibold mb-2">Key Compounds</p>
                <div className="flex flex-wrap gap-2">
                  {synergyData.recommended_compounds.map((compound: string) => (
                    <span key={compound} className="bg-white/10 border border-white/10 px-3 py-1 rounded-full text-xs">
                      {compound}
                    </span>
                  ))}
                </div>
              </div>
            )}

            {synergyData.expected_reduction && (
              <div className="grid gap-3 md:grid-cols-2">
                <div className="bg-slate-900/60 rounded-xl p-3 border border-white/10">
                  <p className="text-xs text-white/50 uppercase font-semibold mb-2">Projected Biomarker Reduction</p>
                  <ul className="space-y-1 text-sm">
                    {Object.entries(synergyData.expected_reduction).map(([marker, value]) => (
                      <li key={marker} className="flex justify-between text-white/80">
                        <span className="uppercase tracking-wide text-xs">{marker}</span>
                        <span>{value}</span>
                      </li>
                    ))}
                  </ul>
                </div>
                <div className="bg-slate-900/60 rounded-xl p-3 border border-white/10">
                  <p className="text-xs text-white/50 uppercase font-semibold mb-2">Care Notes</p>
                  <ul className="text-sm space-y-1 text-white/80">
                    {synergyData.dosing_guidance?.notes?.length
                      ? synergyData.dosing_guidance.notes.map((note: string, idx: number) => (
                        <li key={idx} className="flex gap-2">
                          <span className="text-emerald-300">•</span>
                          {note}
                        </li>
                      ))
                      : <li>No additional cautions.</li>
                    }
                    {synergyData.dosing_guidance?.history_summary && (
                      <li className="text-white/60 text-xs mt-2">History focus: {synergyData.dosing_guidance.history_summary}</li>
                    )}
                  </ul>
                </div>
              </div>
            )}

            {synergyData.warning && (
              <div className="info-box-warning text-sm">
                {synergyData.warning}
              </div>
            )}
          </div>
        )}
      </section>
    )
  }

  const renderRecommendationBody = () => {
    if (loading) {
      return (
        <div className="text-center py-12">
          <div className="spinner mx-auto mb-4" />
          <p className="text-white/80 font-medium">Analyzing therapeutic profile...</p>
          <p className="text-white/50 text-sm mt-1">Consulting 505+ clinical studies</p>
        </div>
      )
    }

    if (recommendations.length === 0) {
      return (
        <div className="text-center py-8">
          <div className="text-4xl mb-3" aria-hidden="true">🌿</div>
          <p className="text-white/60">Select conditions to see therapeutic insights</p>
        </div>
      )
    }

    return (
      <div className="space-y-4">
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
              {product.recommended_cannabinoids && product.recommended_cannabinoids.length > 0 ? (
                <div className="w-full">
                  <div className="text-xs text-gray-500 mb-2">Recommended Cannabinoids</div>
                  <div className="flex flex-wrap gap-2">
                    {product.recommended_cannabinoids.map((c, i) => (
                      <div key={i} className="bg-white rounded-lg p-3 border border-gray-100 text-sm">
                        <div className="font-semibold text-gray-800">{c.cannabinoid}</div>
                        <div className="text-gray-500 text-xs">{c.evidence_count} studies • avg {Math.round((c.avg_confidence || 0) * 100)}%</div>
                      </div>
                    ))}
                  </div>
                </div>
              ) : (
                <>
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
                </>
              )}
            </div>
          </div>

          {/* Evidence Summary (from D1) */}
          {product.evidence_summary && (
            <div className="mb-4 bg-gradient-to-r from-emerald-50 to-teal-50 rounded-xl p-4 border border-emerald-200">
              <div className="flex items-center gap-2 font-bold text-emerald-800 mb-2">
                <span>📊</span> Evidence Summary
              </div>
              <div className="text-emerald-700 text-sm">
                <strong>{product.evidence_summary.total_studies}</strong> clinical studies • {product.evidence_summary.study_types.join(', ')}
              </div>
              <div className="text-emerald-600 text-sm mt-1">
                Average confidence: <strong>{Math.round(product.evidence_summary.avg_confidence * 100)}%</strong>
              </div>
            </div>
          )}

          {/* Recommended Ratio / Delivery / Dosing */}
          {(product.recommended_ratio || product.delivery_methods || product.dosing_guidance) && (
            <div className="mb-4 bg-gradient-to-r from-gray-50 to-gray-100 rounded-xl p-4">
              {product.recommended_ratio && (
                <div className="text-sm text-gray-700 mb-2">
                  <span className="font-semibold text-gray-600">Recommended Ratio:</span> <strong className="text-gray-900">{product.recommended_ratio}</strong>
                </div>
              )}
              {product.delivery_methods && (
                <div className="text-sm text-gray-700 mb-2">
                  <span className="font-semibold text-gray-600">Delivery:</span> {product.delivery_methods.join(', ')}
                </div>
              )}
              {product.dosing_guidance && (
                <div className="text-sm text-gray-700">
                  <span className="font-semibold text-gray-600">Dosing:</span> {typeof product.dosing_guidance === 'string' ? product.dosing_guidance : (product.dosing_guidance as DosingGuidance).primary}
                </div>
              )}
            </div>
          )}

          {/* Citations */}
          {product.citations && product.citations.length > 0 && (
            <div className="mb-4 bg-gradient-to-r from-blue-50 to-indigo-50 rounded-xl p-4 border border-blue-200">
              <div className="flex items-center gap-2 font-bold text-blue-800 mb-3">
                <span>📚</span> Key Citations
              </div>
              <ul className="space-y-3">
                {product.citations.slice(0,5).map((c, i) => (
                  <li key={c.study_id || i} className="bg-white rounded-lg p-3 border border-blue-100">
                    <div className="font-semibold text-gray-800 text-sm">{c.citation}</div>
                    <div className="text-xs text-gray-500 mt-1">
                      {c.study_type} • Confidence <strong className="text-blue-600">{Math.round((c.confidence_score || 0) * 100)}%</strong>
                    </div>
                    {c.key_findings && c.key_findings.length > 0 && (
                      <div className="text-xs text-gray-600 mt-2 pl-3 border-l-2 border-blue-200">
                        {c.key_findings.slice(0,2).join(' • ')}
                      </div>
                    )}
                  </li>
                ))}
              </ul>
            </div>
          )}

          {/* Disclaimer (from D1) */}
          {product.disclaimer && (
            <div className="mb-4 bg-gradient-to-r from-amber-50 to-yellow-50 rounded-xl p-4 border border-amber-200">
              <div className="flex items-center gap-2 font-bold text-amber-800 mb-2">
                <span>⚠️</span> Important Information
              </div>
              <p className="text-amber-700 text-sm leading-relaxed">{product.disclaimer}</p>
            </div>
          )}

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

  return (
    <div className="space-y-5" aria-live="polite">
      {renderSynergyPanel()}
      {renderRecommendationBody()}
    </div>
  )
}
