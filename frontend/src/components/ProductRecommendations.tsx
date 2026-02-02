import { useState, useEffect, useCallback } from 'react'
import { dispensaryAPI } from '@/utils/api'
import { useNeuroBotanicaAnalysis } from '@/hooks/useNeuroBotanicaAnalysis'
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

const CONDITION_OPTIONS = [
  { value: 'chronic_pain', label: 'Chronic Pain', emoji: '🔥' },
  { value: 'anxiety', label: 'Anxiety', emoji: '🧘' },
  { value: 'insomnia', label: 'Insomnia', emoji: '🌙' },
  { value: 'ptsd', label: 'PTSD', emoji: '💚' },
  { value: 'nausea', label: 'Nausea', emoji: '🤢' },
  { value: 'inflammation', label: 'Inflammation', emoji: '🔴' },
  { value: 'depression', label: 'Depression', emoji: '☀️' },
  { value: 'weight_management', label: 'Weight/Metabolism', emoji: '⚖️' },
  { value: 'muscle_spasms', label: 'Muscle Spasms', emoji: '💪' },
  { value: 'seizures', label: 'Seizures', emoji: '⚡' },
]

const TERPENE_INFO: Record<string, { emoji: string; color: string; effect: string; description: string }> = {
  'Myrcene': { emoji: '🥭', color: 'terpene-myrcene', effect: 'Relaxation', description: 'Enhances THC absorption, promotes sedation' },
  'Limonene': { emoji: '🍋', color: 'terpene-limonene', effect: 'Mood', description: 'Elevates mood, relieves stress and anxiety' },
  'Pinene': { emoji: '🌲', color: 'terpene-pinene', effect: 'Focus', description: 'Improves alertness, counteracts THC memory effects' },
  'Linalool': { emoji: '💜', color: 'terpene-linalool', effect: 'Calm', description: 'Reduces anxiety, promotes relaxation' },
  'Caryophyllene': { emoji: '🌶️', color: 'terpene-caryophyllene', effect: 'Pain', description: 'Anti-inflammatory, activates CB2 receptors' },
  'Beta-Caryophyllene': { emoji: '🌶️', color: 'terpene-caryophyllene', effect: 'Pain', description: 'Anti-inflammatory, activates CB2 receptors' },
}

export function ProductRecommendations({
  customer,
  recommendations,
  onRecommendationsUpdate,
}: ProductRecommendationsProps) {
  const [loading, setLoading] = useState(false)
  const [loadingMessage, setLoadingMessage] = useState('Loading recommendations...')
  const [error, setError] = useState<string | null>(null)
  const [expandedAdjuvant, setExpandedAdjuvant] = useState<string | null>(null)
  const [synergyData, setSynergyData] = useState<SynergyResponse | null>(null)
  const [synergyLoading, setSynergyLoading] = useState(false)
  const [synergyError, setSynergyError] = useState<string | null>(null)
  const [hasRunAnalysis, setHasRunAnalysis] = useState(false)

  // NeuroBotanica analysis hook
  const { analyze: analyzeNeuroBotanica, result: neurobotanicaResult, loading: neurobotanicaLoading, error: neurobotanicaError, tkConsentRequired } = useNeuroBotanicaAnalysis()

  useEffect(() => {
    setHasRunAnalysis(false)
    setSynergyData(null)
    setSynergyError(null)
    setSynergyLoading(false)
  }, [customer])

  const fetchRecommendations = useCallback(async () => {
    if (!customer || customer.isNew) {
      return
    }
    setLoading(true)
    setError(null)
    setLoadingMessage('Connecting to clinical evidence database...')

    const startTime = Date.now()

    try {
      const conditions = customer.conditions || []

      // Production: call D1-backed recommendations endpoint
      if (conditions.length > 0) {
        const allRecommendations: Recommendation[] = []
        const failedConditions: string[] = []

        setLoadingMessage(`Analyzing ${conditions.length} condition${conditions.length > 1 ? 's' : ''}...`)

        // Fetch recommendation for EACH condition in parallel for speed
        const promises = conditions.map(async (condition) => {
          try {
            const resp = await dispensaryAPI.getRecommendations({
              condition: condition.toUpperCase(),
              severity: 'moderate'
            })
            const data = resp.data

            return {
              success: true,
              condition,
              recommendation: {
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
              } as Recommendation
            }
          } catch (err) {
            console.error(`Failed to fetch recommendation for ${condition}:`, err)
            return { success: false, condition, error: err }
          }
        })

        const results = await Promise.all(promises)

        results.forEach(result => {
          if (result.success && result.recommendation) {
            allRecommendations.push(result.recommendation)
          } else {
            failedConditions.push(result.condition)
          }
        })

        const elapsed = Date.now() - startTime
        if (elapsed > 5000) {
          console.warn(`Recommendations fetch took ${elapsed}ms (target: <10s)`)
        }

        if (allRecommendations.length > 0) {
          onRecommendationsUpdate(allRecommendations)

          // Show partial error if some conditions failed
          if (failedConditions.length > 0) {
            setError(`Some conditions could not be loaded: ${failedConditions.join(', ')}`)
          }
          setLoading(false)
          return
        }

        // All conditions failed
        if (failedConditions.length === conditions.length) {
          setError('Unable to load recommendations. Please check your connection and try again.')
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
      const errorMessage = error instanceof Error ? error.message : 'An unexpected error occurred'
      setError(errorMessage)
    } finally {
      setLoading(false)
    }
  }, [customer, onRecommendationsUpdate])

  // Retry handler for failed requests
  const handleRetry = useCallback(() => {
    setError(null)
    fetchRecommendations()
  }, [fetchRecommendations])

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
      setSynergyLoading(false)
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

  const runNeuroBotanicaAnalysis = useCallback(async () => {
    if (!customer) return

    setHasRunAnalysis(true)
    await fetchSynergyInsights()

    // Extract compound IDs from customer data (assuming they have selected compounds)
    const compoundIds = customer.selected_compounds || ['cbd', 'thc'] // Default fallback
    
    // Extract demographics
    const demographics = {
      age: customer.age || 30,
      gender: customer.gender || 'unknown',
      weight: customer.weight,
      ...customer.demographics
    }

    // Determine tier based on customer preferences or consent
    const tier = customer.tier || 'computational_only'

    await analyzeNeuroBotanica(compoundIds, demographics, tier)
  }, [customer, analyzeNeuroBotanica, fetchSynergyInsights])

  useEffect(() => {
    fetchRecommendations()
  }, [fetchRecommendations])


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
            {hasRunAnalysis
              ? 'Add biomarkers or conditions to unlock personalized cross-kingdom guidance.'
              : 'Run Analysis to generate TS-PS-001 insights.'}
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
        <div className="text-center py-12" role="status" aria-live="polite">
          <div className="spinner mx-auto mb-4" aria-hidden="true" />
          <p className="text-white/80 font-medium">{loadingMessage}</p>
          <p className="text-white/50 text-sm mt-1">Consulting 505+ clinical studies</p>
          <p className="text-white/40 text-xs mt-3">Target: &lt;10 seconds</p>
        </div>
      )
    }

    if (error && recommendations.length === 0) {
      return (
        <div className="text-center py-8" role="alert">
          <div className="text-4xl mb-3" aria-hidden="true">⚠️</div>
          <p className="text-red-400 font-medium mb-2">Unable to load recommendations</p>
          <p className="text-white/60 text-sm mb-4">{error}</p>
          <button
            onClick={handleRetry}
            className="px-4 py-2 bg-emerald-600 hover:bg-emerald-500 text-white rounded-lg font-medium transition-colors"
          >
            Try Again
          </button>
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

    // Show partial error banner if some conditions failed but we have results
    const errorBanner = error && recommendations.length > 0 ? (
      <div className="bg-amber-500/20 border border-amber-500/30 rounded-lg p-3 mb-4" role="alert">
        <div className="flex items-center gap-2">
          <span className="text-amber-400">⚠️</span>
          <p className="text-amber-200 text-sm">{error}</p>
          <button
            onClick={handleRetry}
            className="ml-auto text-amber-300 hover:text-amber-200 text-sm underline"
          >
            Retry
          </button>
        </div>
      </div>
    ) : null

    // Show all selected conditions at the top
    const selectedConditions = customer?.conditions || []

    return (
      <div className="space-y-4">
        {errorBanner}
        {selectedConditions.length > 0 && (
          <div className="bg-white/5 border border-white/10 rounded-xl p-4">
            <div className="flex items-center gap-2 mb-3">
              <span className="text-xl" aria-hidden="true">📋</span>
              <div>
                <p className="text-white font-semibold">Selected Conditions</p>
                <p className="text-white/60 text-sm">All conditions displayed below</p>
              </div>
            </div>
            <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-3">
              {selectedConditions.map((condition) => {
                const conditionOption = CONDITION_OPTIONS.find(opt => opt.value === condition)
                return (
                  <div
                    key={condition}
                    className="bg-emerald-500/20 border border-emerald-500/30 p-3 rounded-lg text-sm text-emerald-200"
                  >
                    <div className="flex items-center gap-2 mb-1">
                      <span>{conditionOption?.emoji || '🌿'}</span>
                      <span className="font-semibold">{conditionOption?.label || condition}</span>
                    </div>
                    <p className="text-emerald-200/80 text-xs">Condition selected for analysis</p>
                  </div>
                )
              })}
            </div>
          </div>
        )}

        {recommendations.map((product: Recommendation, index: number) => (
          <div
            key={product.product_id || index}
            className="product-card overflow-hidden"
          >
          {/* Product Header */}
          <div className="flex justify-between items-start mb-3">
            <div className="flex-1">
              <div className="flex items-center gap-2 flex-wrap">
                <h3 className="font-bold text-gray-800 text-lg">{product.product_name || 'Recommended Product'}</h3>
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

  const renderNeuroBotanicaPanel = () => {
    if (!customer) return null

    return (
      <section className="bg-gradient-to-br from-emerald-900 via-teal-900 to-cyan-900 rounded-2xl p-6 border border-emerald-500/20">
        <div className="flex items-center justify-between mb-6">
          <div className="flex items-center gap-3">
            <div className="w-10 h-10 bg-emerald-500/20 rounded-xl flex items-center justify-center">
              <span className="text-2xl">🧬</span>
            </div>
            <div>
              <h3 className="text-xl font-bold text-white">NeuroBotanica Analysis</h3>
              <p className="text-emerald-200/80 text-sm">AI-powered cannabis optimization</p>
            </div>
          </div>
          <button
            onClick={runNeuroBotanicaAnalysis}
            disabled={neurobotanicaLoading}
            className="px-4 py-2 bg-emerald-600 hover:bg-emerald-500 disabled:bg-emerald-600/50 text-white rounded-lg font-medium transition-colors"
          >
            {neurobotanicaLoading ? 'Analyzing...' : 'Run Analysis'}
          </button>
        </div>

        {tkConsentRequired && (
          <div className="bg-amber-500/20 border border-amber-500/30 rounded-lg p-4 mb-4">
            <div className="flex items-start gap-3">
              <span className="text-amber-400 text-xl">⚠️</span>
              <div>
                <h4 className="text-amber-200 font-semibold mb-1">Traditional Knowledge Access Required</h4>
                <p className="text-amber-200/80 text-sm">
                  Enhanced analysis requires consent for traditional knowledge integration. 
                  Please provide consent to access TK-enhanced features.
                </p>
              </div>
            </div>
          </div>
        )}

        {neurobotanicaError && !tkConsentRequired && (
          <div className="bg-red-500/20 border border-red-500/30 rounded-lg p-4 mb-4">
            <p className="text-red-200">{neurobotanicaError}</p>
          </div>
        )}

        {neurobotanicaResult && (
          <div className="grid gap-4 md:grid-cols-2 lg:grid-cols-3">
            {/* Interactions */}
            <div className="bg-slate-900/60 rounded-xl p-4 border border-white/10">
              <div className="flex items-center gap-2 mb-3">
                <span className="text-xl">💊</span>
                <h4 className="text-white font-semibold">Drug Interactions</h4>
              </div>
              <div className="text-2xl font-bold text-white mb-1">{neurobotanicaResult.interactions.total_warnings}</div>
              <p className="text-white/60 text-sm">warnings detected</p>
            </div>

            {/* Bias Correction */}
            <div className="bg-slate-900/60 rounded-xl p-4 border border-white/10">
              <div className="flex items-center gap-2 mb-3">
                <span className="text-xl">⚖️</span>
                <h4 className="text-white font-semibold">Bias Correction</h4>
              </div>
              <div className="text-2xl font-bold text-white mb-1">{neurobotanicaResult.bias_correction.adjusted_dose_mg}mg</div>
              <p className="text-white/60 text-sm">adjusted dose</p>
            </div>

            {/* Synergy */}
            <div className="bg-slate-900/60 rounded-xl p-4 border border-white/10">
              <div className="flex items-center gap-2 mb-3">
                <span className="text-xl">🔗</span>
                <h4 className="text-white font-semibold">Compound Synergy</h4>
              </div>
              <div className="text-2xl font-bold text-white mb-1">{(neurobotanicaResult.synergy.synergy_score * 100).toFixed(1)}%</div>
              <p className="text-white/60 text-sm">
                {neurobotanicaResult.synergy.tk_enhanced ? 'TK-enhanced' : 'computational'}
              </p>
            </div>

            {/* Polysaccharides */}
            <div className="bg-slate-900/60 rounded-xl p-4 border border-white/10 md:col-span-2 lg:col-span-3">
              <div className="flex items-center gap-2 mb-3">
                <span className="text-xl">🦠</span>
                <h4 className="text-white font-semibold">Microbiome Effects</h4>
              </div>
              <div className="flex items-center justify-between">
                <div>
                  <div className="text-lg font-bold text-white mb-1">{neurobotanicaResult.polysaccharide_effects.effects}</div>
                  <p className="text-white/60 text-sm">predicted modulation</p>
                </div>
                <div className="text-right">
                  <div className="text-2xl font-bold text-emerald-400">{(neurobotanicaResult.polysaccharide_effects.confidence * 100).toFixed(1)}%</div>
                  <p className="text-white/60 text-sm">confidence</p>
                </div>
              </div>
            </div>
          </div>
        )}

        <div className="mt-4 text-xs text-white/40">
          Processing time: {neurobotanicaResult?.processing_time_ms || 0}ms • TK Governance: Active
        </div>
      </section>
    )
  }

  return (
    <div className="space-y-5" aria-live="polite">
      {renderSynergyPanel()}
      {renderNeuroBotanicaPanel()}
      {renderRecommendationBody()}
    </div>
  )
}
