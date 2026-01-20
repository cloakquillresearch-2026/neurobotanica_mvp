// filepath: src/hooks/useNeuroBotanicaAnalysis.ts
import { useState, useCallback } from 'react'
import { neurobotanicaAPI } from '@/utils/api'

interface AnalysisResult {
  interactions: {
    warnings: string[]
    total_warnings: number
  }
  bias_correction: {
    adjusted_dose_mg: number
    factors_applied: Record<string, unknown>
  }
  synergy: {
    synergy_score: number
    tk_enhanced: boolean
  }
  plant_profile: Record<string, unknown>
  polysaccharide_effects: {
    effects: string
    confidence: number
  }
  processing_time_ms: number
}

interface UseNeuroBotanicaAnalysisReturn {
  analyze: (compoundIds: string[], demographics?: Record<string, unknown>, tier?: string, plantId?: string) => Promise<void>
  result: AnalysisResult | null
  loading: boolean
  error: string | null
  tkConsentRequired: boolean
}

export function useNeuroBotanicaAnalysis(): UseNeuroBotanicaAnalysisReturn {
  const [result, setResult] = useState<AnalysisResult | null>(null)
  const [loading, setLoading] = useState(false)
  const [error, setError] = useState<string | null>(null)
  const [tkConsentRequired, setTkConsentRequired] = useState(false)

  const analyze = useCallback(async (
    compoundIds: string[],
    demographics: Record<string, unknown> = {},
    tier: string = 'computational_only',
    plantId?: string
  ) => {
    setLoading(true)
    setError(null)
    setTkConsentRequired(false)

    try {
      const analysisResult = await neurobotanicaAPI.analyzeCompounds(compoundIds, demographics, tier, plantId)
      
      if (analysisResult.error) {
        if (analysisResult.error.includes('TK consent denied')) {
          setTkConsentRequired(true)
          setError('Traditional knowledge access requires consent. Please provide consent to access enhanced features.')
        } else {
          setError(analysisResult.details || analysisResult.error)
        }
        setResult(null)
      } else {
        setResult(analysisResult)
      }
    } catch (err) {
      setError(err instanceof Error ? err.message : 'Analysis failed')
      setResult(null)
    } finally {
      setLoading(false)
    }
  }, [])

  return {
    analyze,
    result,
    loading,
    error,
    tkConsentRequired
  }
}