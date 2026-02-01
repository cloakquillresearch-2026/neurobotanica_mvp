import axios, { AxiosError } from 'axios'
import type {
  CustomerProfilePayload,
  TransactionPayload,
  InflammatorySynergyPayload,
  PatientProfilePayload,
} from '@/types/customer'

const DEV_FALLBACK_API = 'http://127.0.0.1:8000'

// API URLs for different services
const API_BASE_URL =
  process.env.NEXT_PUBLIC_API_URL ||
  (process.env.NODE_ENV !== 'production' ? DEV_FALLBACK_API : '')

// Cloudflare Worker URLs
const BUDTENDER_API_URL =
  process.env.NEXT_PUBLIC_BUDTENDER_API_URL ||
  'https://neurobotanica-api-production.contessapetrini.workers.dev'

export const api = axios.create({
  baseURL: API_BASE_URL,
  timeout: 10000,
})

// Budtender API proxy instance (for D1-backed endpoints)
export const budtenderApi = axios.create({
  baseURL: BUDTENDER_API_URL,
  timeout: 10000,
  headers: {
    'Content-Type': 'application/json',
  },
})

// Request interceptor for adding auth headers if needed
budtenderApi.interceptors.request.use(
  (config) => {
    // Add any auth headers here if needed
    return config
  },
  (error) => Promise.reject(error)
)

// Response interceptor for better error handling
budtenderApi.interceptors.response.use(
  (response) => response,
  (error: AxiosError) => {
    if (error.response) {
      // Server responded with error status
      const status = error.response.status
      const data = error.response.data as { error?: string; details?: string }

      if (status === 404) {
        throw new Error(data?.error || 'Resource not found')
      } else if (status === 400) {
        throw new Error(data?.error || 'Invalid request')
      } else if (status === 500) {
        throw new Error(data?.details || data?.error || 'Server error - please try again')
      } else if (status === 503) {
        throw new Error('Service temporarily unavailable - database may be down')
      }
    } else if (error.request) {
      // Network error
      throw new Error('Network error - please check your connection')
    }
    throw error
  }
)

// Dispensary API endpoints
export const dispensaryAPI = {
  // Customer management (Cloudflare Worker - D1 backed)
  searchProfiles: (query: string) => budtenderApi.get(`/api/dispensary/search?q=${encodeURIComponent(query)}`),
  createProfile: (profile: CustomerProfilePayload) => budtenderApi.post('/api/dispensary/profile', profile),
  getProfile: (profileId: string) => budtenderApi.get(`/api/dispensary/profile/${profileId}`),
  updateProfile: (profileId: string, profile: CustomerProfilePayload) => budtenderApi.put(`/api/dispensary/profile/${profileId}`, profile),

  // Transactions (local API)
  createTransaction: (transaction: TransactionPayload) => api.post('/api/dispensary/transaction', transaction),

  // D1-backed Recommendations (Cloudflare Worker)
  getRecommendations: async (data: Record<string, unknown>) => {
    try {
      // Use the D1-backed worker endpoint
      return await budtenderApi.post('/api/recommendations', data)
    } catch (error) {
      console.error('Recommendations fetch failed:', error)
      // Return a structured error response
      throw error
    }
  },

  // D1-backed Inflammatory Synergy (Cloudflare Worker)
  predictInflammatorySynergy: async (payload: InflammatorySynergyPayload) => {
    try {
      return await budtenderApi.post('/api/dispensary/inflammatory-synergy', payload)
    } catch (error) {
      console.error('Inflammatory synergy fetch failed:', error)
      throw error
    }
  },

  // Feedback (local API)
  submitFeedback: (feedback: Record<string, unknown>) => api.post('/api/dispensary/feedback', feedback),

  // Analytics (local API)
  getStatistics: () => api.get('/api/dispensary/statistics'),
}

// Adjuvant optimization
export const adjuvantAPI = {
  optimize: (
    compound: string,
    target: string,
    kingdom: string = 'cannabis',
    patientProfile?: PatientProfilePayload,
  ) =>
    api.post('/api/adjuvant/optimize', {
      primary_compound: compound,
      therapeutic_target: target,
      kingdom: kingdom,
      patient_profile: patientProfile,
    }),
}

// NeuroBotanica API endpoints (Cloudflare Workers)
const NEUROBOTANICA_API_BASE =
  process.env.NEXT_PUBLIC_NEUROBOTANICA_API_URL ||
  'https://terpene-api.contessapetrini.workers.dev'

// Retry configuration
const MAX_RETRIES = 3
const RETRY_DELAY_MS = 1000

async function fetchWithRetry(
  url: string,
  options: RequestInit,
  retries = MAX_RETRIES
): Promise<Response> {
  try {
    const response = await fetch(url, {
      ...options,
      signal: AbortSignal.timeout(10000), // 10 second timeout
    })
    return response
  } catch (error) {
    if (retries > 0 && error instanceof TypeError) {
      // Network error, retry after delay
      await new Promise(resolve => setTimeout(resolve, RETRY_DELAY_MS))
      return fetchWithRetry(url, options, retries - 1)
    }
    throw error
  }
}

export const neurobotanicaAPI = {
  checkHealth: async (): Promise<{
    status: string
    engines: string[]
    database?: { healthy: boolean; latency_ms: number }
  }> => {
    try {
      const response = await fetchWithRetry(
        `${NEUROBOTANICA_API_BASE}/api/neurobotanica/health`,
        { method: 'GET' }
      )

      if (!response.ok) {
        const errorData = await response.json().catch(() => ({}))
        throw new Error(
          (errorData as { error?: string }).error ||
          `Health check failed: ${response.status} ${response.statusText}`
        )
      }

      return await response.json()
    } catch (error) {
      console.error('NeuroBotanica health check failed:', error)
      if (error instanceof TypeError && error.message.includes('fetch')) {
        throw new Error('Network error - unable to connect to NeuroBotanica API')
      }
      throw new Error(
        error instanceof Error
          ? error.message
          : 'Unable to connect to NeuroBotanica API'
      )
    }
  },

  analyzeCompounds: async (
    compoundIds: string[],
    demographics: Record<string, unknown> = {},
    tier: string = 'computational_only',
    plantId?: string
  ): Promise<{
    interactions: { warnings: string[]; total_warnings: number }
    bias_correction: {
      adjusted_dose_mg: number
      factors_applied: Record<string, unknown>
      evidence: string
      confidence: number
      warnings?: string[]
    }
    synergy: { synergy_score: number; tk_enhanced: boolean; evidence: string }
    plant_profile: Record<string, unknown>
    polysaccharide_effects: { effects: string; confidence: number }
    processing_time_ms: number
    error?: string
    details?: string
  }> => {
    const startTime = Date.now()

    try {
      const response = await fetchWithRetry(
        `${NEUROBOTANICA_API_BASE}/api/neurobotanica/analyze`,
        {
          method: 'POST',
          headers: { 'Content-Type': 'application/json' },
          body: JSON.stringify({
            compound_ids: compoundIds,
            demographics,
            customer_tier: tier,
            plant_id: plantId,
          }),
        }
      )

      if (!response.ok) {
        let errorData: { error?: string; details?: string } = {}
        try {
          errorData = await response.json()
        } catch {
          // Ignore JSON parse errors
        }

        if (response.status === 400) {
          throw new Error(errorData.error || 'Invalid request parameters')
        } else if (response.status === 500) {
          throw new Error(errorData.details || errorData.error || 'Analysis server error')
        } else if (response.status === 503) {
          throw new Error('NeuroBotanica service temporarily unavailable')
        }

        throw new Error(
          errorData.error || `Analysis failed: ${response.status} ${response.statusText}`
        )
      }

      const result = await response.json()

      // Log performance for monitoring
      const elapsed = Date.now() - startTime
      if (elapsed > 5000) {
        console.warn(`NeuroBotanica analysis took ${elapsed}ms (target: <10s)`)
      }

      return result
    } catch (error) {
      console.error('NeuroBotanica analysis failed:', error)

      if (error instanceof TypeError && error.message.includes('fetch')) {
        throw new Error(
          'Network error - unable to connect to NeuroBotanica API. Please check your internet connection.'
        )
      }

      if (error instanceof DOMException && error.name === 'AbortError') {
        throw new Error('Request timed out - NeuroBotanica API is taking too long to respond')
      }

      throw error
    }
  },

  // New: Standalone bias correction endpoint
  getBiasCorrection: async (
    baseDose: number,
    compoundId: string,
    demographics: Record<string, unknown>,
    tier?: string
  ): Promise<{
    adjusted_dose_mg: number
    factors_applied: Record<string, unknown>
    evidence: string
    confidence: number
    warnings: string[]
  }> => {
    try {
      const response = await fetchWithRetry(
        `${NEUROBOTANICA_API_BASE}/api/neurobotanica/bias-correction`,
        {
          method: 'POST',
          headers: { 'Content-Type': 'application/json' },
          body: JSON.stringify({
            base_dose: baseDose,
            compound_id: compoundId,
            demographics,
            tier: tier || 'computational_only',
          }),
        }
      )

      if (!response.ok) {
        const errorData = await response.json().catch(() => ({}))
        throw new Error(
          (errorData as { error?: string }).error || 'Bias correction failed'
        )
      }

      return await response.json()
    } catch (error) {
      console.error('Bias correction failed:', error)
      throw error
    }
  },
}

// General API utilities
export const handleApiError = (error: unknown) => {
  if (error && typeof error === 'object' && 'response' in error && error.response) {
    const detail = (error as { response: { data?: { detail?: string } } }).response?.data?.detail
    // Server responded with error status
    throw new Error(detail || 'API request failed')
  }

  if (error && typeof error === 'object' && 'request' in error && (error as { request?: unknown }).request) {
    // Network error
    throw new Error('Network error - check your connection')
  }

  if (error instanceof Error) {
    // Other error
    throw new Error(error.message || 'Unknown error occurred')
  }

  throw new Error('Unknown error occurred')
}