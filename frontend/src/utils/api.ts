import axios from 'axios'
import type {
  CustomerProfilePayload,
  TransactionPayload,
  InflammatorySynergyPayload,
  PatientProfilePayload,
} from '@/types/customer'

const DEV_FALLBACK_API = 'http://127.0.0.1:8000'

const API_BASE_URL =
  process.env.NEXT_PUBLIC_API_URL ||
  (process.env.NODE_ENV !== 'production' ? DEV_FALLBACK_API : '')

export const api = axios.create({
  baseURL: API_BASE_URL,
  timeout: 10000,
})

// Dispensary API endpoints
export const dispensaryAPI = {
  // Customer management
  createProfile: (profile: CustomerProfilePayload) => api.post('/api/dispensary/profile', profile),
  getProfile: (profileId: string) => api.get(`/api/dispensary/profile/${profileId}`),
  updateProfile: (profileId: string, profile: CustomerProfilePayload) => api.put(`/api/dispensary/profile/${profileId}`, profile),

  // Transactions
  createTransaction: (transaction: TransactionPayload) => api.post('/api/dispensary/transaction', transaction),

  // Recommendations
  getRecommendations: (data: Record<string, unknown>) => api.post('/api/dispensary/recommend', data),
  predictInflammatorySynergy: (payload: InflammatorySynergyPayload) =>
    api.post('/api/dispensary/inflammatory-synergy', payload),

  // Feedback
  submitFeedback: (feedback: Record<string, unknown>) => api.post('/api/dispensary/feedback', feedback),

  // Analytics
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
const NEUROBOTANICA_API_BASE = 'https://neurobotanica-api.contessapetrini.workers.dev'

export const neurobotanicaAPI = {
  checkHealth: async () => {
    try {
      const response = await fetch(`${NEUROBOTANICA_API_BASE}/api/neurobotanica/health`);
      if (!response.ok) {
        throw new Error(`Health check failed: ${response.status} ${response.statusText}`);
      }
      return await response.json();
    } catch (error) {
      console.error('NeuroBotanica health check failed:', error);
      throw new Error('Unable to connect to NeuroBotanica API. Please check your internet connection.');
    }
  },
  
  analyzeCompounds: async (compoundIds: string[], demographics: Record<string, unknown> = {}, tier: string = 'computational_only', plantId?: string) => {
    try {
      const response = await fetch(`${NEUROBOTANICA_API_BASE}/api/neurobotanica/analyze`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ compound_ids: compoundIds, demographics, customer_tier: tier, plant_id: plantId })
      });
      
      if (!response.ok) {
        const errorText = await response.text();
        throw new Error(`Analysis failed: ${response.status} ${response.statusText} - ${errorText}`);
      }
      
      return await response.json();
    } catch (error) {
      console.error('NeuroBotanica analysis failed:', error);
      if (error instanceof TypeError && error.message.includes('fetch')) {
        throw new Error('Network error - unable to connect to NeuroBotanica API. Please check your internet connection.');
      }
      throw error;
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