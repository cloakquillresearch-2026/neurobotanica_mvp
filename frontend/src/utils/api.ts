import axios from 'axios'

const API_BASE_URL = process.env.NEXT_PUBLIC_API_URL || ''

export const api = axios.create({
  baseURL: API_BASE_URL,
  timeout: 10000,
})

// Dispensary API endpoints
export const dispensaryAPI = {
  // Customer management
  createProfile: (profile: any) => api.post('/api/dispensary/profile', profile),
  getProfile: (profileId: string) => api.get(`/api/dispensary/profile/${profileId}`),
  updateProfile: (profileId: string, profile: any) => api.put(`/api/dispensary/profile/${profileId}`, profile),

  // Transactions
  createTransaction: (transaction: any) => api.post('/api/dispensary/transaction', transaction),

  // Recommendations
  getRecommendations: (data: any) => api.post('/api/dispensary/recommend', data),

  // Feedback
  submitFeedback: (feedback: any) => api.post('/api/dispensary/feedback', feedback),

  // Analytics
  getStatistics: () => api.get('/api/dispensary/statistics'),
}

// Adjuvant optimization
export const adjuvantAPI = {
  optimize: (compound: string, target: string, kingdom: string = 'cannabis', patientProfile?: any) =>
    api.post('/api/adjuvant/optimize', {
      primary_compound: compound,
      therapeutic_target: target,
      kingdom: kingdom,
      patient_profile: patientProfile,
    }),
}

// General API utilities
export const handleApiError = (error: any) => {
  if (error.response) {
    // Server responded with error status
    throw new Error(error.response.data.detail || 'API request failed')
  } else if (error.request) {
    // Network error
    throw new Error('Network error - check your connection')
  } else {
    // Other error
    throw new Error(error.message || 'Unknown error occurred')
  }
}