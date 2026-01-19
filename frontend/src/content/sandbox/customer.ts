// Minimal sandbox customer used to satisfy imports during local dev
import type { ExperienceLevel } from '@/types/customer'

export const SANDBOX_CUSTOMER = {
  customer_id: 'sandbox-customer',
  id: 'sandbox-customer',
  name: 'Test User',
  isSandbox: true,
  isNew: false,
  experience_level: 'beginner' as ExperienceLevel,
  conditions: ['anxiety'],
  biomarkers: {},
  preferences: {},
}

export default SANDBOX_CUSTOMER
