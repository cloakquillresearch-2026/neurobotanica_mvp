export type ExperienceLevel = 'naive' | 'beginner' | 'intermediate' | 'regular' | 'experienced'

export interface BiomarkerValues {
  tnf_alpha?: number
  il6?: number
  crp?: number
  il1b?: number
  [key: string]: number | undefined
}

export interface CustomerProfileData {
  customer_id: string
  first_name?: string
  last_name?: string
  phone?: string
  email?: string
  age?: number | null
  notes?: string
  conditions: string[]
  experience_level: ExperienceLevel
  biomarkers?: BiomarkerValues
  last_visit?: string
  isNew?: boolean
  isSandbox?: boolean
  selected_compounds?: string[]
  tier?: string
  gender?: string
  weight?: number
  demographics?: Record<string, unknown>
}

export interface ConditionPayload {
  name: string
  severity: number
  is_primary: boolean
}

export interface CustomerProfilePayload {
  first_name?: string
  last_name?: string
  phone?: string
  email?: string
  age: number | null
  biological_sex: string
  weight_kg: number | null
  conditions: ConditionPayload[]
  experience_level: ExperienceLevel
  administration_preferences: string[]
  primary_goal: string
  biomarkers: Record<string, number>
  notes?: string
}

export interface TransactionItemPayload {
  product_id: string
  quantity: number
  price: number
  thc_percent?: number
  cbd_percent?: number
  product_name?: string
}

export interface TransactionPayload {
  customer_id: string
  items: TransactionItemPayload[]
  total_amount: number
  timestamp: string
  notes?: string
  status?: string
}

export interface ConditionProfilePayload {
  conditions: ConditionPayload[]
  experience_level: ExperienceLevel | string
  administration_preferences: string[]
  primary_goal: string
}

export interface InflammatorySynergyPayload {
  biomarkers: Record<string, number>
  condition_profile: ConditionProfilePayload
  available_kingdoms: string[]
}

export interface PatientProfilePayload {
  age?: number
  sex?: string
  conditions?: ConditionPayload[]
  notes?: string
}
