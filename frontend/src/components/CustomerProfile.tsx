import { useState, useEffect, useMemo } from 'react'
import { dispensaryAPI } from '@/utils/api'
import type {
  CustomerProfileData,
  ExperienceLevel,
  CustomerProfilePayload,
} from '@/types/customer'

interface CustomerProfileProps {
  customer: CustomerProfileData
  onProfileUpdate: (customer: CustomerProfileData) => void
}

interface EditFormState {
  first_name: string
  last_name: string
  age: string
  conditions: string[]
  experience_level: ExperienceLevel
  notes: string
  biomarkers: Record<string, string>
}

interface BiomarkerField {
  key: keyof CustomerProfilePayload['biomarkers'] | string
  label: string
  unit: string
  helper: string
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

const EXPERIENCE_LEVELS: Array<{
  value: ExperienceLevel
  label: string
  emoji: string
  color: string
}> = [
  { value: 'naive', label: 'First Time', emoji: '🌱', color: 'from-green-400 to-emerald-500' },
  { value: 'beginner', label: 'Beginner', emoji: '🌿', color: 'from-teal-400 to-cyan-500' },
  { value: 'intermediate', label: 'Occasional', emoji: '🍃', color: 'from-blue-400 to-indigo-500' },
  { value: 'regular', label: 'Regular', emoji: '🌳', color: 'from-purple-400 to-violet-500' },
  { value: 'experienced', label: 'Expert', emoji: '🏆', color: 'from-amber-400 to-orange-500' },
]

export function CustomerProfile({ customer, onProfileUpdate }: CustomerProfileProps) {
  const [saving, setSaving] = useState(false)
  const [editForm, setEditForm] = useState<EditFormState>(() => ({
    first_name: customer.first_name || '',
    last_name: customer.last_name || '',
    age: customer.age?.toString() || '',
    conditions: customer.conditions || [],
    experience_level: customer.experience_level || 'beginner',
    notes: customer.notes || '',
    biomarkers: {
      tnf_alpha: customer.biomarkers?.tnf_alpha?.toString() || '',
      il6: customer.biomarkers?.il6?.toString() || '',
      crp: customer.biomarkers?.crp?.toString() || '',
      il1b: customer.biomarkers?.il1b?.toString() || '',
    },
  }))

  useEffect(() => {
    setEditForm({
      first_name: customer.first_name || '',
      last_name: customer.last_name || '',
      age: customer.age?.toString() || '',
      conditions: customer.conditions || [],
      experience_level: customer.experience_level || 'beginner',
      notes: customer.notes || '',
      biomarkers: {
        tnf_alpha: customer.biomarkers?.tnf_alpha?.toString() || '',
        il6: customer.biomarkers?.il6?.toString() || '',
        crp: customer.biomarkers?.crp?.toString() || '',
        il1b: customer.biomarkers?.il1b?.toString() || '',
      },
    })
  }, [customer])

  const biomarkerFields = useMemo<BiomarkerField[]>(() => ([
    { key: 'tnf_alpha', label: 'TNF-α', unit: 'pg/mL', helper: 'Cytokine load linked to severe inflammation' },
    { key: 'il6', label: 'IL-6', unit: 'pg/mL', helper: 'Autoimmune and flare sensitivity marker' },
    { key: 'crp', label: 'CRP', unit: 'mg/L', helper: 'General inflammation + cardiovascular risk' },
    { key: 'il1b', label: 'IL-1β', unit: 'pg/mL', helper: 'Neuroinflammation + immune response' },
  ]), [])

  const handleBiomarkerChange = (field: string, value: string) => {
    setEditForm((prev) => ({
      ...prev,
      biomarkers: {
        ...prev.biomarkers,
        [field]: value,
      },
    }))
  }

  const handleSave = async () => {
    const biomarkerPayload = Object.entries(editForm.biomarkers)
      .reduce<Record<string, number>>((acc, [key, value]) => {
        const parsed = parseFloat(value)
        if (!Number.isNaN(parsed)) {
          acc[key] = parsed
        }
        return acc
      }, {})

    const buildUpdatedCustomer = (profileId?: string): CustomerProfileData => ({
      ...customer,
      ...editForm,
      age: parseInt(editForm.age || '', 10) || undefined,
      customer_id: profileId ?? customer.customer_id,
      isNew: false,
      conditions: editForm.conditions,
      biomarkers: biomarkerPayload,
    })

    try {
      setSaving(true)

      if (customer.isSandbox) {
        onProfileUpdate(buildUpdatedCustomer(customer.customer_id))
        return
      }

      // Prepare profile data for API
      const profileData: CustomerProfilePayload = {
        age: parseInt(editForm.age || '', 10) || null,
        biological_sex: 'unspecified', // Default for dispensary
        weight_kg: null, // Not collected in tablet interface
        conditions: editForm.conditions.map((condition, index) => ({
          name: condition,
          severity: 7, // Default medium severity
          is_primary: index === 0, // First condition is primary
        })),
        experience_level: editForm.experience_level,
        administration_preferences: ['inhalation'], // Default
        primary_goal: 'pain_relief', // Could be enhanced to map from conditions
        biomarkers: biomarkerPayload,
      }

      let savedProfile
      if (customer.customer_id && !customer.customer_id.startsWith('temp_')) {
        // Update existing profile
        savedProfile = await dispensaryAPI.updateProfile(customer.customer_id, profileData)
      } else {
        // Create new profile
        savedProfile = await dispensaryAPI.createProfile(profileData)
      }

      onProfileUpdate(buildUpdatedCustomer(savedProfile.data.profile_id))
    } catch (error) {
      console.error('Failed to save customer profile:', error)
      // Still update local state even if API fails
      onProfileUpdate(buildUpdatedCustomer())
    } finally {
      setSaving(false)
    }
  }

  const toggleCondition = (condition: string) => {
    const current = editForm.conditions
    if (current.includes(condition)) {
      setEditForm({ ...editForm, conditions: current.filter(c => c !== condition) })
    } else {
      setEditForm({ ...editForm, conditions: [...current, condition] })
    }
  }

  return (
    <div className="space-y-5">
      <div role="status" aria-live="polite" className="sr-only">
        {saving ? 'Saving profile' : 'Profile ready'}
      </div>
      {/* Experience Level Selection */}
      <div>
        <label className="block text-white/80 text-sm font-medium mb-3">
          Experience Level
        </label>
        <div className="grid grid-cols-5 gap-2">
          {EXPERIENCE_LEVELS.map((level) => (
            <button
              key={level.value}
              onClick={() => setEditForm({ ...editForm, experience_level: level.value })}
              className={`relative overflow-hidden rounded-xl p-0.5 transition-all duration-200 ${
                editForm.experience_level === level.value
                  ? `bg-gradient-to-r ${level.color} shadow-lg scale-105`
                  : 'bg-white/20 hover:bg-white/30'
              }`}
            >
              <div className={`rounded-[10px] p-3 text-center ${
                editForm.experience_level === level.value
                  ? 'bg-gray-900/80'
                  : 'bg-gray-900/60'
              }`}>
                <div className="text-2xl mb-1">{level.emoji}</div>
                <div className="text-white text-xs font-medium">{level.label}</div>
              </div>
            </button>
          ))}
        </div>
      </div>

      {/* Conditions Selection */}
      <div>
        <label className="block text-white/80 text-sm font-medium mb-3">
          Primary Conditions <span className="text-white/40">(tap to select)</span>
        </label>
        <div className="grid grid-cols-2 gap-2">
          {CONDITION_OPTIONS.map((condition) => {
            const isSelected = editForm.conditions.includes(condition.value)
            return (
              <button
                key={condition.value}
                onClick={() => toggleCondition(condition.value)}
                className={`flex items-center gap-3 p-3 rounded-xl transition-all duration-200 ${
                  isSelected
                    ? 'bg-gradient-to-r from-emerald-500 to-teal-500 text-white shadow-lg shadow-emerald-500/30'
                    : 'bg-white/10 text-white/80 hover:bg-white/20 border border-white/20'
                }`}
              >
                <span className="text-xl">{condition.emoji}</span>
                <span className="font-medium text-sm">{condition.label}</span>
                {isSelected && (
                  <svg className="w-5 h-5 ml-auto" fill="currentColor" viewBox="0 0 20 20">
                    <path fillRule="evenodd" d="M16.707 5.293a1 1 0 010 1.414l-8 8a1 1 0 01-1.414 0l-4-4a1 1 0 011.414-1.414L8 12.586l7.293-7.293a1 1 0 011.414 0z" clipRule="evenodd" />
                  </svg>
                )}
              </button>
            )
          })}
        </div>
      </div>

      {/* Age Input */}
      <div>
        <label className="block text-white/80 text-sm font-medium mb-2">
          Age <span className="text-white/40">(optional)</span>
        </label>
        <input
          type="number"
          value={editForm.age}
          onChange={(e) => setEditForm({ ...editForm, age: e.target.value })}
          placeholder="Enter age"
          className="input-modern"
        />
      </div>

      {/* Consultation Notes */}
      <div>
        <label className="block text-white/80 text-sm font-medium mb-2">
          Notes <span className="text-white/40">(optional)</span>
        </label>
        <textarea
          value={editForm.notes}
          onChange={(e) => setEditForm({ ...editForm, notes: e.target.value })}
          placeholder="Preferences, sensitivities, previous experiences..."
          rows={2}
          className="input-modern resize-none"
        />
      </div>

        {/* Biomarker Inputs */}
        <div className="bg-white/5 border border-white/15 rounded-2xl p-4">
          <div className="flex items-center gap-2 mb-4">
            <span className="text-xl" aria-hidden="true">🧪</span>
            <div>
              <p className="text-white font-semibold">Inflammatory Biomarkers</p>
              <p className="text-white/60 text-xs">Optional labs supercharge TS-PS-001 personalization</p>
            </div>
          </div>
          <div className="grid grid-cols-1 sm:grid-cols-2 gap-4">
            {biomarkerFields.map((field) => (
              <label key={field.key} className="text-sm text-white/80 space-y-1">
                <span className="flex items-center justify-between">
                  <span>{field.label}</span>
                  <span className="text-white/50 text-xs">{field.unit}</span>
                </span>
                <input
                  type="number"
                  inputMode="decimal"
                  min="0"
                  step="0.01"
                  value={editForm.biomarkers[field.key] || ''}
                  onChange={(event) => handleBiomarkerChange(field.key, event.target.value)}
                  className="input-modern"
                  aria-describedby={`${field.key}-helper`}
                />
                <span id={`${field.key}-helper`} className="text-white/50 text-xs block">
                  {field.helper}
                </span>
              </label>
            ))}
          </div>
          <div className="info-box-clinical mt-4 text-xs">
            Provide whatever lab data you have (even approximations). Missing values are fine—TS-PS-001 will use heuristics where necessary.
          </div>
        </div>

      {/* Warning for first-time users */}
      {editForm.experience_level === 'naive' && (
        <div className="info-box-warning">
          <div className="flex items-start gap-3">
            <span className="text-2xl">⚠️</span>
            <div>
              <div className="font-bold text-amber-800">First-Time User</div>
              <div className="text-amber-700 text-sm mt-1">
                Recommend <strong>low THC</strong> products. Consider CBD-dominant options or balanced 1:1 ratios. Start low, go slow.
              </div>
            </div>
          </div>
        </div>
      )}

      {/* Save Button */}
      <button
        onClick={handleSave}
        className="btn-primary w-full flex items-center justify-center gap-2 disabled:opacity-60 disabled:cursor-not-allowed"
        disabled={saving}
      >
        {saving ? (
          <>
            <span className="spinner w-5 h-5 border-2 border-white border-t-transparent" aria-hidden="true" />
            Saving profile…
          </>
        ) : (
          <>
            <span className="text-lg" aria-hidden="true">🌿</span>
            NeuroBotanica Recommends
          </>
        )}
      </button>
    </div>
  )
}
