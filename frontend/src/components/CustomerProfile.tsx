import { useState, useEffect } from 'react'

interface CustomerProfileProps {
  customer: any
  onProfileUpdate: (customer: any) => void
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

const EXPERIENCE_LEVELS = [
  { value: 'naive', label: 'First Time', emoji: '🌱', color: 'from-green-400 to-emerald-500' },
  { value: 'beginner', label: 'Beginner', emoji: '🌿', color: 'from-teal-400 to-cyan-500' },
  { value: 'intermediate', label: 'Occasional', emoji: '🍃', color: 'from-blue-400 to-indigo-500' },
  { value: 'regular', label: 'Regular', emoji: '🌳', color: 'from-purple-400 to-violet-500' },
  { value: 'experienced', label: 'Expert', emoji: '🏆', color: 'from-amber-400 to-orange-500' },
]

export function CustomerProfile({ customer, onProfileUpdate }: CustomerProfileProps) {
  const [isEditing, setIsEditing] = useState(customer?.isNew || false)
  const [editForm, setEditForm] = useState({
    first_name: customer?.first_name || '',
    last_name: customer?.last_name || '',
    age: customer?.age || '',
    conditions: customer?.conditions || [],
    experience_level: customer?.experience_level || 'beginner',
    notes: customer?.notes || '',
  })

  useEffect(() => {
    if (customer?.isNew) {
      setIsEditing(true)
    }
    setEditForm({
      first_name: customer?.first_name || '',
      last_name: customer?.last_name || '',
      age: customer?.age || '',
      conditions: customer?.conditions || [],
      experience_level: customer?.experience_level || 'beginner',
      notes: customer?.notes || '',
    })
  }, [customer])

  const handleSave = async () => {
    const updatedCustomer = {
      ...customer,
      ...editForm,
      age: parseInt(editForm.age as string) || undefined,
      isNew: false,
    }
    onProfileUpdate(updatedCustomer)
    setIsEditing(false)
  }

  const toggleCondition = (condition: string) => {
    const current = editForm.conditions as string[]
    if (current.includes(condition)) {
      setEditForm({ ...editForm, conditions: current.filter(c => c !== condition) })
    } else {
      setEditForm({ ...editForm, conditions: [...current, condition] })
    }
  }

  if (!customer) return null

  return (
    <div className="space-y-5">
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
            const isSelected = (editForm.conditions as string[]).includes(condition.value)
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
        className="btn-primary w-full flex items-center justify-center gap-2"
      >
        <span className="text-lg">🌿</span>
        NeuroBotanica Recommends
      </button>
    </div>
  )
}
