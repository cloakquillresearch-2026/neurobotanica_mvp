import React from 'react'

interface Option {
  value: string
  label: string
}

interface Props {
  options: Option[]
  value?: string
  onChange: (value: string) => void
}

export function ConditionSelector({ options, value, onChange }: Props) {
  return (
    <div>
      <label className="block text-white/80 text-sm font-medium mb-2">Select Condition</label>
      <select
        value={value}
        onChange={(e) => onChange(e.target.value)}
        className="input-modern w-full"
        aria-label="Condition"
      >
        <option value="">-- Select condition --</option>
        {options.map((opt) => (
          <option key={opt.value} value={opt.value}>{opt.label}</option>
        ))}
      </select>
    </div>
  )
}

export default ConditionSelector
