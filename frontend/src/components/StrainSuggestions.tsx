import React, { useState } from 'react'

interface Strain {
  id: string
  name: string
  thc: number
  cbd: number
  confidence: number
  rationale?: string
}

interface Props {
  productId: string
  productName: string
}

// Simple mock strain lookup — will be replaced by real matching logic later
const MOCK_STRAINS: Strain[] = [
  { id: 'STRAIN_001', name: 'Calm CBD', thc: 2, cbd: 18, confidence: 0.94, rationale: 'High CBD with linalool dominant terpene profile matches anxiolytic evidence.' },
  { id: 'STRAIN_002', name: 'Gentle Relief', thc: 6, cbd: 12, confidence: 0.86, rationale: 'Balanced profile with myrcene for relaxation and sustained relief.' },
  { id: 'STRAIN_003', name: 'Night Ease', thc: 12, cbd: 4, confidence: 0.72, rationale: 'Higher THC + myrcene helps with sleep latency; use with caution for anxiety.' },
]

export function StrainSuggestions({ productName }: Props) {
  const strains = MOCK_STRAINS
  const [expanded, setExpanded] = useState<string | null>(null)

  return (
    <div className="mt-4 bg-white/5 border border-white/10 rounded-lg p-3">
      <div className="flex items-center justify-between mb-3">
        <div className="font-semibold text-white">Strain Matches for {productName}</div>
        <div className="text-xs text-white/60">Based on cannabinoid & terpene profile</div>
      </div>
      <div className="space-y-2">
        {strains.map((s) => (
          <div key={s.id} className="bg-gray-900/60 rounded-md p-3 border border-white/6">
            <div className="flex items-start justify-between">
              <div>
                <div className="flex items-center gap-2">
                  <div className="font-semibold text-white">{s.name}</div>
                  <div className="text-xs text-white/60">THC {s.thc}% • CBD {s.cbd}%</div>
                </div>
                <div className="text-white/60 text-sm mt-1">{s.rationale}</div>
              </div>
              <div className="text-right">
                <div className={`text-sm font-bold ${s.confidence >= 0.9 ? 'text-emerald-300' : s.confidence >= 0.75 ? 'text-amber-300' : 'text-white/60'}`}>
                  {Math.round(s.confidence * 100)}% match
                </div>
                <button
                  className="mt-2 btn-secondary text-xs"
                  onClick={() => setExpanded(expanded === s.id ? null : s.id)}
                >
                  {expanded === s.id ? 'Hide why' : 'Why this?'}
                </button>
              </div>
            </div>

            {expanded === s.id && (
              <div className="mt-3 bg-white/3 p-3 rounded-md text-sm text-white/70">
                <strong>Why this match:</strong>
                <p className="mt-1">{s.rationale} This score factors cannabinoid ratio, terpene overlap, and clinical confidence.</p>
              </div>
            )}
          </div>
        ))}
      </div>
    </div>
  )
}

export default StrainSuggestions
