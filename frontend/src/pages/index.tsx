import { useState, useMemo } from 'react'
import Head from 'next/head'
import { CustomerSearch } from '@/components/CustomerSearch'
import { ProductRecommendations } from '@/components/ProductRecommendations'
import { CustomerProfile } from '@/components/CustomerProfile'

// Condition-specific guidance data
const CONDITION_GUIDANCE: Record<string, {
  title: string
  emoji: string
  color: string
  cannabinoid: string
  terpenes: string[]
  tips: string[]
}> = {
  chronic_pain: {
    title: 'Pain Management',
    emoji: '🔥',
    color: 'from-orange-500 to-red-500',
    cannabinoid: 'High THC (15-25%) with low CBD',
    terpenes: ['Myrcene', 'Beta-Caryophyllene', 'Linalool'],
    tips: ['Indica strains preferred', 'Consider edibles for long-lasting relief', 'Start low, go slow']
  },
  anxiety: {
    title: 'Anxiety Relief',
    emoji: '🧘',
    color: 'from-blue-500 to-cyan-500',
    cannabinoid: 'High CBD (15-25%) with minimal THC',
    terpenes: ['Linalool', 'Limonene', 'Pinene'],
    tips: ['Avoid high THC strains', 'CBD-dominant products recommended', 'Consider tinctures for precise dosing']
  },
  insomnia: {
    title: 'Sleep Support',
    emoji: '🌙',
    color: 'from-purple-500 to-indigo-500',
    cannabinoid: 'Moderate THC (10-18%) + CBN if available',
    terpenes: ['Myrcene', 'Linalool', 'Terpinolene'],
    tips: ['Take 1-2 hours before bed', 'Indica strains work best', 'Avoid caffeine after 2pm']
  },
  ptsd: {
    title: 'PTSD Support',
    emoji: '💚',
    color: 'from-teal-500 to-emerald-500',
    cannabinoid: 'Balanced THC:CBD (1:1 or 2:1)',
    terpenes: ['Linalool', 'Limonene', 'Beta-Caryophyllene'],
    tips: ['Start with low THC', 'Consider daytime CBD, evening THC', 'Track symptoms and dosage']
  },
  depression: {
    title: 'Mood Elevation',
    emoji: '☀️',
    color: 'from-amber-500 to-yellow-500',
    cannabinoid: 'Moderate THC (10-15%) with CBD',
    terpenes: ['Limonene', 'Pinene', 'Terpinolene'],
    tips: ['Sativa strains for daytime', 'Avoid indica during the day', 'Combine with physical activity']
  },
  inflammation: {
    title: 'Inflammation',
    emoji: '🔴',
    color: 'from-rose-500 to-pink-500',
    cannabinoid: 'CBD-dominant with some THC',
    terpenes: ['Beta-Caryophyllene', 'Myrcene', 'Pinene'],
    tips: ['Topicals for localized inflammation', 'Consider anti-inflammatory diet', 'Consistent dosing important']
  },
  nausea: {
    title: 'Nausea Relief',
    emoji: '🤢',
    color: 'from-green-500 to-lime-500',
    cannabinoid: 'THC-dominant (CB1 activation)',
    terpenes: ['Limonene', 'Ginger terpenes'],
    tips: ['Fast-acting methods preferred', 'Small doses often work better', 'Avoid on empty stomach']
  },
  weight_management: {
    title: 'Weight & Metabolism',
    emoji: '⚖️',
    color: 'from-cyan-500 to-blue-500',
    cannabinoid: 'THCV-rich or CBD-dominant (avoid high THC)',
    terpenes: ['Limonene', 'Pinene', 'Beta-Caryophyllene'],
    tips: ['THCV suppresses appetite (opposite of THC)', 'Avoid indica strains that cause munchies', 'CBD promotes fat browning', 'Combine with Berberine for GLP-1-like effects']
  }
}

const TERPENE_DATABASE: Record<string, {
  emoji: string
  name: string
  effect: string
  found_in: string
  benefits: string[]
}> = {
  'Myrcene': {
    emoji: '🥭',
    name: 'Myrcene',
    effect: 'Relaxation & Sedation',
    found_in: 'Mangoes, Hops, Lemongrass',
    benefits: ['Enhances THC absorption', 'Muscle relaxation', 'Anti-inflammatory']
  },
  'Limonene': {
    emoji: '🍋',
    name: 'Limonene',
    effect: 'Mood & Energy',
    found_in: 'Citrus fruits, Juniper',
    benefits: ['Elevates mood', 'Reduces stress', 'Anti-anxiety']
  },
  'Linalool': {
    emoji: '💜',
    name: 'Linalool',
    effect: 'Calming & Anti-Anxiety',
    found_in: 'Lavender, Birch bark',
    benefits: ['Reduces anxiety', 'Promotes sleep', 'Anti-convulsant']
  },
  'Pinene': {
    emoji: '🌲',
    name: 'Pinene',
    effect: 'Focus & Alertness',
    found_in: 'Pine needles, Rosemary',
    benefits: ['Improves memory', 'Opens airways', 'Counteracts THC fog']
  },
  'Beta-Caryophyllene': {
    emoji: '🌶️',
    name: 'Beta-Caryophyllene',
    effect: 'Pain & Inflammation',
    found_in: 'Black pepper, Cloves',
    benefits: ['Activates CB2 receptors', 'Anti-inflammatory', 'Stress relief']
  },
  'THCV': {
    emoji: '⚡',
    name: 'THCV (Tetrahydrocannabivarin)',
    effect: 'Metabolism & Energy',
    found_in: 'African sativas, Durban Poison',
    benefits: ['Suppresses appetite', 'Improves insulin sensitivity', 'Energizing without anxiety', 'May aid weight management']
  }
}

export default function BudtenderAssistant() {
  const [customer, setCustomer] = useState<any>(null)
  const [recommendations, setRecommendations] = useState<any[]>([])

  // Dynamic guidance based on selected conditions
  const activeGuidance = useMemo(() => {
    if (!customer?.conditions?.length) return []
    return customer.conditions
      .map((c: string) => CONDITION_GUIDANCE[c])
      .filter(Boolean)
  }, [customer?.conditions])

  // Get unique terpenes from active conditions
  const relevantTerpenes = useMemo(() => {
    if (!activeGuidance.length) return Object.values(TERPENE_DATABASE).slice(0, 3)
    const terpeneSet = new Set<string>()
    activeGuidance.forEach((g: any) => g.terpenes.forEach((t: string) => terpeneSet.add(t)))
    return Array.from(terpeneSet)
      .map(t => TERPENE_DATABASE[t])
      .filter(Boolean)
      .slice(0, 4)
  }, [activeGuidance])

  return (
    <div className="min-h-screen animated-bg">
      <Head>
        <title>NeuroBotanica - Budtender Education Assistant</title>
        <meta name="viewport" content="width=device-width, initial-scale=1, maximum-scale=1, user-scalable=no" />
        <meta name="description" content="AI-powered cannabis education for budtenders" />
        <link rel="icon" href="/favicon.ico" />
        <link href="https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700;800&display=swap" rel="stylesheet" />
      </Head>

      {/* Header */}
      <header className="bg-black/30 backdrop-blur-lg border-b border-white/10 sticky top-0 z-50">
        <div className="container mx-auto px-4 py-4">
          <div className="flex justify-between items-center">
            <div className="flex items-center gap-4">
              <div className="w-12 h-12 bg-gradient-to-br from-emerald-400 to-teal-500 rounded-xl flex items-center justify-center shadow-lg shadow-emerald-500/30 animate-pulse-glow">
                <svg className="w-7 h-7 text-white" fill="currentColor" viewBox="0 0 24 24">
                  <path d="M12 2C6.5 2 2 6.5 2 12s4.5 10 10 10 10-4.5 10-10S17.5 2 12 2zm-2 15l-5-5 1.41-1.41L10 14.17l7.59-7.59L19 8l-9 9z"/>
                </svg>
              </div>
              <div>
                <h1 className="text-2xl font-bold text-white text-glow-white">
                  NeuroBotanica
                </h1>
                <p className="text-emerald-300 text-sm font-medium">
                  🌿 Budtender Education Assistant
                </p>
              </div>
            </div>
            
            <div className="hidden md:flex items-center gap-3">
              <div className="stat-card">
                <div className="text-xl font-bold text-white">505+</div>
                <div className="text-xs text-emerald-300">Clinical Studies</div>
              </div>
              <div className="stat-card">
                <div className="text-xl font-bold text-white">9</div>
                <div className="text-xs text-purple-300">Adjuvants</div>
              </div>
              <div className="stat-card">
                <div className="text-xl font-bold text-white">63</div>
                <div className="text-xs text-amber-300">Compounds</div>
              </div>
            </div>
          </div>
        </div>
      </header>

      {/* Main Content */}
      <main className="container mx-auto px-4 py-6">
        <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
          
          {/* Left Column */}
          <div className="space-y-6">
            {/* Customer Search */}
            <div className="glass-card p-6">
              <div className="flex items-center gap-3 mb-5">
                <div className="w-10 h-10 bg-gradient-to-r from-emerald-500 to-teal-500 rounded-full flex items-center justify-center text-white font-bold shadow-lg">
                  1
                </div>
                <h2 className="text-xl font-bold text-white">Customer Consultation</h2>
              </div>
              <CustomerSearch onCustomerSelect={setCustomer} />
            </div>

            {/* Customer Profile */}
            {customer && (
              <div className="glass-card p-6">
                <div className="flex items-center gap-3 mb-5">
                  <div className="w-10 h-10 bg-gradient-to-r from-purple-500 to-indigo-500 rounded-full flex items-center justify-center text-white font-bold shadow-lg">
                    2
                  </div>
                  <h2 className="text-xl font-bold text-white">Customer Needs</h2>
                </div>
                <CustomerProfile customer={customer} onProfileUpdate={setCustomer} />
              </div>
            )}

            {/* Dynamic Quick Reference - Updates based on conditions */}
            <div className="glass-card p-6">
              <h3 className="text-lg font-bold text-white mb-4 flex items-center gap-2">
                <span className="text-2xl">📚</span>
                {activeGuidance.length > 0 ? 'Personalized Guidance' : 'Quick Reference Guide'}
              </h3>
              
              {activeGuidance.length > 0 ? (
                <div className="space-y-4">
                  {activeGuidance.map((guide: any, idx: number) => (
                    <div 
                      key={idx}
                      className={`bg-gradient-to-r ${guide.color} p-0.5 rounded-xl`}
                    >
                      <div className="bg-gray-900/95 rounded-[10px] p-4">
                        <div className="flex items-center gap-2 mb-3">
                          <span className="text-2xl">{guide.emoji}</span>
                          <span className="font-bold text-white text-lg">{guide.title}</span>
                        </div>
                        <div className="space-y-3">
                          <div>
                            <div className="text-white/60 text-xs uppercase font-semibold mb-1">Cannabinoid Profile</div>
                            <div className="text-white font-medium">{guide.cannabinoid}</div>
                          </div>
                          <div>
                            <div className="text-white/60 text-xs uppercase font-semibold mb-1">Key Terpenes</div>
                            <div className="flex flex-wrap gap-2">
                              {guide.terpenes.map((t: string, i: number) => (
                                <span key={i} className="bg-white/10 text-white px-2 py-1 rounded-full text-sm">
                                  {TERPENE_DATABASE[t]?.emoji || '🌿'} {t}
                                </span>
                              ))}
                            </div>
                          </div>
                          <div>
                            <div className="text-white/60 text-xs uppercase font-semibold mb-1">Tips</div>
                            <ul className="text-white/80 text-sm space-y-1">
                              {guide.tips.map((tip: string, i: number) => (
                                <li key={i} className="flex items-start gap-2">
                                  <span className="text-emerald-400">•</span>
                                  {tip}
                                </li>
                              ))}
                            </ul>
                          </div>
                        </div>
                      </div>
                    </div>
                  ))}
                </div>
              ) : (
                <div className="grid grid-cols-2 gap-3">
                  <div className="bg-gradient-to-r from-orange-500/20 to-amber-500/20 rounded-xl p-4 border border-orange-400/30">
                    <div className="text-orange-300 font-bold mb-1">🔥 Pain Relief</div>
                    <div className="text-white/80 text-sm">High THC + Myrcene</div>
                    <div className="text-orange-200/60 text-xs mt-1">Indica strains preferred</div>
                  </div>
                  <div className="bg-gradient-to-r from-blue-500/20 to-cyan-500/20 rounded-xl p-4 border border-blue-400/30">
                    <div className="text-blue-300 font-bold mb-1">🧘 Anxiety</div>
                    <div className="text-white/80 text-sm">CBD + Linalool</div>
                    <div className="text-blue-200/60 text-xs mt-1">Low THC recommended</div>
                  </div>
                  <div className="bg-gradient-to-r from-purple-500/20 to-violet-500/20 rounded-xl p-4 border border-purple-400/30">
                    <div className="text-purple-300 font-bold mb-1">🌙 Sleep</div>
                    <div className="text-white/80 text-sm">Indica + CBN + Myrcene</div>
                    <div className="text-purple-200/60 text-xs mt-1">Evening use only</div>
                  </div>
                  <div className="bg-gradient-to-r from-green-500/20 to-emerald-500/20 rounded-xl p-4 border border-green-400/30">
                    <div className="text-green-300 font-bold mb-1">⚡ Focus</div>
                    <div className="text-white/80 text-sm">Sativa + Pinene</div>
                    <div className="text-green-200/60 text-xs mt-1">Morning/day use</div>
                  </div>
                </div>
              )}
            </div>
          </div>

          {/* Right Column */}
          <div className="space-y-6">
            {customer ? (
              <div className="glass-card p-6">
                <div className="flex items-center gap-3 mb-5">
                  <div className="w-10 h-10 bg-gradient-to-r from-amber-500 to-orange-500 rounded-full flex items-center justify-center text-white font-bold shadow-lg animate-float">
                    3
                  </div>
                  <h2 className="text-xl font-bold text-white">NeuroBotanica Insights</h2>
                  <span className="ml-auto bg-emerald-500/20 text-emerald-300 px-3 py-1 rounded-full text-xs font-semibold border border-emerald-400/30">
                    ✨ Evidence-Based
                  </span>
                </div>
                <ProductRecommendations
                  customer={customer}
                  recommendations={recommendations}
                  onRecommendationsUpdate={setRecommendations}
                />
              </div>
            ) : (
              <div className="glass-card p-8 text-center">
                <div className="w-24 h-24 mx-auto mb-6 bg-gradient-to-br from-emerald-400 to-teal-500 rounded-full flex items-center justify-center shadow-lg shadow-emerald-500/30 animate-float">
                  <svg className="w-12 h-12 text-white" fill="currentColor" viewBox="0 0 24 24">
                    <path d="M17.66 11.2C17.43 10.9 17.15 10.64 16.89 10.38C16.22 9.78 15.46 9.35 14.82 8.72C13.33 7.26 13 4.85 13.95 3C13 3.23 12.17 3.75 11.46 4.32C8.87 6.4 7.85 10.07 9.07 13.22C9.11 13.32 9.15 13.42 9.15 13.55C9.15 13.77 9 13.97 8.8 14.05C8.57 14.15 8.33 14.09 8.14 13.93C8.08 13.88 8.04 13.83 8 13.76C6.87 12.33 6.69 10.28 7.45 8.64C5.78 10 4.87 12.3 5 14.47C5.06 14.97 5.12 15.47 5.29 15.97C5.43 16.57 5.7 17.17 6 17.7C7.08 19.43 8.95 20.67 10.96 20.92C13.1 21.19 15.39 20.8 17.03 19.32C18.86 17.66 19.5 15 18.56 12.72L18.43 12.46C18.22 12 17.66 11.2 17.66 11.2Z"/>
                  </svg>
                </div>
                <h3 className="text-2xl font-bold text-white mb-3 text-glow-white">
                  Ready to Educate
                </h3>
                <p className="text-white/70 mb-6 max-w-md mx-auto">
                  Start a consultation to receive personalized therapeutic insights backed by <span className="text-emerald-400 font-semibold">505+ clinical studies</span>.
                </p>
                
                <div className="grid grid-cols-3 gap-4 mb-6">
                  <div className="bg-white/5 rounded-xl p-4 border border-white/10">
                    <div className="text-3xl mb-2">🧬</div>
                    <div className="text-white text-sm font-medium">Terpene Analysis</div>
                    <div className="text-white/50 text-xs">Personalized profiles</div>
                  </div>
                  <div className="bg-white/5 rounded-xl p-4 border border-white/10">
                    <div className="text-3xl mb-2">💊</div>
                    <div className="text-white text-sm font-medium">Adjuvant Synergy</div>
                    <div className="text-white/50 text-xs">Enhanced effects</div>
                  </div>
                  <div className="bg-white/5 rounded-xl p-4 border border-white/10">
                    <div className="text-3xl mb-2">📊</div>
                    <div className="text-white text-sm font-medium">Clinical Evidence</div>
                    <div className="text-white/50 text-xs">Research-backed</div>
                  </div>
                </div>
              </div>
            )}

            {/* Dynamic Terpene Education Panel */}
            <div className="glass-card p-6">
              <h3 className="text-lg font-bold text-white mb-4 flex items-center gap-2">
                <span className="text-2xl">🧬</span>
                {activeGuidance.length > 0 ? 'Recommended Terpenes' : 'Terpene Education'}
              </h3>
              <div className="space-y-3">
                {relevantTerpenes.map((terpene, idx) => (
                  <div key={idx} className="flex items-start gap-3 bg-white/5 rounded-xl p-4 border border-white/10">
                    <div className="w-12 h-12 bg-gradient-to-br from-emerald-500/30 to-teal-500/30 rounded-full flex items-center justify-center flex-shrink-0">
                      <span className="text-2xl">{terpene.emoji}</span>
                    </div>
                    <div className="flex-1">
                      <div className="flex items-center justify-between mb-1">
                        <div className="text-white font-semibold">{terpene.name}</div>
                        <span className="text-emerald-400 text-xs bg-emerald-500/20 px-2 py-0.5 rounded-full">
                          {terpene.effect}
                        </span>
                      </div>
                      <div className="text-white/50 text-xs mb-2">Found in: {terpene.found_in}</div>
                      <div className="flex flex-wrap gap-1">
                        {terpene.benefits.map((b, i) => (
                          <span key={i} className="text-xs bg-white/10 text-white/70 px-2 py-0.5 rounded-full">
                            {b}
                          </span>
                        ))}
                      </div>
                    </div>
                  </div>
                ))}
              </div>
            </div>
          </div>
        </div>
      </main>

      {/* Footer */}
      <footer className="bg-black/40 backdrop-blur-lg border-t border-white/10 py-4 mt-8">
        <div className="container mx-auto px-4 text-center">
          <p className="text-white/70 text-sm">
            <span className="text-emerald-400 font-semibold">NeuroBotanica</span> • Evidence-Based Cannabis Education
          </p>
          <p className="text-white/40 text-xs mt-1">
            Powered by Cloak and Quill Research 501(c)(3) | Henderson, Nevada
          </p>
        </div>
      </footer>
    </div>
  )
}
