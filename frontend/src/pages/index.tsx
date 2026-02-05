import { useState, useMemo, useEffect, useCallback } from 'react'
import Head from 'next/head'
import { useRouter } from 'next/router'
import Link from 'next/link'
import { CustomerSearch } from '@/components/CustomerSearch'
import { ProductRecommendations } from '@/components/ProductRecommendations'
import type { Recommendation } from '@/components/ProductRecommendations'
import { CustomerProfile } from '@/components/CustomerProfile'
import type { CustomerProfileData } from '@/types/customer'
import { conditionNames } from '@/utils/conditions'
import { SANDBOX_CUSTOMER } from '@/content/sandbox/customer'

// Condition-specific guidance data
type ConditionGuidance = {
  title: string
  emoji: string
  color: string
  cannabinoid: string
  terpenes: string[]
  tips: string[]
}

const CONDITION_GUIDANCE: Record<string, ConditionGuidance> = {
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

type TerpeneInfo = {
  emoji: string
  name: string
  effect: string
  found_in: string
  benefits: string[]
}

const TERPENE_DATABASE: Record<string, TerpeneInfo> = {
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
    emoji: '🌸',
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

type ConditionKey = keyof typeof CONDITION_GUIDANCE
type TerpeneKey = keyof typeof TERPENE_DATABASE

const NAV_LINKS = [
  { label: 'Consultation', href: '#consultation' },
  { label: 'Customer Needs', href: '#customer-needs' },
  { label: 'Insights', href: '#insights' },
  { label: 'Terpenes', href: '#terpene-education' },
  { label: 'Orientation', href: '/orientation' }
]

const HEADER_METRICS = [
  { label: 'Clinical Studies', value: '505+', accent: 'text-emerald-300' },
  { label: 'Adjuvants', value: '15', accent: 'text-purple-300' },
  { label: 'Compounds', value: '63', accent: 'text-amber-300' }
]

const SANDBOX_FLAG = 'nb_sandbox_mode'

function OrientationPrompt() {
  const [isComplete, setIsComplete] = useState<boolean | null>(null);
  
  // Check completion status (client-side only)
  useEffect(() => {
    const completed = localStorage.getItem('orientation_completed') === 'true';
    setIsComplete(completed);
  }, []);
  
  // Don't show anything during SSR
  if (isComplete === null) return null;
  
  // Don't show if orientation is complete
  if (isComplete) return null;
  
  return (
    <div className="bg-gradient-to-r from-purple-600 to-blue-600 border-2 border-purple-400 rounded-lg p-6 mb-6">
      <div className="flex items-start gap-4">
        <div className="flex-shrink-0">
          <svg className="w-12 h-12 text-white" fill="none" stroke="currentColor" viewBox="0 0 24 24">
            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M12 6.253v13m0-13C10.832 5.477 9.246 5 7.5 5S4.168 5.477 3 6.253v13C4.168 18.477 5.754 18 7.5 18s3.332.477 4.5 1.253m0-13C13.168 5.477 14.754 5 16.5 5c1.747 0 3.332.477 4.5 1.253v13C19.832 18.477 18.247 18 16.5 18c-1.746 0-3.332.477-4.5 1.253" />
          </svg>
        </div>
        
        <div className="flex-1">
          <h3 className="text-2xl font-bold text-white mb-2">
            Complete Orientation First
          </h3>
          <p className="text-white/90 mb-4">
            This professional certification tool requires completion of a 15-minute orientation covering evidence-based recommendations, biomarker interpretation, and Nevada compliance requirements.
          </p>
          
          <Link 
            href="/orientation"
            className="inline-flex items-center gap-2 bg-white text-purple-600 font-semibold px-6 py-3 rounded-lg hover:bg-purple-50 transition-colors"
          >
            <svg className="w-5 h-5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M14 5l7 7m0 0l-7 7m7-7H3" />
            </svg>
            Start Orientation
          </Link>
        </div>
      </div>
    </div>
  );
}

export default function BudtenderAssistant() {
  // ✅ MOVE ALL useState HOOKS TO THE TOP (before any returns)
  const [isHydrated, setIsHydrated] = useState(false)
  const [customer, setCustomer] = useState<CustomerProfileData | null>(null)
  const [recommendations, setRecommendations] = useState<Recommendation[]>([])
  const [isSandboxMode, setIsSandboxMode] = useState(false)
  const [isMobileNavOpen, setIsMobileNavOpen] = useState(false)
  const [customerSearchKey, setCustomerSearchKey] = useState(0)
  const router = useRouter()

  const conditionList = useMemo(() => conditionNames(customer?.conditions), [customer?.conditions])

  // Mark as hydrated after first render
  useEffect(() => {
    setIsHydrated(true)
  }, [])

  // Check orientation completion status and redirect if needed (only after hydration)
  useEffect(() => {
    if (!isHydrated) return

    // NO automatic redirects - just check status
    // OrientationPrompt component will handle the UI
  }, [isHydrated])

  // Dynamic guidance based on selected conditions
  const activeGuidance = useMemo<ConditionGuidance[]>(() => {
    if (!conditionList.length) return []
    return conditionList
      .map((condition) => CONDITION_GUIDANCE[condition as ConditionKey])
      .filter((guide): guide is ConditionGuidance => Boolean(guide))
  }, [conditionList])

  // Get unique terpenes from active conditions
  const relevantTerpenes = useMemo<TerpeneInfo[]>(() => {
    if (!activeGuidance.length) return Object.values(TERPENE_DATABASE).slice(0, 3)
    const terpeneSet = new Set<string>()
    activeGuidance.forEach((guidance) => guidance.terpenes.forEach((terpene) => terpeneSet.add(terpene)))
    return Array.from(terpeneSet)
      .map((terpene) => TERPENE_DATABASE[terpene as TerpeneKey])
      .filter((entry): entry is TerpeneInfo => Boolean(entry))
      .slice(0, 4)
  }, [activeGuidance])

  const enterSandbox = useCallback(() => {
    const { ...sandboxCustomerWithoutFlag } = SANDBOX_CUSTOMER
    setCustomer(sandboxCustomerWithoutFlag)
    setRecommendations([])
    setIsSandboxMode(true)
    if (typeof window !== 'undefined') {
      window.localStorage.setItem(SANDBOX_FLAG, '1')
    }
  }, [])

  const exitSandbox = useCallback(() => {
    setCustomer(null)
    setRecommendations([])
    setIsSandboxMode(false)
    if (typeof window !== 'undefined') {
      window.localStorage.removeItem(SANDBOX_FLAG)
    }
  }, [])

  const handleNewClient = useCallback(() => {
    setCustomer(null)
    setRecommendations([])
    setCustomerSearchKey((prev) => prev + 1)

    if (isSandboxMode) {
      enterSandbox()
      return
    }

    const newCustomer: CustomerProfileData = {
      customer_id: `temp_${Date.now()}`,
      first_name: '',
      last_name: '',
      age: undefined,
      gender: '',
      weight: undefined,
      conditions: [],
      experience_level: 'beginner',
      notes: '',
      biomarkers: {},
      isNew: true,
      isSandbox: false
    }

    setCustomer(newCustomer)
  }, [enterSandbox, isSandboxMode])

  useEffect(() => {
    if (!router.isReady) return
    const sandboxQueryParam = router.query.sandbox
    const querySandbox = Array.isArray(sandboxQueryParam)
      ? sandboxQueryParam.some((value) => value === '1' || value === 'true')
      : sandboxQueryParam === '1' || sandboxQueryParam === 'true'
    if (querySandbox && !isSandboxMode) {
      enterSandbox()
      if (querySandbox) {
        router.replace('/', undefined, { shallow: true })
      }
    }
  }, [router, enterSandbox, isSandboxMode])

  useEffect(() => {
    if (typeof document === 'undefined') return
    if (isMobileNavOpen) {
      document.body.style.setProperty('overflow', 'hidden')
    } else {
      document.body.style.removeProperty('overflow')
    }
    return () => {
      document.body.style.removeProperty('overflow')
    }
  }, [isMobileNavOpen])

  // Show loading state while checking orientation or during SSR
  if (!isHydrated) {
    return (
      <div className="min-h-screen bg-gradient-to-br from-slate-900 via-slate-800 to-slate-900 flex items-center justify-center">
        <div className="text-center space-y-4">
          <div className="w-8 h-8 border-4 border-emerald-500 border-t-transparent rounded-full animate-spin mx-auto"></div>
          <p className="text-white/70">
            {!isHydrated ? 'Loading NeuroBotanica...' : 'Checking certification status...'}
          </p>
        </div>
      </div>
    );
  }

  const closeMobileNav = () => setIsMobileNavOpen(false)

  return (
    <div className="min-h-screen animated-bg">
      <Head>
        <title>NeuroBotanica - Budtender Education Assistant</title>
        <meta name="viewport" content="width=device-width, initial-scale=1, maximum-scale=1, user-scalable=no" />
        <meta name="description" content="AI-powered cannabis education for budtenders" />
        <link rel="icon" href="/favicon.ico" />
      </Head>

      <a href="#main-content" className="skip-link">
        Skip to main content
      </a>

      {/* ADD THIS LINE - Orientation prompt at the top */}
      <OrientationPrompt />

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
            <div className="flex items-center gap-3">
              <nav id="primary-navigation" aria-label="Primary" className="hidden md:flex items-center gap-2">
                {NAV_LINKS.map((link) => {
                  const className = 'text-white/80 text-sm font-semibold px-3 py-2 rounded-full hover:bg-white/10 focus-visible:bg-white/10 transition'
                  const key = link.href
                  if (link.href.startsWith('/')) {
                    return (
                      <Link key={key} href={link.href} className={className}>
                        {link.label}
                      </Link>
                    )
                  }
                  return (
                    <a key={key} href={link.href} className={className}>
                      {link.label}
                    </a>
                  )
                })}
              </nav>
              <button
                type="button"
                className="md:hidden inline-flex items-center justify-center rounded-xl bg-white/10 px-3 py-2 text-white shadow-lg focus-visible:ring-4 focus-visible:ring-emerald-500/50"
                onClick={() => setIsMobileNavOpen(true)}
                aria-label="Open navigation menu"
                aria-controls="mobile-navigation"
                aria-expanded={isMobileNavOpen}
              >
                <svg className="w-6 h-6" fill="none" stroke="currentColor" strokeWidth="1.5" viewBox="0 0 24 24" aria-hidden="true">
                  <path strokeLinecap="round" strokeLinejoin="round" d="M3.75 7.5h16.5M3.75 12h16.5M3.75 16.5h10.5" />
                </svg>
              </button>
            </div>
          </div>

          <div className="hidden md:flex items-center gap-3" aria-label="Evidence metrics">
            {HEADER_METRICS.map((metric) => (
              <div key={metric.label} className="stat-card" role="status">
                <div className="text-xl font-bold text-white">{metric.value}</div>
                <div className={`text-xs ${metric.accent}`}>{metric.label}</div>
              </div>
            ))}
          </div>
        </div>
      </header>

      {isMobileNavOpen && (
        <div className="mobile-nav-backdrop" role="dialog" aria-modal="true" aria-label="Navigation menu">
          <div className="mobile-nav-panel" id="mobile-navigation">
            <div className="flex items-center justify-between">
              <div>
                <p className="text-white/60 text-xs uppercase tracking-wide">Navigate</p>
                <p className="text-white font-semibold">Budtender Assistant</p>
              </div>
              <button
                type="button"
                className="text-white/70 hover:text-white"
                onClick={closeMobileNav}
                aria-label="Close navigation menu"
              >
                <svg className="w-6 h-6" fill="none" stroke="currentColor" strokeWidth="1.5" viewBox="0 0 24 24">
                  <path strokeLinecap="round" strokeLinejoin="round" d="M6 18L18 6M6 6l12 12" />
                </svg>
              </button>
            </div>
            <nav className="space-y-2" aria-label="Mobile primary">
              {NAV_LINKS.map((link) => {
                const className = 'block w-full bg-white/5 border border-white/10 rounded-xl px-4 py-3 text-white/90 font-semibold hover:bg-white/10'
                const clickHandler = () => closeMobileNav()
                if (link.href.startsWith('/')) {
                  return (
                    <Link key={link.href} href={link.href} className={className} onClick={clickHandler}>
                      {link.label}
                    </Link>
                  )
                }
                return (
                  <a key={link.href} href={link.href} className={className} onClick={clickHandler}>
                    {link.label}
                  </a>
                )
              })}
            </nav>
            <div className="space-y-2" aria-label="Evidence metrics">
              {HEADER_METRICS.map((metric) => (
                <div key={metric.label} className="bg-white/5 rounded-xl p-3 border border-white/10">
                  <div className="text-white text-lg font-bold">{metric.value}</div>
                  <div className={`text-xs ${metric.accent}`}>{metric.label}</div>
                </div>
              ))}
            </div>
            <p className="text-white/60 text-xs leading-relaxed">
              NeuroBotanica delivers sub-2s evidence lookups with <span className="text-emerald-300 font-semibold">85% cache hit rate</span> across Cloudflare edge.
            </p>
          </div>
        </div>
      )}

      {/* Main Content */}
      <main id="main-content" className="container mx-auto px-4 py-6">
        <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
          
            {/* Left Column */}
            <div className="space-y-6">
              <section className="glass-card p-4" aria-live="polite">
                {isSandboxMode ? (
                  <div className="flex flex-col gap-3">
                    <div>
                      <p className="text-sm uppercase tracking-wide text-emerald-300 font-semibold">Sandbox practice active</p>
                      <p className="text-white text-lg font-bold">You are in guided training mode.</p>
                      <p className="text-white/70 text-sm">Any edits stay on this tablet only. Use this mode to rehearse without touching production data.</p>
                    </div>
                    <div className="flex gap-3 flex-wrap">
                      <button
                        onClick={handleNewClient}
                        className="flex items-center gap-2 bg-purple-500 hover:bg-purple-600 text-white font-semibold px-4 py-2 rounded-lg transition-colors"
                      >
                        <svg className="w-5 h-5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                          <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M12 4v16m8-8H4" />
                        </svg>
                        New Client
                      </button>
                      <button type="button" className="btn-secondary" onClick={exitSandbox}>
                        Exit Sandbox
                      </button>
                      <Link href="/orientation" className="text-sm text-white/80 underline-offset-4 hover:underline">
                        Review lessons
                      </Link>
                    </div>
                  </div>
                ) : (
                  <div className="flex flex-col gap-3">
                    <div>
                      <p className="text-sm uppercase tracking-wide text-white/60">Need practice?</p>
                      <p className="text-white text-lg font-bold">Practice Mode</p>
                      <p className="text-white/70 text-sm">Try the recommendation tool with example patients (no real data).</p>
                    </div>
                    <div className="flex gap-3 flex-wrap">
                      <button
                        onClick={handleNewClient}
                        className="flex items-center gap-2 bg-purple-500 hover:bg-purple-600 text-white font-semibold px-4 py-2 rounded-lg transition-colors"
                      >
                        <svg className="w-5 h-5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                          <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M12 4v16m8-8H4" />
                        </svg>
                        New Client
                      </button>
                      <button type="button" className="btn-primary" onClick={enterSandbox}>
                        Start Practice Session
                      </button>
                    </div>
                  </div>
                )}
              </section>
 
              {/* Customer Search */}
              <section id="consultation" aria-labelledby="consultation-heading" className="glass-card p-6">
                <div className="flex items-center gap-3 mb-5">
                  <div className="w-10 h-10 bg-gradient-to-r from-emerald-500 to-teal-500 rounded-full flex items-center justify-center text-white font-bold shadow-lg">
                    1
                  </div>
                  <h2 id="consultation-heading" className="text-xl font-bold text-white">Customer Consultation</h2>
                </div>
                <CustomerSearch key={customerSearchKey} onCustomerSelect={setCustomer} isSandboxMode={isSandboxMode} />
              </section>

              {/* Customer Profile */}
              {customer && (
                <section id="customer-needs" aria-labelledby="customer-needs-heading" className="glass-card p-6">
                  <div className="flex items-center gap-3 mb-5">
                    <div className="w-10 h-10 bg-gradient-to-r from-purple-500 to-indigo-500 rounded-full flex items-center justify-center text-white font-bold shadow-lg">
                      2
                    </div>
                    <h2 id="customer-needs-heading" className="text-xl font-bold text-white">Customer Needs</h2>
                  </div>
                  <CustomerProfile customer={customer} onProfileUpdate={setCustomer} />
                </section>
              )}

              {/* Dynamic Quick Reference - Updates based on conditions */}
              <section className="glass-card p-6" aria-live="polite">
                <h3 className="text-lg font-bold text-white mb-4 flex items-center gap-2">
                  <span className="text-2xl">📚</span>
                  {activeGuidance.length > 0 ? 'Personalized Guidance' : 'Quick Reference Guide'}
                </h3>
                
                {activeGuidance.length > 0 ? (
                  <div className="space-y-4">
                    {activeGuidance.map((guide, idx) => (
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
                                {guide.terpenes.map((terpene, terpeneIndex) => (
                                  <span key={terpeneIndex} className="bg-white/10 text-white px-2 py-1 rounded-full text-sm">
                                    {TERPENE_DATABASE[terpene as TerpeneKey]?.emoji || '🌿'} {terpene}
                                  </span>
                                ))}
                              </div>
                            </div>
                            <div>
                              <div className="text-white/60 text-xs uppercase font-semibold mb-1">Tips</div>
                              <ul className="text-white/80 text-sm space-y-1">
                                {guide.tips.map((tip, tipIndex) => (
                                  <li key={tipIndex} className="flex items-start gap-2">
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
              </section>
            </div>

            {/* Right Column */}
            <div className="space-y-6">
              {customer ? (
                <section id="insights" aria-labelledby="insights-heading" className="glass-card p-6">
                  <div className="flex items-center gap-3 mb-5">
                    <div className="w-10 h-10 bg-gradient-to-r from-amber-500 to-orange-500 rounded-full flex items-center justify-center text-white font-bold shadow-lg animate-float">
                      3
                    </div>
                    <h2 id="insights-heading" className="text-xl font-bold text-white">NeuroBotanica Insights</h2>
                    <span className="ml-auto bg-emerald-500/20 text-emerald-300 px-3 py-1 rounded-full text-xs font-semibold border border-emerald-400/30">
                      ✨ Evidence-Based
                    </span>
                  </div>
                  <ProductRecommendations
                    customer={customer}
                    recommendations={recommendations}
                    onRecommendationsUpdate={setRecommendations}
                  />
                </section>
              ) : (
                <section className="glass-card p-8 text-center" aria-live="polite">
                  <div className="w-24 h-24 mx-auto mb-6 bg-gradient-to-br from-emerald-400 to-teal-500 rounded-full flex items-center justify-center shadow-lg shadow-emerald-500/30 animate-float">
                    <svg className="w-12 h-12 text-white" fill="currentColor" viewBox="0 0 24 24" role="img" aria-label="Cannabis leaf icon">
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
                </section>
              )}

              {/* Dynamic Terpene Education Panel */}
              <section id="terpene-education" aria-labelledby="terpene-heading" className="glass-card p-6">
                <h3 id="terpene-heading" className="text-lg font-bold text-white mb-4 flex items-center gap-2">
                  <span className="text-2xl">🧬</span>
                  {activeGuidance.length > 0 ? 'Recommended Terpenes' : 'Terpene Education'}
                </h3>
                <div className="space-y-3">
                  {relevantTerpenes.map((terpene, idx) => (
                    <div key={idx} className="flex items-start gap-3 bg-white/5 rounded-xl p-4 border border-white/10">
                      <div className="w-12 h-12 bg-gradient-to-br from-emerald-500/30 to-teal-500/30 rounded-full flex items-center justify-center flex-shrink-0">
                        <span className="text-2xl" aria-hidden="true">{terpene.emoji}</span>
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
              </section>
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
