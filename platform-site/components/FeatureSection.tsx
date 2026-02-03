import { ArrowRight, Database, Shield, Activity, Dna, Layers, Zap, Brain, Heart } from 'lucide-react';

const kingdoms = [
  {
    name: 'Cannabis (The Foundation)',
    focus: 'Dimeric cannabinoids, terpene synergy, and 16 neuropsychiatric conditions.',
    status: 'Live · 505+ indexed studies',
    icon: Activity
  },
  {
    name: 'Marine Polysaccharides',
    focus: 'Fucoidan and carrageenan for immunomodulation and oncology support.',
    status: 'Coming Q2 2026',
    icon: Zap
  },
  {
    name: 'Fungal Therapeutics',
    focus: 'Beta-glucans and hericenones for neuroprotection and immune enhancement.',
    status: 'In development',
    icon: Database
  },
  {
    name: 'Plant Polysaccharides',
    focus: 'Galactomannans for glycemic control and gut health.',
    status: 'In development',
    icon: Layers
  }
];

const infrastructures = [
  {
    title: 'Evidence-Based',
    description: 'ChemPath™ + BioPath™ execute TK validation, CoA checks, and cultural provenance scoring.',
    color: 'bg-blue-50 text-blue-700'
  },
  {
    title: 'AI-Powered',
    description: 'GenomePath™ inference blends Morgan fingerprints with multilingual embeddings at the edge.',
    color: 'bg-red-50 text-red-700'
  },
  {
    title: 'Comprehensive',
    description: 'RegPath™ unifies terpene, cannabinoid, and polysaccharide dimers with regulatory telemetry.',
    color: 'bg-purple-50 text-purple-700'
  },
  {
    title: 'Cloak & Quill Research',
    description: 'Cloak & Quill Research — a 501(c)(3) nonprofit stewarding NeuroBotanica for public benefit.',
    color: 'bg-teal-50 text-teal-700'
  }
];

const fingerprintBits = ['3A', '7F', '12', 'C8', '0E', '5D', '1B', '9C', '2F', 'AA', '44', '6B'];

const tkEmbeddings = [
  {
    label: 'TRPV1 Downshift',
    detail: '“Cooling breath” intents spanning Mapuche + Yawanawá pharmacopeia'
  },
  {
    label: 'CB2 Analgesic Arc',
    detail: 'Philippines + Yucatán chronic pain rituals clustered via multilingual embeddings'
  },
  {
    label: 'Neuroimmune Shield',
    detail: 'Marine polysaccharide rites for post-operative recovery with beta-glucan co-signals'
  }
];

function MorganTKVisualization() {
  return (
    <div className="relative h-full w-full p-6 sm:p-10">
      <div className="absolute inset-0 bg-gradient-to-br from-emerald-500/10 via-slate-900/60 to-indigo-900/70" aria-hidden />
      <svg className="absolute inset-0 h-full w-full opacity-40" viewBox="0 0 600 400" aria-hidden>
        <defs>
          <linearGradient id="fusion" x1="0%" y1="0%" x2="100%" y2="0%">
            <stop offset="0%" stopColor="#6EE7B7" />
            <stop offset="100%" stopColor="#A5B4FC" />
          </linearGradient>
        </defs>
        <path d="M80 80 C 200 160, 280 120, 360 180" stroke="url(#fusion)" strokeWidth="3" fill="none" strokeLinecap="round" />
        <path d="M100 220 C 220 220, 260 280, 380 260" stroke="url(#fusion)" strokeWidth="2" fill="none" strokeLinecap="round" strokeDasharray="8 10" />
        <path d="M60 320 C 220 280, 300 360, 420 320" stroke="url(#fusion)" strokeWidth="2" fill="none" strokeLinecap="round" opacity="0.6" />
      </svg>
      <div className="relative grid gap-6 md:grid-cols-[1fr_auto_1fr] items-center">
        <div className="grid grid-cols-3 sm:grid-cols-4 gap-3">
          {fingerprintBits.map((bit) => (
            <div
              key={bit}
              className="h-12 w-12 rounded-2xl border border-emerald-300/40 bg-emerald-400/10 backdrop-blur flex items-center justify-center text-xs font-mono text-emerald-100 shadow-[0_10px_25px_rgba(0,0,0,0.35)]"
            >
              {bit}
            </div>
          ))}
        </div>
        <div className="hidden md:flex flex-col items-center gap-2 text-white/80">
          <div className="w-20 h-20 rounded-full border border-white/40 flex items-center justify-center bg-white/10 backdrop-blur">
            <span className="text-[10px] font-mono tracking-[0.4em] text-white/70">ECFP4</span>
          </div>
          <span className="text-[10px] uppercase tracking-[0.5em] text-white/40">Fusion</span>
        </div>
        <div className="space-y-3">
          {tkEmbeddings.map((token) => (
            <div key={token.label} className="rounded-2xl border border-indigo-300/40 bg-white/5 px-4 py-3 text-sm text-slate-100 shadow-[0_12px_30px_rgba(0,0,0,0.45)]">
              <p className="font-semibold text-emerald-200">{token.label}</p>
              <p className="text-xs text-white/70">{token.detail}</p>
            </div>
          ))}
        </div>
      </div>
    </div>
  );
}

const specializedLayers = [
  {
    title: 'PsychePath™',
    description: 'Psychological therapeutics addressing mental health care disparities.',
    icon: Brain,
    accent: 'text-indigo-600'
  },
  {
    title: 'SereniPath™',
    description: 'Cross-cultural wellness and anxiety frameworks with equity guardrails.',
    icon: Heart,
    accent: 'text-rose-500'
  }
];

const roadmap = [
  {
    quarter: 'Q1 2026 (Now)',
    bullets: ['Nevada dispensary pilot live'],
  },
  {
    quarter: 'Q2 2026 – Specialized Therapeutics',
    bullets: ['PsychePath™ + SereniPath™ integrations', 'Enterprise API beta for NeuroBotanica Discovery endpoints'],
  },
  {
    quarter: 'Q3 2026',
    bullets: ['Cross-kingdom synergy engines live', 'Geographic expansion into California'],
  }
];

export default function FeatureSection() {
  return (
    <div className="bg-slate-50 rounded-[32px] border border-white/60 shadow-parchment py-24 sm:py-32 px-6 lg:px-12">
      <div className="mx-auto max-w-6xl">
        {/* Section 1: Four Kingdoms */}
        <div className="mx-auto max-w-2xl text-center mb-16">
          <h2 className="text-base font-semibold leading-7 text-emerald-600">Cross-Kingdom Intelligence</h2>
          <p className="mt-2 text-3xl font-bold tracking-tight text-gray-900 sm:text-4xl">Synergy Beyond Boundaries</p>
          <p className="mt-6 text-lg leading-8 text-gray-600">
            True therapeutic breakthroughs happen at the intersection of kingdoms. NeuroBotanica&apos;s inference engine predicts high-value dimeric combinations across four distinct domains.
          </p>
        </div>

        <dl className="grid gap-10 lg:grid-cols-4">
          {kingdoms.map((feature) => (
            <div key={feature.name} className="flex flex-col rounded-3xl bg-white p-6 shadow-parchment border border-sand/60">
              <dt className="flex items-center gap-3 text-base font-semibold leading-7 text-gray-900">
                <feature.icon className="h-5 w-5 flex-none text-emerald-600" aria-hidden="true" />
                <span>{feature.name}</span>
              </dt>
              <dd className="mt-4 flex flex-auto flex-col text-base leading-7 text-gray-600">
                <p className="flex-auto">{feature.focus}</p>
                <span className="mt-4 inline-flex items-center justify-start rounded-full bg-slate-100 px-3 py-1 text-xs font-semibold text-slate-600">
                  {feature.status}
                </span>
              </dd>
            </div>
          ))}
        </dl>

        {/* Section 2: Enterprise Infrastructure */}
        <div className="mx-auto max-w-2xl text-center mt-28 mb-12">
          <h2 className="text-base font-semibold leading-7 text-blue-600">Powered by OmniPath™</h2>
          <p className="mt-2 text-3xl font-bold tracking-tight text-gray-900 sm:text-4xl">The Infrastructure of Compliance</p>
          <p className="mt-6 text-lg leading-8 text-gray-600">
            Discovery is only the first step. Powered by OmniPath™, our proprietary trade-secret engines automate the journey from molecular identification to regulatory submission.
          </p>
        </div>

        <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-4 gap-6">
          {infrastructures.map((item) => (
            <div key={item.title} className={`p-8 rounded-2xl text-center ${item.color}`}>
              <h3 className="text-xl font-bold mb-2">{item.title}</h3>
              <p className="text-sm opacity-90">{item.description}</p>
            </div>
          ))}
        </div>

        {/* Section 3: AI */}
        <div className="mt-28 overflow-hidden bg-gray-900 rounded-3xl shadow-2xl lg:grid lg:grid-cols-2 lg:gap-4">
          <div className="px-6 pb-12 pt-10 sm:px-16 sm:pt-16 lg:py-16 lg:pr-0 xl:py-20 xl:px-20">
            <div className="lg:self-center">
              <h2 className="text-3xl font-bold tracking-tight text-white sm:text-4xl">
                <span className="block text-emerald-400 mb-2">GenomePath™ AI</span>
                Bridging Ancient Wisdom & Modern Genomics
              </h2>
              <p className="mt-4 text-lg leading-6 text-gray-300">
                The vocabulary of Traditional Knowledge and the language of modern genomics have never spoken the same dialect—until now.
              </p>
              <ul className="mt-8 space-y-4 text-gray-400">
                <li className="flex items-center"><Dna className="w-5 h-5 mr-3 text-emerald-500" /> 98.11% accuracy linking TK intents to genomic targets (e.g., TRPV1 downregulation).</li>
                <li className="flex items-center"><Shield className="w-5 h-5 mr-3 text-emerald-500" /> Bias-aware and built to route 70% of commercial value to Indigenous stakeholders.</li>
                <li className="flex items-center"><ArrowRight className="w-5 h-5 mr-3 text-emerald-500" /> Hybrid AI combining Morgan ECFP4 fingerprints with multilingual text embeddings.</li>
              </ul>
            </div>
          </div>
          <MorganTKVisualization />
        </div>

        {/* Section 4: Specialized Therapeutics */}
        <div className="mx-auto mt-28 max-w-2xl text-center mb-10">
          <h2 className="text-base font-semibold leading-7 text-indigo-600">Coming Q2 2026</h2>
          <p className="mt-2 text-3xl font-bold tracking-tight text-gray-900">Specialized Therapeutic Layers</p>
          <p className="mt-4 text-lg leading-7 text-gray-600">Announcing PsychePath™ and SereniPath™ to extend NeuroBotanica&apos;s operating system into mental wellness and cross-cultural care.</p>
        </div>
        <div className="grid grid-cols-1 md:grid-cols-2 gap-8 max-w-4xl mx-auto">
          {specializedLayers.map((layer) => (
            <div key={layer.title} className="flex flex-col items-center text-center p-8 bg-white rounded-3xl shadow-parchment border border-slate-200">
              <layer.icon className={`w-10 h-10 mb-4 ${layer.accent}`} />
              <h3 className="text-xl font-bold text-slate-900">{layer.title}</h3>
              <p className="mt-2 text-slate-600">{layer.description}</p>
            </div>
          ))}
        </div>

        {/* Section 5: Roadmap */}
        <div className="mt-24">
          <h3 className="text-lg font-semibold text-gray-900 mb-6 text-center">2026 Platform Roadmap</h3>
          <div className="grid gap-6 md:grid-cols-3">
            {roadmap.map((phase) => (
              <div key={phase.quarter} className="rounded-2xl border border-sand/60 bg-white p-6 shadow-parchment">
                <p className="text-sm font-mono uppercase tracking-[0.3em] text-indigo/60">{phase.quarter}</p>
                <ul className="mt-4 space-y-2 text-sm text-gray-700">
                  {phase.bullets.map((item) => (
                    <li key={item} className="flex items-start gap-2">
                      <ArrowRight className="w-4 h-4 text-emerald-500 mt-0.5" />
                      <span>{item}</span>
                    </li>
                  ))}
                </ul>
              </div>
            ))}
          </div>
        </div>
      </div>
    </div>
  );
}
