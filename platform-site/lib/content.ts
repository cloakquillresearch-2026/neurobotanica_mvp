// Landing page content primitives
export const landingContent = {
  hero: {
    title: 'NeuroBotanica',
    subtitle: 'AI-powered botanical therapeutic analysis',
    cta: 'Get Started'
  },
  features: [
    {
      title: 'Terpene Analysis',
      description: 'Advanced analysis of botanical terpene profiles'
    },
    {
      title: 'Therapeutic Insights',
      description: 'Evidence-based therapeutic recommendations'
    },
    {
      title: 'Research-Backed',
      description: 'Powered by peer-reviewed scientific research'
    }
  ]
};

export type HeroMetric = {
  value: string;
  label: string;
  detail: string;
};

export const heroMetrics: HeroMetric[] = [
  { value: '505+', label: 'Studies Parsed', detail: 'Peer-reviewed + ethnobotanical corpus' },
  { value: '49.8%', label: 'Terpene Lift', detail: 'Accuracy boost over legacy tooling' },
  { value: '<10s', label: 'Edge Latency', detail: 'Cloudflare Workers global average' }
];

export const partnerSeals: string[] = ['GuideStar · Platinum', 'Cloudflare Impact', 'Stripe Nonprofit'];

export type PlatformHighlight = {
  title: string;
  description: string;
  bullets: string[];
  engines: string;
};

export const platformHighlights: PlatformHighlight[] = [
  {
    title: 'Evidence-Based',
    description: 'Backed by peer-reviewed research with cultural provenance and ChemPath™ + BioPath™ validation.',
    bullets: ['TK validation engine', 'Clinical tiering workflows', 'Bias-aware scoring models'],
    engines: 'ChemPath™ · BioPath™'
  },
  {
    title: 'AI-Powered',
    description: 'GenomePath™ inference models tuned on cannabis-specific corpora and cross-kingdom embeddings.',
    bullets: ['Cloudflare Workers edge lattice', 'Federated learning ready'],
    engines: 'GenomePath™'
  },
  {
    title: 'Comprehensive',
    description: 'RegPath™ keeps terpene, cannabinoid, and polysaccharide dimers in one cockpit with regulatory telemetry.',
    bullets: ['Budtender tablet mode', 'POS telemetry ingestion', 'Regulatory snapshots'],
    engines: 'RegPath™'
  },
  {
    title: "Cloak & Quill Research",
    description: '501(c)(3) nonprofit stewarding NeuroBotanica for public benefit and open scientific access.',
    bullets: ['Public-benefit stewardship', 'Nonprofit research oversight', 'Open evidence distribution to partners'],
    engines: 'Public benefit charter'
  }
];

export type JourneyStep = {
  stage: string;
  title: string;
  description: string;
  kpi: string;
};

export const journeySteps: JourneyStep[] = [
  {
    stage: 'Stage 01',
    title: 'Ingest & Normalize',
    description: 'Upload manuscripts, terpene panels, and TK metadata for harmonization.',
    kpi: '11s ingestion median'
  },
  {
    stage: 'Stage 02',
    title: 'Edge Inference',
    description: 'Workers multi-region lattice predicts viable dimers and terpene deltas.',
    kpi: '235ms edge median'
  },
  {
    stage: 'Stage 03',
    title: 'Bias-Aware Review',
    description: 'Human-in-the-loop review with cultural equity scoring and audit logs.',
    kpi: '92% faster than manual'
  },
  {
    stage: 'Stage 04',
    title: 'Pilot & Deploy',
    description: 'Stripe-backed retainers, budtender tablet rollout, and compliance packet.',
    kpi: '4-week MVP pilots'
  }
];

export type TerpeneMetric = {
  title: string;
  delta: string;
  description: string;
};

export const terpeneMetrics: TerpeneMetric[] = [
  { title: 'Myrcene Sedation Index', delta: '+18%', description: 'Improved nighttime regimen adherence.' },
  { title: 'Limonene Mood Band', delta: '+26%', description: 'Floor staff accuracy on mood pairings.' },
  { title: 'Pinene Focus Arc', delta: '-12%', description: 'Reduced cognitive haze complaints.' }
];

export type CaseBrief = {
  title: string;
  subtitle: string;
  topic: string;
  keyInsight: string;
  solution: string;
  citation: string;
};

export const caseBriefs: CaseBrief[] = [
  {
    title: 'Cannabis Therapeutics (JAMA Review 2026)',
    subtitle: 'Evidence gap across U.S. clinicians',
    topic: 'Therapeutic use of cannabis and cannabinoids',
    keyInsight: 'Only 33% of clinicians feel confident in cannabis knowledge while 86% are requesting structured education.',
    solution: 'NeuroBotanica’s tiered evidence engine (ChemPath™ + BioPath™) delivers clinician-ready packets with TK provenance.',
    citation: 'Hsu M, et al. JAMA. 2026;335(4):345-359.'
  },
  {
    title: 'Psilocybin Analogs & Cellular Longevity',
    subtitle: 'Emory University · Discover Magazine 2026',
    topic: 'Psilocybin analogs and anti-aging outcomes',
    keyInsight: 'Metabolized psilocin extended human cell lifespan by 50%+ and preserved telomeres in vitro.',
    solution: 'Our Compliant Analog Discovery workflow flags non-hallucinogenic variants for regulated therapeutic pipelines.',
    citation: 'Source: Discover Magazine / Emory University research notes.'
  }
];
