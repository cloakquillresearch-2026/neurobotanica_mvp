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
  { value: '500+', label: 'Studies Parsed', detail: 'Peer-reviewed + ethnobotanical corpus' },
  { value: '49.8%', label: 'Terpene Lift', detail: 'Accuracy boost over legacy tooling' },
  { value: '<10s', label: 'Edge Latency', detail: 'Cloudflare Workers global average' }
];

export const partnerSeals: string[] = ['GuideStar · Silver', 'Cloudflare Impact', 'Stripe Nonprofit'];

export type PlatformHighlight = {
  title: string;
  description: string;
  bullets: string[];
  cta: string;
};

export const platformHighlights: PlatformHighlight[] = [
  {
    title: 'Evidence-Based',
    description: 'Backed by peer-reviewed research with cultural provenance.',
    bullets: ['TK validation engine', 'Clinical tiering', 'Bias-aware scoring'],
    cta: 'See validation matrix'
  },
  {
    title: 'AI-Powered',
    description: 'Cross-kingdom inference models tuned on cannabis-specific corpora.',
    bullets: ['Cloudflare Workers edge', 'Federated learning ready', 'HIPAA-aligned logging'],
    cta: 'Inspect architecture'
  },
  {
    title: 'Comprehensive',
    description: 'Terpene, cannabinoid, and polysaccharide dimers in one cockpit.',
    bullets: ['Budtender tablet mode', 'POS telemetry ingestion', 'Regulatory snapshots'],
    cta: 'Preview cockpit'
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
  description: string;
  results: string;
  tags: string[];
  outcome: string;
  detail: string;
};

export const caseBriefs: CaseBrief[] = [
  {
    title: 'Clinical Study',
    description: 'Analyzed therapeutic profiles',
    results: '95% accuracy in predictions',
    tags: ['Clinical', 'Research', 'Validated'],
    outcome: '95% accuracy achieved',
    detail: 'Comprehensive analysis of botanical therapeutic profiles'
  }
];
