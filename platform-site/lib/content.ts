// Landing page content
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

export const heroMetrics = [
  { value: '500+', label: 'Studies' },
  { value: '50+', label: 'Terpenes' },
  { value: '99%', label: 'Accuracy' }
];

export const partnerSeals = [
  { name: 'GuideStar', level: 'Silver' }
];

export const platformHighlights = [
  { title: 'Evidence-Based', description: 'Backed by peer-reviewed research' },
  { title: 'AI-Powered', description: 'Advanced machine learning analysis' },
  { title: 'Comprehensive', description: 'Full terpene profile analysis' }
];

export const journeySteps = [
  { step: 1, title: 'Upload', description: 'Upload your botanical data' },
  { step: 2, title: 'Analyze', description: 'AI analyzes terpene profiles' },
  { step: 3, title: 'Insights', description: 'Get therapeutic recommendations' }
];

export const terpeneMetrics = [
  { name: 'Myrcene', description: 'Sedative effects' },
  { name: 'Limonene', description: 'Mood elevation' },
  { name: 'Pinene', description: 'Focus and alertness' }
];

export type CaseBrief = {
  title: string;
  description: string;
  results: string;
  tags: string[];
};

export const caseBriefs: CaseBrief[] = [
  {
    title: 'Clinical Study',
    description: 'Analyzed therapeutic profiles',
    results: '95% accuracy in predictions',
    tags: ['Clinical', 'Research', 'Validated']
  }
];
