export type OrientationStep = {
  title: string
  description: string
  tip?: string
}

export type OrientationMedia = {
  type: 'image' | 'video' | 'iframe'
  src: string
  alt: string
  caption?: string
  poster?: string
  fit?: 'cover' | 'contain'
}

export type OrientationLesson = {
  id: string
  title: string
  duration: string
  summary: string
  objectives: string[]
  steps: OrientationStep[]
  media: OrientationMedia
}

export const ORIENTATION_LESSONS: OrientationLesson[] = [
  // LESSON 1: NeuroBotanica Workflow with Canva Video
  {
    id: 'workflow',
    title: 'NeuroBotanica in 90 Seconds',
    duration: '2 min',
    summary: 'High-level walkthrough covering evidence, terpene intelligence, and why TS-PS-001 matters at the counter.',
    objectives: [
      'Understand the mission and Nevada compliance guardrails',
      'See how terpene + biomarker pairing shortens consultations',
      'Remember the three talking points for customers'
    ],
    steps: [
      {
        title: 'Watch the platform overview',
        description: 'This 90-second video explains the NeuroBotanica value proposition, TS-PS-001 overview, and Nevada market impact.'
      },
      {
        title: 'Review key statistics',
        description: 'Card highlights 505+ studies, accuracy lifts, and pilot KPIs for credibility.'
      },
      {
        title: 'Commitment check',
        description: 'Simple question: "Can you explain NeuroBotanica in one line?" to reinforce retention.'
      }
    ],
    media: {
      type: 'iframe',
      src: 'https://www.canva.com/design/DAG-oazPd7M/Iyp4VI6UeyLtaYNBfAbyEQ/watch?embed',
      alt: 'NeuroBotanica platform workflow overview',
      caption: 'Watch how clinical evidence powers personalized recommendations'
    }
  },

  // LESSON 2: Case Study Sprint - ONE CARD with 3 Personas
  {
    id: 'case-studies',
    title: 'Case Study Sprint',
    duration: '3 min',
    summary: 'Scenario-based exercise using three personas for pain, anxiety, and nausea management.',
    objectives: [
      'Map conditions to terpene/adjuvant narratives',
      'Practice phrasing recommendations with clinical language',
      'Capture confidence ratings for analytics'
    ],
    steps: [
      {
        title: 'Choose a persona',
        description: 'Three personas available: Sarah (chronic pain), Marcus (anxiety), Janet (nausea). Cards include symptom cues, preferred consumption methods, and concerns.'
      },
      {
        title: 'Answer branching prompts',
        description: 'Learners select best-fit guidance and see instant feedback with rationale. Sarah needs non-psychoactive pain relief (β-caryophyllene + CBD). Marcus wants natural anxiety management (linalool + limonene). Janet requires anti-nausea support (THC + CBDA + ginger).'
      },
      {
        title: 'Submit confidence pulse',
        description: 'Thumbs up/down plus optional note stored for the manager dashboard.'
      }
    ],
    media: {
      type: 'image',
      src: 'data:image/svg+xml,%3Csvg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 800 400"%3E%3Crect fill="%2323394d" width="800" height="400"/%3E%3Cg%3E%3Crect x="20" y="80" width="220" height="280" fill="%231a2332" rx="8"/%3E%3Ccircle cx="130" cy="140" r="30" fill="%2300A86B"/%3E%3Ctext x="130" y="190" font-family="Arial" font-size="18" fill="%23fff" text-anchor="middle" font-weight="bold"%3ESarah, 45%3C/text%3E%3Ctext x="130" y="215" font-family="Arial" font-size="14" fill="%23a0aec0" text-anchor="middle"%3EChronic Pain%3C/text%3E%3Ctext x="130" y="240" font-family="Arial" font-size="12" fill="%2300A86B" text-anchor="middle"%3Eβ-caryophyllene%3C/text%3E%3Ctext x="130" y="260" font-family="Arial" font-size="12" fill="%2300A86B" text-anchor="middle"%3E+ CBD%3C/text%3E%3C/g%3E%3Cg%3E%3Crect x="290" y="80" width="220" height="280" fill="%231a2332" rx="8"/%3E%3Ccircle cx="400" cy="140" r="30" fill="%2300A86B"/%3E%3Ctext x="400" y="190" font-family="Arial" font-size="18" fill="%23fff" text-anchor="middle" font-weight="bold"%3EMarcus, 32%3C/text%3E%3Ctext x="400" y="215" font-family="Arial" font-size="14" fill="%23a0aec0" text-anchor="middle"%3EAnxiety%3C/text%3E%3Ctext x="400" y="240" font-family="Arial" font-size="12" fill="%2300A86B" text-anchor="middle"%3ELinalool%3C/text%3E%3Ctext x="400" y="260" font-family="Arial" font-size="12" fill="%2300A86B" text-anchor="middle"%3E+ Limonene%3C/text%3E%3C/g%3E%3Cg%3E%3Crect x="560" y="80" width="220" height="280" fill="%231a2332" rx="8"/%3E%3Ccircle cx="670" cy="140" r="30" fill="%2300A86B"/%3E%3Ctext x="670" y="190" font-family="Arial" font-size="18" fill="%23fff" text-anchor="middle" font-weight="bold"%3EJanet, 58%3C/text%3E%3Ctext x="670" y="215" font-family="Arial" font-size="14" fill="%23a0aec0" text-anchor="middle"%3ENausea%3C/text%3E%3Ctext x="670" y="240" font-family="Arial" font-size="12" fill="%2300A86B" text-anchor="middle"%3ETHC + CBDA%3C/text%3E%3Ctext x="670" y="260" font-family="Arial" font-size="12" fill="%2300A86B" text-anchor="middle"%3E+ Ginger%3C/text%3E%3C/g%3E%3Ctext x="400" y="40" font-family="Arial" font-size="20" fill="%23fff" text-anchor="middle" font-weight="bold"%3EThree Case Studies%3C/text%3E%3C/svg%3E',
      alt: 'Three case study personas: Sarah (pain), Marcus (anxiety), Janet (nausea)',
      caption: 'Practice with three real-world scenarios and get instant feedback',
      fit: 'contain'
    }
  },

  // LESSON 3: Terpenes & Adjuvants Education
  {
    id: 'terpenes-adjuvants',
    title: 'Terpenes & Adjuvants',
    duration: '3 min',
    summary: 'Learn the science behind terpenes and how adjuvants enhance bioavailability.',
    objectives: [
      'Understand what terpenes are and their therapeutic effects',
      'Learn the entourage effect — how terpenes synergize with cannabinoids',
      'Know key terpenes: linalool, β-caryophyllene, myrcene, limonene',
      'Understand adjuvants and their role in enhancing efficacy'
    ],
    steps: [
      {
        title: 'What are terpenes?',
        description: 'Aromatic compounds in cannabis that contribute to smell and therapeutic effects. Over 100 terpenes identified in cannabis, each with unique properties.',
        tip: 'Think of terpenes like essential oils — they give plants their scent and medicinal qualities.'
      },
      {
        title: 'The entourage effect',
        description: 'Terpenes work synergistically with cannabinoids (THC, CBD) to enhance or modify effects. For example: myrcene increases THC absorption across blood-brain barrier.',
        tip: 'This is why whole-plant extracts often work better than isolated cannabinoids.'
      },
      {
        title: 'Key terpenes to know',
        description: 'Linalool (calming, anti-anxiety), β-caryophyllene (anti-inflammatory, CB2 agonist), Myrcene (sedation, muscle relaxation), Limonene (mood elevation, anti-anxiety).'
      },
      {
        title: 'What are adjuvants?',
        description: 'Non-cannabinoid compounds that enhance bioavailability or potentiate effects. Examples: Piperine (black pepper) increases cannabinoid absorption, MCT oil improves fat-soluble cannabinoid delivery, Ginger provides anti-nausea synergy.'
      },
      {
        title: 'Practical application',
        description: 'How to read a terpene profile: Look at lab results showing terpene percentages. Match dominant terpenes to customer goals (e.g., high linalool for anxiety). Explain the "why" behind recommendations using terpene science.'
      }
    ],
    media: {
      type: 'image',
      src: 'data:image/svg+xml,%3Csvg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 600 400"%3E%3Crect fill="%2323394d" width="600" height="400"/%3E%3Ctext x="300" y="40" font-family="Arial" font-size="24" fill="%2300A86B" text-anchor="middle" font-weight="bold"%3ETerpenes %26 Adjuvants%3C/text%3E%3Cg%3E%3Crect x="40" y="80" width="250" height="120" fill="%231a2332" rx="8"/%3E%3Ctext x="165" y="105" font-family="Arial" font-size="16" fill="%23fff" text-anchor="middle" font-weight="bold"%3EKey Terpenes%3C/text%3E%3Ctext x="165" y="130" font-family="Arial" font-size="13" fill="%2300A86B" text-anchor="middle"%3ELinalool - Calming%3C/text%3E%3Ctext x="165" y="150" font-family="Arial" font-size="13" fill="%2300A86B" text-anchor="middle"%3EMyrcene - Sedation%3C/text%3E%3Ctext x="165" y="170" font-family="Arial" font-size="13" fill="%2300A86B" text-anchor="middle"%3ELimonene - Mood%3C/text%3E%3Ctext x="165" y="190" font-family="Arial" font-size="13" fill="%2300A86B" text-anchor="middle"%3Eβ-caryophyllene - Pain%3C/text%3E%3C/g%3E%3Cg%3E%3Crect x="310" y="80" width="250" height="120" fill="%231a2332" rx="8"/%3E%3Ctext x="435" y="105" font-family="Arial" font-size="16" fill="%23fff" text-anchor="middle" font-weight="bold"%3EAdjuvants%3C/text%3E%3Ctext x="435" y="130" font-family="Arial" font-size="13" fill="%2300A86B" text-anchor="middle"%3EPiperine - Absorption%3C/text%3E%3Ctext x="435" y="150" font-family="Arial" font-size="13" fill="%2300A86B" text-anchor="middle"%3EMCT Oil - Bioavailability%3C/text%3E%3Ctext x="435" y="170" font-family="Arial" font-size="13" fill="%2300A86B" text-anchor="middle"%3EGinger - Anti-nausea%3C/text%3E%3C/g%3E%3Cg%3E%3Crect x="40" y="220" width="520" height="140" fill="%231a2332" rx="8"/%3E%3Ctext x="300" y="245" font-family="Arial" font-size="16" fill="%23fff" text-anchor="middle" font-weight="bold"%3EThe Entourage Effect%3C/text%3E%3Ctext x="300" y="270" font-family="Arial" font-size="13" fill="%23a0aec0" text-anchor="middle"%3ETerpenes + Cannabinoids work together synergistically%3C/text%3E%3Ctext x="300" y="295" font-family="Arial" font-size="13" fill="%23a0aec0" text-anchor="middle"%3EWhole-plant extracts often more effective than isolates%3C/text%3E%3Ctext x="300" y="320" font-family="Arial" font-size="13" fill="%23a0aec0" text-anchor="middle"%3EExample: Myrcene increases THC blood-brain barrier crossing%3C/text%3E%3Ctext x="300" y="345" font-family="Arial" font-size="13" fill="%2300A86B" text-anchor="middle"%3EThis is why terpene profiles matter for recommendations%3C/text%3E%3C/g%3E%3C/svg%3E',
      alt: 'Terpenes and adjuvants educational overview',
      caption: 'Master the science behind personalized cannabis recommendations',
      fit: 'contain'
    }
  }
]
