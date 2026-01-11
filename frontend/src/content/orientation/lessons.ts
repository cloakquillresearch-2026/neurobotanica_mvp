export type OrientationStep = {
  title: string
  description: string
  tip?: string
}

export type OrientationMedia = {
  type: 'image' | 'video'
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
  {
    id: 'welcome',
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
        title: 'Watch the brand story clip',
        description: '45-second video that frames cultural preservation plus modern compliance.',
        tip: 'Cache the MP4 so it plays smoothly on low bandwidth.'
      },
      {
        title: 'Review outcome stats',
        description: 'Card highlights 505+ studies, accuracy lifts, and pilot KPIs for credibility.'
      },
      {
        title: 'Commitment check',
        description: 'Simple question: "Can you explain NeuroBotanica in one line?" to reinforce retention.'
      }
    ],
    media: {
      type: 'video',
      src: 'https://storage.googleapis.com/gtv-videos-bucket/sample/ForBiggerJoyrides.mp4',
      poster: '/orientation/welcome-story.svg',
      alt: 'NeuroBotanica brand spot covering cultural preservation and compliance',
      caption: 'Clip + talking points to explain tradition-to-tech pipeline.'
    }
  },
  {
    id: 'tablet-walkthrough',
    title: 'Tablet Workflow Practice',
    duration: '3 min',
    summary: 'Hands-on overlay guiding budtenders through search, biomarker entry, and interpreting TS-PS-001 outputs.',
    objectives: [
      'Locate each panel quickly during live consultations',
      'Know what to do when biomarkers are missing',
      'Feel confident explaining the synergy card'
    ],
    steps: [
      {
        title: 'Launch sandbox mode',
        description: 'Creates a mock customer record with preset conditions so trainees cannot impact production data.'
      },
      {
        title: 'Follow the glow path',
        description: 'Guided overlay highlights each component with short copy blocks and keyboard focus states.'
      },
      {
        title: 'Reset + repeat',
        description: 'One tap reset clears sandbox data allowing unlimited reps.'
      }
    ],
    media: {
      type: 'image',
      src: '/orientation/tablet-walkthrough.svg',
      alt: 'Tablet walkthrough illustration',
      caption: 'Overlay-driven mock consultation that highlights each panel.',
      fit: 'contain'
    }
  },
  {
    id: 'case-sprint',
    title: 'Case Study Sprint',
    duration: '3 min',
    summary: 'Scenario-based exercise using personas for pain, anxiety, and metabolic goals.',
    objectives: [
      'Map conditions to terpene/adjuvant narratives',
      'Practice phrasing recommendations with clinical language',
      'Capture confidence ratings for analytics'
    ],
    steps: [
      {
        title: 'Choose a persona',
        description: 'Cards include symptom cues, preferred consumption methods, and concerns.'
      },
      {
        title: 'Answer branching prompts',
        description: 'Learners select best-fit guidance and see instant feedback with rationale.'
      },
      {
        title: 'Submit confidence pulse',
        description: 'Thumbs up/down plus optional note stored for the manager dashboard.'
      }
    ],
    media: {
      type: 'image',
      src: '/orientation/case-sprint.svg',
      alt: 'Case study sprint persona cards',
      caption: 'Three personas illustrate how to tailor terpene + adjuvant guidance.',
      fit: 'contain'
    }
  }
]
