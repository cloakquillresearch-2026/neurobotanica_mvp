import Image from 'next/image'
import Link from 'next/link'
import type { OrientationLesson } from '@/content/orientation/lessons'
import { useHeroTimeline } from '@/hooks/useHeroTimeline'
import { useScrollScene } from '@/hooks/useScrollScene'

interface OrientationMediaHighlightProps {
  lessons: OrientationLesson[]
}

const HERO_STATS = [
  { label: 'Apothecary runs/day', value: '48', note: 'Edge simulated consults' },
  { label: 'Sensor inputs', value: '16', note: 'Biomarkers + terpene cues' },
  { label: 'Latency budget', value: '<200 ms', note: 'Cloudflare lab orchestration' },
]

export function OrientationMediaHighlight({ lessons }: OrientationMediaHighlightProps) {
  const labLesson = lessons.find((lesson) => lesson.id === 'tablet-walkthrough')
  useHeroTimeline('apothecary-hero')
  useScrollScene({ triggerId: 'apothecary-hero', targetSelector: '[data-animate="stagger"] > *' })

  return (
    <section className="space-y-10">
      <article
        className="relative overflow-hidden glass-card border border-emerald-400/20 bg-gradient-to-br from-[#02131f]/95 via-[#022e27]/90 to-[#0c5f41]/85 p-6 lg:p-10"
        aria-labelledby="apothecary-hero-title"
        data-animate-id="apothecary-hero"
      >
        <div className="absolute -right-24 -bottom-24 w-[28rem] h-[28rem] bg-emerald-500/15 blur-[160px]" aria-hidden="true" />
        <div className="absolute -left-16 -top-20 w-72 h-72 bg-sky-500/15 blur-[130px]" aria-hidden="true" />
        <div className="absolute inset-0 opacity-40" aria-hidden="true">
          <div className="absolute inset-0 bg-[radial-gradient(circle_at_top,_rgba(255,255,255,0.08),_transparent_55%)]" />
        </div>

        <div className="grid gap-10 lg:grid-cols-[minmax(0,1.2fr)_minmax(280px,1fr)] items-center">
          <div className="space-y-8 max-w-2xl" data-animate="fade-up" data-animate-target="apothecary-hero">
            <div className="space-y-4">
              <p className="inline-flex items-center gap-2 text-xs font-semibold uppercase tracking-[0.4em] text-emerald-200/90">
                Apothecary Lab Mode
                <span className="h-px w-10 bg-emerald-200/60" aria-hidden="true" />
              </p>
              <h2 id="apothecary-hero-title" className="text-4xl lg:text-5xl font-semibold text-white leading-tight">
                NeuroBotanica unites apothecary ritual with compliance intelligence.
              </h2>
              <p className="text-white/80 text-base">
                {labLesson?.summary || 'Interactive walkthrough that pairs biomarker cues with dispensary instrumentation.'}
              </p>
            </div>

            <div className="grid gap-4 sm:grid-cols-3" data-animate="stagger" data-animate-target="apothecary-hero">
              {HERO_STATS.map((stat) => (
                <div key={stat.label} className="rounded-2xl border border-white/15 bg-white/5 p-4 backdrop-blur">
                  <dt className="text-[0.65rem] uppercase tracking-[0.3em] text-white/60">{stat.label}</dt>
                  <dd className="text-white text-3xl font-semibold leading-tight mt-2">{stat.value}</dd>
                  <p className="text-white/60 text-xs mt-1">{stat.note}</p>
                </div>
              ))}
            </div>

            <div className="grid gap-3 sm:grid-cols-3" data-animate="fade-up" data-animate-delay="100" data-animate-target="apothecary-hero">
              {(labLesson?.objectives || [
                'Blend ethnobotanical context with TS-PS-001 instrumentation cues.',
                'Practice the four-step ritual: scan, interpret, blend, certify.',
                'Document sensory notes for Nevada compliance audits.'
              ]).map((objective, index) => (
                <div key={objective} className="rounded-xl border border-white/10 bg-black/20 p-3">
                  <p className="text-emerald-200/90 text-xs mb-1">Step 0{index + 1}</p>
                  <p className="text-white/80 text-sm leading-relaxed">{objective}</p>
                </div>
              ))}
            </div>

            <div className="flex flex-col sm:flex-row gap-3" data-animate="fade-up" data-animate-delay="160" data-animate-target="apothecary-hero">
              <Link href="/?sandbox=1" className="btn-primary shadow-lg shadow-emerald-500/30 px-6 py-3 text-base">
                Launch Lab Sandbox
              </Link>
              <button
                type="button"
                className="btn-secondary bg-white/5 border-white/30 px-6 py-3 text-base"
                data-animate="pulse"
                data-animate-delay="220"
              >
                Queue instrumentation reel
              </button>
            </div>
          </div>

          <div
            className="relative aspect-[4/3] w-full max-w-md mx-auto"
            data-animate="float"
            data-animate-delay="140"
            data-animate-target="apothecary-hero"
          >
            <div className="absolute inset-0 rounded-[40px] bg-gradient-to-br from-white/15 to-white/5" aria-hidden="true" />
            <Image
              src="/orientation/apothecary-lab.svg"
              alt="Hybrid DNA helix made of circuitry and botanical glyphs"
              fill
              className="object-contain drop-shadow-[0_30px_60px_rgba(15,255,200,0.35)]"
              priority
            />
            <div className="absolute inset-4 rounded-[32px] border border-white/40" aria-hidden="true" />
            <div className="absolute inset-y-6 inset-x-8 border border-emerald-200/30 rounded-[28px] animate-pulse" aria-hidden="true" />
          </div>
        </div>
      </article>

      <div className="grid gap-6 lg:grid-cols-3">
        {lessons.map((lesson) => (
          <figure key={lesson.id} className="glass-card p-4 space-y-4 border border-white/10" aria-labelledby={`${lesson.id}-media-title`}>
            <div className="relative w-full h-44 rounded-2xl overflow-hidden bg-black/30">
              {lesson.media.type === 'image' ? (
                <Image
                  src={lesson.media.src}
                  alt={lesson.media.alt}
                  fill
                  sizes="(min-width: 1024px) 30vw, 90vw"
                  className="object-cover scale-[1.02]"
                  priority={lesson.id === 'welcome'}
                />
              ) : (
                <video
                  className="w-full h-full object-cover"
                  controls
                  preload="metadata"
                  poster={lesson.media.src}
                >
                  <source src={lesson.media.src} />
                </video>
              )}
            </div>
            <figcaption className="space-y-1">
              <p id={`${lesson.id}-media-title`} className="text-white font-semibold text-lg">
                {lesson.title}
              </p>
              {lesson.media.caption && <p className="text-white/70 text-sm leading-relaxed">{lesson.media.caption}</p>}
            </figcaption>
          </figure>
        ))}
      </div>
    </section>
  )
}
