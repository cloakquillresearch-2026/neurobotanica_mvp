import { useCallback, useEffect } from 'react'
import { getFirebaseAnalytics } from '@/lib/firebaseClient'
import Head from 'next/head'
import Link from 'next/link'
import { useRouter } from 'next/router'
import { OrientationLayout } from '@/components/orientation/OrientationLayout'
import { OrientationLessonCard } from '@/components/orientation/OrientationLessonCard'
import { OrientationMediaHighlight } from '@/components/orientation/OrientationMediaHighlight'
import { OrientationProgressSummary } from '@/components/orientation/OrientationProgressSummary'
import { ORIENTATION_LESSONS } from '@/content/orientation/lessons'
import { useOrientationProgress } from '@/hooks/useOrientationProgress'

const LESSON_IDS = ORIENTATION_LESSONS.map((lesson) => lesson.id)

export default function OrientationPage() {
  const router = useRouter()
  const { progress, completedCount, total, isComplete, markLessonComplete, resetProgress } = useOrientationProgress(LESSON_IDS)
  const handleResumeLessons = useCallback(() => {
    if (typeof window === 'undefined') return
    window.scrollTo({ top: 0, behavior: 'smooth' })
  }, [])
  const handleLaunchSandbox = useCallback(() => {
    if (typeof window !== 'undefined') {
      window.localStorage.setItem('nb_sandbox_mode', '1')
    }
    router.push('/?sandbox=1')
  }, [router])

  useEffect(() => {
    getFirebaseAnalytics().catch(() => {
      /* analytics optional */
    })
  }, [])

  return (
    <OrientationLayout>
      <Head>
        <title>NeuroBotanica Orientation</title>
        <meta name="description" content="Budtender onboarding course for NeuroBotanica" />
      </Head>

      <OrientationProgressSummary completedCount={completedCount} total={total} onReset={resetProgress} />

      <OrientationMediaHighlight lessons={ORIENTATION_LESSONS} />

      <div className="grid gap-6 lg:grid-cols-2">
        {ORIENTATION_LESSONS.map((lesson) => (
          <OrientationLessonCard
            key={lesson.id}
            lesson={lesson}
            completed={Boolean(progress[lesson.id]?.completed)}
            onComplete={markLessonComplete}
          />
        ))}
      </div>

      <section className="glass-card p-8 flex flex-col lg:flex-row gap-8 items-start">
        <div className="flex-1 space-y-3">
          <p className="text-sm uppercase tracking-wide text-white/60">Next Steps</p>
          <h2 className="text-3xl font-bold text-white leading-snug">
            {isComplete ? 'Orientation complete. Ready to launch a consultation?' : 'Keep going to unlock sandbox mode.'}
          </h2>
          <p className="text-white/70 text-base max-w-2xl">
            Your progress automatically saves on this tablet. Once you finish the lessons, you can access certification quizzes and
            sandbox practice directly from the Budtender app.
          </p>
        </div>
        <div className="flex flex-col gap-3 w-full max-w-sm lg:w-auto">
          <button type="button" className="btn-primary text-base py-3" onClick={handleLaunchSandbox}>
            Launch Sandbox Practice
          </button>
          <Link href="/" className="btn-secondary text-center text-base py-3">
            Return to Budtender App
          </Link>
          {!isComplete && (
            <button type="button" className="btn-secondary text-base py-3" onClick={handleResumeLessons}>
              Resume Lessons
            </button>
          )}
        </div>
      </section>
    </OrientationLayout>
  )
}
