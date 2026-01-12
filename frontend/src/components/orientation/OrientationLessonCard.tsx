import type { OrientationLesson } from '@/content/orientation/lessons'

interface OrientationLessonCardProps {
  lesson: OrientationLesson
  completed?: boolean
  onComplete?: (lessonId: string) => void
}

export function OrientationLessonCard({ lesson, completed, onComplete }: OrientationLessonCardProps) {
  return (
    <article className="glass-card p-6 space-y-5" aria-live="polite">
      <header className="flex flex-col gap-3 sm:flex-row sm:items-start sm:justify-between">
        <div className="space-y-1">
          <p className="text-[0.65rem] uppercase tracking-[0.4em] text-emerald-200/80">Lesson</p>
          <h3 className="text-2xl font-bold text-white leading-tight">{lesson.title}</h3>
          <p className="text-white/70 text-sm max-w-prose">{lesson.summary}</p>
        </div>
        <div className="flex items-center gap-2">
          <span className="px-3 py-1 rounded-full text-xs font-semibold bg-black/40 text-white/80 border border-white/10">
            {lesson.duration}
          </span>
          <span className={`px-3 py-1 rounded-full text-xs font-semibold ${completed ? 'bg-emerald-500/20 text-emerald-200 border border-emerald-400/30' : 'bg-white/5 text-white/60 border border-white/10'}`}>
            {completed ? 'Completed' : 'Pending'}
          </span>
        </div>
      </header>

      <div className="grid gap-3 lg:grid-cols-[minmax(0,1fr)_minmax(0,1fr)]">
        <section className="bg-black/25 border border-white/5 rounded-2xl p-4 space-y-3">
          <p className="text-sm font-semibold text-white tracking-wide">Objectives</p>
          <ul className="text-white/80 text-sm space-y-2">
            {lesson.objectives.map((objective) => (
              <li key={objective} className="flex gap-2">
                <span className="mt-[6px] h-1.5 w-1.5 rounded-full bg-emerald-300" aria-hidden="true" />
                <span>{objective}</span>
              </li>
            ))}
          </ul>
        </section>

        <section className="border border-white/10 rounded-2xl p-4 bg-gradient-to-br from-white/5 via-transparent to-white/0">
          <p className="text-sm font-semibold text-white mb-3">Key Actions</p>
          <ol className="space-y-3">
            {lesson.steps.map((step, index) => (
              <li key={step.title} className="flex gap-3">
                <div className="w-9 h-9 rounded-full bg-white/10 text-white flex items-center justify-center font-semibold text-sm">
                  {index + 1}
                </div>
                <div className="space-y-1">
                  <p className="text-white font-semibold text-sm">{step.title}</p>
                  <p className="text-white/70 text-sm leading-relaxed">{step.description}</p>
                  {step.tip && <p className="text-emerald-300 text-xs">Tip: {step.tip}</p>}
                </div>
              </li>
            ))}
          </ol>
        </section>
      </div>

      {onComplete && !completed && (
        <div className="flex flex-wrap items-center gap-3">
          <button
            type="button"
            className="btn-primary text-sm px-5 py-2"
            onClick={() => onComplete(lesson.id)}
          >
            Mark Lesson Complete
          </button>
          <p className="text-white/60 text-xs">Progress syncs automatically to the budtender tablet.</p>
        </div>
      )}
    </article>
  )
}
