interface OrientationProgressSummaryProps {
  completedCount: number
  total: number
  onReset: () => void
}

export function OrientationProgressSummary({ completedCount, total, onReset }: OrientationProgressSummaryProps) {
  return (
    <div className="mb-8">
      <h2 className="text-2xl font-bold text-white mb-4">Orientation Progress</h2>
      <p className="text-white/70">
        Completed {completedCount} of {total} lessons
      </p>
      <button onClick={onReset} className="btn-secondary mt-2">
        Reset Progress
      </button>
    </div>
  )
}