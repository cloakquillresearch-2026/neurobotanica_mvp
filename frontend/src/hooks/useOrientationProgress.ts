import { useState, useCallback, useEffect } from 'react'

interface LessonProgress {
  completed: boolean
  completedAt?: string
}

interface OrientationProgress {
  [lessonId: string]: LessonProgress
}

export function useOrientationProgress(lessonIds: string[]) {
  const [progress, setProgress] = useState<OrientationProgress>({})

  const completedCount = Object.values(progress).filter(p => p.completed).length
  const total = lessonIds.length
  const isComplete = completedCount === total

  // Automatically set completion flag and redirect when all lessons are complete
  useEffect(() => {
    // Only run on client side
    if (typeof window === 'undefined') return;
    
    if (isComplete) {
      localStorage.setItem('orientation_completed', 'true');
      // Show success message and redirect after a brief delay
      setTimeout(() => {
        window.location.href = '/';
      }, 2000);
    }
  }, [isComplete]);

  const markLessonComplete = useCallback((lessonId: string) => {
    setProgress(prev => ({
      ...prev,
      [lessonId]: {
        completed: true,
        completedAt: new Date().toISOString()
      }
    }))
  }, [])

  const resetProgress = useCallback(() => {
    setProgress({})
  }, [])

  return {
    progress,
    completedCount,
    total,
    isComplete,
    markLessonComplete,
    resetProgress
  }
}