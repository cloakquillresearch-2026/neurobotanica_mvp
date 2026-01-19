import { useEffect } from 'react'

export function useHeroTimeline(id: string) {
  useEffect(() => {
    // Stub implementation - add timeline logic here if needed
    console.log(`Hero timeline for ${id}`)
  }, [id])
}