import { useEffect } from 'react'

interface ScrollSceneOptions {
  triggerId: string
  targetSelector: string
}

export function useScrollScene(options: ScrollSceneOptions) {
  useEffect(() => {
    // Stub implementation - add scroll scene logic here if needed
    console.log(`Scroll scene for ${options.triggerId} with selector ${options.targetSelector}`)
  }, [options.triggerId, options.targetSelector])
}