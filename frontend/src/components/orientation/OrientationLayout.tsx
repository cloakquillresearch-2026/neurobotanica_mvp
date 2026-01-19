import { ReactNode } from 'react'

interface OrientationLayoutProps {
  children: ReactNode
}

export function OrientationLayout({ children }: OrientationLayoutProps) {
  return (
    <div className="min-h-screen bg-gradient-to-br from-slate-900 via-slate-800 to-slate-900">
      <div className="container mx-auto px-4 py-8 max-w-6xl">
        {children}
      </div>
    </div>
  )
}