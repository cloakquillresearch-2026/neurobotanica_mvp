import { useEffect } from 'react'
import { useRouter } from 'next/router'
import Head from 'next/head'
import Link from 'next/link'

export default function BudtenderRedirect() {
  const router = useRouter()
  useEffect(() => {
    // Safely redirect old /budtender path to main app
    router.replace('/')
  }, [router])

  return (
    <div className="min-h-screen flex items-center justify-center">
      <Head>
        <title>Redirecting…</title>
      </Head>
      <div className="text-center">
        <h1 className="text-xl font-semibold">Redirecting to Budtender App</h1>
        <p className="text-sm text-gray-400 mt-2">If you are not redirected, <Link href="/" className="underline">click here</Link>.</p>
      </div>
    </div>
  )
}
