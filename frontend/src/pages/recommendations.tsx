import Head from 'next/head'
import { useState } from 'react'
import { CustomerSearch } from '@/components/CustomerSearch'
import { CustomerProfile } from '@/components/CustomerProfile'
import { ProductRecommendations } from '@/components/ProductRecommendations'
import type { Recommendation } from '@/components/ProductRecommendations'
import type { CustomerProfileData } from '@/types/customer'

export default function RecommendationsPage() {
  const [customer, setCustomer] = useState<CustomerProfileData | null>(null)
  const [recommendations, setRecommendations] = useState<Recommendation[]>([])

  return (
    <div className="min-h-screen">
      <Head>
        <title>Recommendations — NeuroBotanica</title>
      </Head>
      <main className="container mx-auto px-4 py-8">
        <h1 className="text-2xl font-bold text-white mb-4">Recommendation Results</h1>

        <div className="grid grid-cols-1 lg:grid-cols-3 gap-6">
          <div className="lg:col-span-1 glass-card p-6">
            <h2 className="font-semibold text-white mb-3">1. Select Customer</h2>
            <CustomerSearch onCustomerSelect={setCustomer} />

            {customer && (
              <div className="mt-6">
                <h3 className="font-semibold text-white mb-2">Customer Profile</h3>
                <CustomerProfile customer={customer} onProfileUpdate={setCustomer} />
              </div>
            )}
          </div>

          <div className="lg:col-span-2 glass-card p-6">
            <h2 className="font-semibold text-white mb-3">2. Recommendations</h2>
            {customer ? (
              <ProductRecommendations customer={customer} recommendations={recommendations} onRecommendationsUpdate={setRecommendations} />
            ) : (
              <p className="text-white/60">Select a customer to view recommendations.</p>
            )}
          </div>
        </div>
      </main>
    </div>
  )
}
