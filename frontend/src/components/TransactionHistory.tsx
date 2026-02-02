import { useEffect, useState } from 'react'
import { dispensaryAPI } from '@/utils/api'

interface TransactionItem {
  product_name?: string
  quantity?: number
  [key: string]: unknown
}

interface Transaction {
  transaction_id: string
  created_at: string
  total_amount?: number
  items?: TransactionItem[]
  notes?: string
}

interface TransactionHistoryProps {
  customerId: string
}

export function TransactionHistory({ customerId }: TransactionHistoryProps) {
  const [history, setHistory] = useState<Transaction[]>([])
  const [loading, setLoading] = useState(true)
  const [error, setError] = useState<string | null>(null)

  useEffect(() => {
    if (!customerId) {
      setLoading(false)
      return
    }

    let canceled = false
    setLoading(true)
    setError(null)

    dispensaryAPI
      .getTransactions(customerId)
      .then((response) => {
        if (canceled) {
          return
        }
        setHistory(response.data?.transactions || [])
      })
      .catch((err) => {
        console.error('Failed to load transaction history:', err)
        if (!canceled) {
          setError('Unable to load transaction history right now.')
        }
      })
      .finally(() => {
        if (!canceled) {
          setLoading(false)
        }
      })

    return () => {
      canceled = true
    }
  }, [customerId])

  if (!customerId) {
    return <div className="p-4 text-white/60">Save the profile to view history.</div>
  }

  if (loading) {
    return <div className="p-4 text-white/60">Loading history…</div>
  }

  if (error) {
    return <div className="p-4 text-amber-300 text-sm">{error}</div>
  }

  if (!history.length) {
    return <div className="p-4 text-white/60">No past consultations recorded yet.</div>
  }

  return (
    <div className="space-y-4">
      {history.map((txn) => (
        <div key={txn.transaction_id} className="bg-white/5 border border-white/10 rounded-xl p-4">
          <div className="flex flex-wrap items-center justify-between gap-2 mb-3">
            <div className="text-emerald-300 font-medium">
              {new Date(txn.created_at).toLocaleDateString()}
            </div>
            <div className="text-white/40 text-xs font-mono">{txn.transaction_id}</div>
          </div>

          <div className="space-y-1 text-sm text-white">
            {txn.items?.length ? (
              txn.items.map((item, index) => (
                <div key={`${txn.transaction_id}-${index}`}>• {item.product_name || 'Recommended item'}</div>
              ))
            ) : (
              <div className="text-white/50">No product details recorded.</div>
            )}
          </div>

          {(txn.notes || txn.total_amount) && (
            <div className="mt-3 pt-3 border-t border-white/10 text-sm text-white/70 space-y-1">
              {typeof txn.total_amount === 'number' && (
                <div>Total: ${txn.total_amount.toFixed(2)}</div>
              )}
              {txn.notes && <div className="italic">“{txn.notes}”</div>}
            </div>
          )}
        </div>
      ))}
    </div>
  )
}
