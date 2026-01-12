import { useState } from 'react'
import { ShoppingCart, Minus, Plus, Trash2, CreditCard } from 'lucide-react'
import { dispensaryAPI } from '@/utils/api'
import type { CustomerProfileData, TransactionPayload } from '@/types/customer'

interface CartItem {
  product_id: string
  product_name: string
  price: number
  quantity: number
  thc_percent: number
  cbd_percent: number
}

interface CartProps {
  items: CartItem[]
  onUpdateCart: (items: CartItem[]) => void
  customer: CustomerProfileData | null
}

export function Cart({ items, onUpdateCart, customer }: CartProps) {
  const [isProcessing, setIsProcessing] = useState(false)

  const total = items.reduce((sum, item) => sum + item.price * item.quantity, 0)
  const totalTHC = items.reduce((sum, item) => sum + item.thc_percent * item.quantity, 0)
  const totalCBD = items.reduce((sum, item) => sum + item.cbd_percent * item.quantity, 0)

  const handleUpdateQuantity = (productId: string, quantity: number) => {
    if (quantity === 0) {
      onUpdateCart(items.filter(item => item.product_id !== productId))
    } else {
      onUpdateCart(items.map(item =>
        item.product_id === productId
          ? { ...item, quantity }
          : item
      ))
    }
  }

  const handleCheckout = async () => {
    if (!customer || items.length === 0) return

    setIsProcessing(true)
    try {
      const transaction: TransactionPayload = {
        customer_id: customer.customer_id,
        items: items.map(({ product_id, price, quantity, thc_percent, cbd_percent }) => ({
          product_id,
          price,
          quantity,
          thc_percent,
          cbd_percent,
        })),
        total_amount: total,
        timestamp: new Date().toISOString(),
      }

      await dispensaryAPI.createTransaction(transaction)

      // Clear cart and show success
      onUpdateCart([])
      alert('Transaction completed successfully!')
    } catch (error) {
      console.error('Checkout failed:', error)
      alert('Checkout failed. Please try again.')
    } finally {
      setIsProcessing(false)
    }
  }

  if (items.length === 0) {
    return (
      <div className="bg-white rounded-lg shadow p-4">
        <div className="text-center text-gray-500 py-8">
          <ShoppingCart size={48} className="mx-auto mb-4 opacity-50" />
          <p>Cart is empty</p>
          <p className="text-sm">Add products to get started</p>
        </div>
      </div>
    )
  }

  return (
    <div className="bg-white rounded-lg shadow p-4">
      <h2 className="text-lg font-semibold mb-4 flex items-center gap-2">
        <ShoppingCart size={20} />
        Cart ({items.length} items)
      </h2>

      <div className="space-y-3 mb-4 max-h-64 overflow-y-auto">
        {items.map((item) => (
          <div key={item.product_id} className="border border-gray-200 rounded-lg p-3">
            <div className="flex justify-between items-start mb-2">
              <div className="flex-1">
                <h3 className="font-medium text-sm">{item.product_name}</h3>
                <div className="text-xs text-gray-500">
                  THC: {item.thc_percent}% | CBD: {item.cbd_percent}%
                </div>
              </div>
              <button
                onClick={() => onUpdateCart(items.filter(cartItem => cartItem.product_id !== item.product_id))}
                className="text-red-500 hover:text-red-700 ml-2"
              >
                <Trash2 size={16} />
              </button>
            </div>

            <div className="flex justify-between items-center">
              <div className="flex items-center gap-2">
                <button
                  onClick={() => handleUpdateQuantity(item.product_id, Math.max(0, item.quantity - 1))}
                  className="w-6 h-6 rounded border border-gray-300 flex items-center justify-center hover:bg-gray-50"
                >
                  <Minus size={12} />
                </button>
                <span className="w-8 text-center text-sm">{item.quantity}</span>
                <button
                  onClick={() => handleUpdateQuantity(item.product_id, item.quantity + 1)}
                  className="w-6 h-6 rounded border border-gray-300 flex items-center justify-center hover:bg-gray-50"
                >
                  <Plus size={12} />
                </button>
              </div>
              <div className="font-semibold">${(item.price * item.quantity).toFixed(2)}</div>
            </div>
          </div>
        ))}
      </div>

      <div className="border-t pt-4 space-y-2">
        <div className="flex justify-between text-sm">
          <span>Total THC:</span>
          <span>{totalTHC.toFixed(1)}%</span>
        </div>
        <div className="flex justify-between text-sm">
          <span>Total CBD:</span>
          <span>{totalCBD.toFixed(1)}%</span>
        </div>
        <div className="flex justify-between font-semibold text-lg">
          <span>Total:</span>
          <span>${total.toFixed(2)}</span>
        </div>
      </div>

      <button
        onClick={handleCheckout}
        disabled={isProcessing}
        className="w-full bg-green-600 text-white py-3 rounded-lg font-semibold hover:bg-green-700 disabled:opacity-50 disabled:cursor-not-allowed mt-4 flex items-center justify-center gap-2"
      >
        <CreditCard size={20} />
        {isProcessing ? 'Processing...' : 'Checkout'}
      </button>
    </div>
  )
}