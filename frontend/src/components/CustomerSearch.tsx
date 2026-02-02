import { useState } from 'react'
import { dispensaryAPI } from '@/utils/api'
import type { CustomerProfileData } from '@/types/customer'

interface CustomerSearchProps {
  onCustomerSelect: (customer: CustomerProfileData) => void
  isSandboxMode?: boolean
}

const COMMON_CONDITIONS = [
  { id: 'chronic_pain', label: 'Chronic Pain', emoji: '🔥', color: 'from-orange-500 to-red-500' },
  { id: 'anxiety', label: 'Anxiety', emoji: '🧘', color: 'from-blue-500 to-cyan-500' },
  { id: 'insomnia', label: 'Insomnia', emoji: '🌙', color: 'from-purple-500 to-indigo-500' },
  { id: 'ptsd', label: 'PTSD', emoji: '💚', color: 'from-teal-500 to-emerald-500' },
  { id: 'nausea', label: 'Nausea', emoji: '🤢', color: 'from-green-500 to-lime-500' },
  { id: 'inflammation', label: 'Inflammation', emoji: '🔴', color: 'from-rose-500 to-pink-500' },
  { id: 'depression', label: 'Depression', emoji: '☀️', color: 'from-amber-500 to-yellow-500' },
  { id: 'weight_management', label: 'Metabolism', emoji: '⚖️', color: 'from-cyan-500 to-blue-500' },
]

export function CustomerSearch({ onCustomerSelect, isSandboxMode = false }: CustomerSearchProps) {
  const [searchTerm, setSearchTerm] = useState('')
  const [searchResults, setSearchResults] = useState<CustomerProfileData[]>([])
  const [isSearching, setIsSearching] = useState(false)
  const [showQuickStart, setShowQuickStart] = useState(true)
  const [showNewClientModal, setShowNewClientModal] = useState(false)
  const [newClientName, setNewClientName] = useState({ first: '', last: '' })

  const handleSearch = async (term: string) => {
    if (!term.trim()) return

    setIsSearching(true)
    setShowQuickStart(false)
    try {
      const response = await dispensaryAPI.searchProfiles(term)
      const data = response.data
      setSearchResults(data.customers || [])
    } catch (error) {
      console.error('Search failed:', error)
      // If API fails, show empty results instead of mock data
      setSearchResults([])
    } finally {
      setIsSearching(false)
    }
  }

  const handleConditionSelect = (conditionId: string) => {
    // This should only update the customer's conditions during consultation
    // For now, just log - this will be handled by the parent component
    console.log('Condition selected:', conditionId)
  }

  const handleNewClient = () => {
    setShowNewClientModal(true)
  }

  const handleSaveNewClient = async () => {
    if (!newClientName.first.trim() || !newClientName.last.trim()) return

    try {
      console.log('Save Client clicked (new client modal)')
      const newCustomer: CustomerProfileData = {
        customer_id: `temp_${Date.now()}`,
        first_name: newClientName.first.trim(),
        last_name: newClientName.last.trim(),
        phone: '',
        conditions: [],
        experience_level: 'beginner',
        age: undefined,
        gender: '',
        weight: undefined,
        notes: '',
        biomarkers: {},
        isNew: true,
        isSandbox: false
      }

      // Save to database first
      const response = await dispensaryAPI.createProfile({
        first_name: newCustomer.first_name,
        last_name: newCustomer.last_name,
        age: newCustomer.age || null,
        biological_sex: newCustomer.gender || '',
        weight_kg: newCustomer.weight || null,
        conditions: [],
        experience_level: newCustomer.experience_level,
        administration_preferences: [],
        primary_goal: '',
        biomarkers: Object.fromEntries(
          Object.entries(newCustomer.biomarkers || {}).filter(([, value]) => value !== undefined)
        ) as Record<string, number>,
        notes: newCustomer.notes
      })

      console.log('Create profile response:', response?.data)

      // Update the customer ID with the one from the database
      newCustomer.customer_id = response.data.customer_id
      newCustomer.isNew = false

      // Close modal and reset form
      setShowNewClientModal(false)
      setNewClientName({ first: '', last: '' })

      // Select the customer to show consultation view
      onCustomerSelect(newCustomer)
      setShowQuickStart(false)
    } catch (error) {
      console.error('Failed to create new client:', error)
      // Show error to user instead of creating local customer
      alert('Failed to save client to database. Please try again.')
    }
  }

  return (
    <div className="space-y-5">
      {/* Search Bar */}
      <div className="relative">
        <input
          type="text"
          placeholder="🔍 Search by phone or name..."
          value={searchTerm}
          onChange={(e) => setSearchTerm(e.target.value)}
          onKeyPress={(e) => e.key === 'Enter' && handleSearch(searchTerm)}
          className="input-modern pr-12"
        />
        <button
          onClick={() => handleSearch(searchTerm)}
          disabled={isSearching}
          className="absolute right-2 top-2 w-9 h-9 bg-gradient-to-r from-emerald-500 to-teal-500 rounded-lg flex items-center justify-center text-white hover:from-emerald-400 hover:to-teal-400 transition-all"
        >
          {isSearching ? (
            <div className="w-4 h-4 border-2 border-white border-t-transparent rounded-full animate-spin" />
          ) : (
            <svg className="w-5 h-5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M21 21l-6-6m2-5a7 7 0 11-14 0 7 7 0 0114 0z" />
            </svg>
          )}
        </button>
      </div>

      {/* Quick Start - Condition Selection */}
      {showQuickStart && (
        <div className="bg-gradient-to-r from-emerald-500/10 to-teal-500/10 rounded-xl p-5 border border-emerald-400/30">
          <div className="flex items-center gap-2 mb-4">
            <span className="text-xl">⚡</span>
            <span className="font-bold text-white">Quick Start: What brings them in today?</span>
          </div>
          
          {/* New Client Button - Hide in sandbox mode */}
          {!isSandboxMode && (
            <button
              onClick={handleNewClient}
              className="w-full mb-4 py-3 bg-gradient-to-r from-purple-500 to-indigo-500 hover:from-purple-400 hover:to-indigo-400 text-white font-semibold rounded-xl transition-all border border-purple-400/30 shadow-lg"
            >
              <div className="flex items-center justify-center gap-2">
                <svg className="w-5 h-5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                  <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M18 9v3m0 0v3m0-3h3m-3 0h-3m-2-5a4 4 0 11-8 0 4 4 0 018 0zM3 20a6 6 0 0112 0v1H3v-1z" />
                </svg>
                + New Client
              </div>
            </button>
          )}
          
          <div className="grid grid-cols-2 sm:grid-cols-4 gap-3">
            {COMMON_CONDITIONS.map((condition) => (
              <button
                key={condition.id}
                onClick={() => handleConditionSelect(condition.id)}
                className={`group relative overflow-hidden bg-gradient-to-r ${condition.color} p-0.5 rounded-xl hover:scale-105 transition-all duration-200 shadow-lg`}
              >
                <div className="bg-gray-900/90 rounded-[10px] p-3 h-full">
                  <div className="text-2xl mb-1">{condition.emoji}</div>
                  <div className="text-white text-sm font-medium">{condition.label}</div>
                </div>
              </button>
            ))}
          </div>
          <button
            onClick={() => handleConditionSelect('general')}
            className="w-full mt-4 py-3 text-emerald-400 hover:text-white font-semibold bg-white/5 hover:bg-white/10 rounded-xl transition-all border border-white/10"
          >
            Start General Consultation →
          </button>
        </div>
      )}

      {/* Search Results */}
      {searchResults.length > 0 && (
        <div className="space-y-3">
          <div className="flex items-center gap-2 text-white/60 text-sm">
            <svg className="w-4 h-4" fill="currentColor" viewBox="0 0 20 20">
              <path d="M9 6a3 3 0 11-6 0 3 3 0 016 0zM17 6a3 3 0 11-6 0 3 3 0 016 0zM12.93 17c.046-.327.07-.66.07-1a6.97 6.97 0 00-1.5-4.33A5 5 0 0119 16v1h-6.07zM6 11a5 5 0 015 5v1H1v-1a5 5 0 015-5z"/>
            </svg>
            <span>Returning Customers</span>
          </div>
          {searchResults.map((customer) => (
            <button
              key={customer.customer_id}
              onClick={() => onCustomerSelect(customer)}
              className="w-full text-left p-4 bg-white/90 rounded-xl hover:shadow-xl hover:shadow-emerald-500/20 transition-all group"
            >
              <div className="flex justify-between items-start">
                <div>
                  <div className="font-bold text-gray-800 group-hover:text-emerald-700 transition-colors">
                    {customer.first_name} {customer.last_name}
                  </div>
                  <div className="text-gray-500 text-sm">{customer.phone}</div>
                </div>
                <div className="text-xs text-gray-400 bg-gray-100 px-2 py-1 rounded-full">
                  Last: {customer.last_visit}
                </div>
              </div>
              {customer.conditions && customer.conditions.length > 0 && (
                <div className="flex flex-wrap gap-2 mt-3">
                  {customer.conditions.map((condition: string, index: number) => {
                    const condInfo = COMMON_CONDITIONS.find(c => c.id === condition)
                    return (
                      <span
                        key={index}
                        className={`text-xs px-2 py-1 rounded-full bg-gradient-to-r ${condInfo?.color || 'from-gray-400 to-gray-500'} text-white`}
                      >
                        {condInfo?.emoji} {condInfo?.label || condition}
                      </span>
                    )
                  })}
                </div>
              )}
            </button>
          ))}
        </div>
      )}

      {/* New Customer Button */}
      {searchTerm && searchResults.length === 0 && !isSearching && !showQuickStart && (
        <button
          onClick={() => handleConditionSelect('new')}
          className="btn-primary w-full flex items-center justify-center gap-2"
        >
          <svg className="w-5 h-5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M18 9v3m0 0v3m0-3h3m-3 0h-3m-2-5a4 4 0 11-8 0 4 4 0 018 0zM3 20a6 6 0 0112 0v1H3v-1z" />
          </svg>
          Start New Consultation
        </button>
      )}

      {/* New Client Modal */}
      {showNewClientModal && (
        <div className="fixed inset-0 bg-black/50 backdrop-blur-sm flex items-center justify-center z-50 p-4">
          <div className="bg-gray-900 rounded-xl p-6 w-full max-w-md border border-white/10">
            <h3 className="text-xl font-bold text-white mb-4">New Client</h3>
            
            <div className="space-y-4">
              <div>
                <label className="block text-white/80 text-sm font-medium mb-2">
                  First Name *
                </label>
                <input
                  type="text"
                  value={newClientName.first}
                  onChange={(e) => setNewClientName(prev => ({ ...prev, first: e.target.value }))}
                  className="w-full px-3 py-2 bg-white/10 border border-white/20 rounded-lg text-white placeholder-white/50 focus:outline-none focus:border-emerald-400"
                  placeholder="Enter first name"
                  autoFocus
                />
              </div>
              
              <div>
                <label className="block text-white/80 text-sm font-medium mb-2">
                  Last Name *
                </label>
                <input
                  type="text"
                  value={newClientName.last}
                  onChange={(e) => setNewClientName(prev => ({ ...prev, last: e.target.value }))}
                  className="w-full px-3 py-2 bg-white/10 border border-white/20 rounded-lg text-white placeholder-white/50 focus:outline-none focus:border-emerald-400"
                  placeholder="Enter last name"
                  onKeyPress={(e) => e.key === 'Enter' && handleSaveNewClient()}
                />
              </div>
            </div>
            
            <div className="flex gap-3 mt-6">
              <button
                onClick={() => {
                  setShowNewClientModal(false)
                  setNewClientName({ first: '', last: '' })
                }}
                className="flex-1 py-2 px-4 bg-white/10 hover:bg-white/20 text-white/80 rounded-lg transition-colors"
              >
                Cancel
              </button>
              <button
                onClick={handleSaveNewClient}
                disabled={!newClientName.first.trim() || !newClientName.last.trim()}
                className="flex-1 py-2 px-4 bg-gradient-to-r from-emerald-500 to-teal-500 hover:from-emerald-400 hover:to-teal-400 text-white font-semibold rounded-lg transition-all disabled:opacity-50 disabled:cursor-not-allowed"
              >
                Save Client
              </button>
            </div>
          </div>
        </div>
      )}
    </div>
  )
}
