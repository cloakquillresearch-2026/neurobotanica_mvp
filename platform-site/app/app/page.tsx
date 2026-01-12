"use client";

import { FormEvent, useState } from 'react';
import Link from 'next/link';
import { signInWithEmailAndPassword } from 'firebase/auth';
import { auth } from '@/lib/firebaseClient';
import { useAuth } from '@/components/AuthProvider';

export default function AppPage() {
  const { user, loading } = useAuth();
  const [email, setEmail] = useState('');
  const [password, setPassword] = useState('');
  const [submitting, setSubmitting] = useState(false);
  const [error, setError] = useState<string | null>(null);

  async function handleEmailSignIn(e: FormEvent<HTMLFormElement>) {
    e.preventDefault();
    setError(null);
    setSubmitting(true);
    try {
      const authInstance = auth;
      if (!authInstance) {
        setError('Authentication is not available. Please try again later.');
        setSubmitting(false);
        return;
      }
      await signInWithEmailAndPassword(authInstance, email.trim(), password);
    } catch (err) {
      const error = err as { code?: string };
      const message = error.code === 'auth/invalid-credential' || error.code === 'auth/wrong-password'
        ? 'Sign-in failed. Please check your email and password.'
        : 'Unable to sign in right now. Please try again.';
      setError(message);
    } finally {
      setSubmitting(false);
    }
  }

  return (
    <main className="min-h-screen bg-sand pb-16">
      <div className="max-w-5xl mx-auto px-6 pt-16 lg:pt-20">
        <div className="mb-4 rounded-2xl border border-sand bg-white/70 px-4 py-3 text-xs sm:text-sm text-indigo/80 font-mono">
          {loading
            ? 'Checking your NeuroBotanica session…'
            : user
            ? `Signed in as ${user.email ?? user.uid}. Experimental workbench surfaces will respect your organization and role context as they are wired in.`
            : 'You are not signed in. Use the form below to sign in with your NeuroBotanica email and password. This workbench will stay read-only until additional tools are enabled.'}
        </div>
        <div className="mb-8 rounded-3xl border border-sand bg-white/80 px-4 py-4 sm:px-6 sm:py-5">
          {user ? (
            <p className="text-xs sm:text-sm text-indigo/80 font-mono">
              You are currently signed in. Future analysis, billing, and pilot surfaces will use this session; for now, the cards below
              remain placeholders.
            </p>
          ) : (
            <form onSubmit={handleEmailSignIn} className="space-y-3">
              <div className="flex flex-col gap-2 sm:flex-row sm:items-center sm:gap-3">
                <input
                  type="email"
                  required
                  value={email}
                  onChange={(e) => setEmail(e.target.value)}
                  placeholder="you@clinic.org"
                  className="flex-1 rounded-2xl border border-sand/80 bg-white/90 px-3 py-2 text-sm focus:outline-none focus:ring-2 focus:ring-teal"
                />
                <input
                  type="password"
                  required
                  value={password}
                  onChange={(e) => setPassword(e.target.value)}
                  placeholder="Password"
                  className="flex-1 rounded-2xl border border-sand/80 bg-white/90 px-3 py-2 text-sm focus:outline-none focus:ring-2 focus:ring-teal"
                />
                <button
                  type="submit"
                  className="mt-1 inline-flex items-center justify-center rounded-full bg-canopy px-4 py-2 text-sm font-medium text-white shadow-copper hover:bg-canopy/90 disabled:opacity-60"
                  disabled={submitting}
                >
                  {submitting ? 'Signing in…' : 'Sign in'}
                </button>
              </div>
              {error && <p className="text-xs text-red-600 font-mono">{error}</p>}
              <p className="text-[11px] sm:text-xs text-indigo/70 font-mono">
                Uses the same Firebase Authentication project as the NeuroBotanica budtender and API flows. Account creation and role
                assignment are still handled out-of-band during pilots.
              </p>
            </form>
          )}
        </div>
        <header className="mb-10">
          <p className="text-sm uppercase tracking-[0.35em] text-indigo/60 mb-3">NeuroBotanica Workbench</p>
          <h1 className="text-3xl sm:text-4xl lg:text-5xl font-display text-charcoal leading-tight">
            NeuroBotanica Workbench
          </h1>
          <p className="mt-4 text-indigo/80 max-w-2xl text-base sm:text-lg">
            Phase 1 shell for authenticated organizations. Pilot status, analysis tools, and dispensary surfaces will plug in here as
            subscription flows and Firebase sign-in are wired.
          </p>
        </header>

        <section className="grid gap-6 md:grid-cols-3">
          <div className="rounded-3xl bg-cream border border-sand shadow-parchment p-6 flex flex-col justify-between">
            <div>
              <h2 className="text-lg font-display text-charcoal mb-2">Pilot Status</h2>
              <p className="text-sm text-indigo/80">
                High-level snapshot of your active NeuroBotanica or VeriTrad pilots. In Phase 1 this is a static summary stub; dynamic
                status wiring will follow Stripe + org setup.
              </p>
            </div>
            <p className="mt-4 text-xs font-mono text-indigo/70">Status feed coming soon · Phase 1 placeholder</p>
          </div>

          <div className="rounded-3xl bg-cream border border-sand shadow-parchment p-6 flex flex-col justify-between">
            <div>
              <h2 className="text-lg font-display text-charcoal mb-2">Run a Test Analysis</h2>
              <p className="text-sm text-indigo/80">
                Entry point for routing a single strain or dossier into the discovery API. During Phase 1 this is a non-interactive stub
                that will later open guided analysis flows.
              </p>
            </div>
            <button
              type="button"
              className="mt-4 inline-flex items-center justify-center rounded-full border border-teal/40 px-4 py-2 text-sm font-medium text-teal/80 bg-white/70 cursor-not-allowed"
            >
              Test analysis launcher (coming soon)
            </button>
          </div>

          <div className="rounded-3xl bg-cream border border-sand shadow-parchment p-6 flex flex-col justify-between">
            <div>
              <h2 className="text-lg font-display text-charcoal mb-2">Dispensary Tools</h2>
              <p className="text-sm text-indigo/80">
                Overview of the Budtender Suite for Nevada dispensaries. This links to informational content about the floor interface
                rather than the live shared-device UI.
              </p>
            </div>
            <Link
              href="/#budtender-suite"
              className="mt-4 inline-flex items-center justify-center rounded-full bg-canopy px-4 py-2 text-sm font-medium text-white shadow-copper hover:bg-canopy/90"
            >
              Learn about the Budtender Suite
            </Link>
          </div>
        </section>

        <section className="mt-10 rounded-3xl border border-sand bg-white/70 px-6 py-5 text-xs text-indigo/70 font-mono">
          <p>
            Auth note: this Phase 1 shell does not yet enforce Firebase Authentication or Stripe-backed org routing. Treat this as a
            visual contract only until the identity and billing flows from the SaaS blueprint are implemented.
          </p>
        </section>
      </div>
    </main>
  );
}
