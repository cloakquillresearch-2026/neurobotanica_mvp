"use client";

import Image from 'next/image';
import Link from 'next/link';
import { motion } from 'framer-motion';
import HeroHalo from './HeroHalo';
import DataChip from './DataChip';
import JourneyCard from './JourneyStep';
import CaseBriefCard from './CaseBriefCard';
import TerpeneDial from './TerpeneDial';
import {
  heroMetrics,
  platformHighlights,
  journeySteps,
  caseBriefs,
  terpeneMetrics,
  partnerSeals
} from '@/lib/content';

const heroHeadline = ['Predicting Cross-Kingdom Therapeutics,', 'Before They Reach the Lab.'];

const stagger = {
  hidden: {},
  show: {
    transition: {
      staggerChildren: 0.12
    }
  }
};

const fadeUp = {
  hidden: { opacity: 0, y: 20 },
  show: { opacity: 1, y: 0, transition: { duration: 0.8, ease: 'easeOut' } }
};

export default function LandingPage() {
  return (
    <div className="max-w-6xl mx-auto px-6 lg:px-10">
      <section className="grid gap-12 lg:grid-cols-[1.1fr_0.9fr] pt-16 lg:pt-28 items-center">
        <motion.div variants={stagger} initial="hidden" animate="show" className="space-y-8">
          <motion.div variants={fadeUp} className="inline-flex items-center gap-2 rounded-full border border-teal/30 px-4 py-1.5 text-sm text-indigo">
            <span className="h-2 w-2 rounded-full bg-canopy"></span>
            Cross-kingdom AI × Bias-aware evidence engine
          </motion.div>
          <motion.h1 variants={fadeUp} className="text-4xl sm:text-5xl lg:text-6xl font-display leading-tight text-charcoal">
            {heroHeadline.map((line) => (
              <span key={line} className="block">
                {line}
              </span>
            ))}
          </motion.h1>
          <motion.p variants={fadeUp} className="text-lg text-indigo/80 max-w-2xl">
            NeuroBotanica is a cross-kingdom botanical intelligence platform that predicts novel dimeric cannabinoids and polysaccharide
            adjuvant synergies before physical synthesis, organizes evidence into tiers, and helps teams move from ideas to IND-ready
            formulations with clearer traditional knowledge attribution.
          </motion.p>
          <motion.div variants={fadeUp} className="flex flex-wrap gap-4">
            <button className="button-primary">Schedule Lab Session</button>
            <button className="button-ghost">See NeuroBotanica in 60s</button>
            <Link href="/app" className="button-ghost">
              Sign In
            </Link>
          </motion.div>
          <motion.div variants={fadeUp} className="grid gap-4 sm:grid-cols-3">
            {heroMetrics.map((metric) => (
              <DataChip key={metric.label} metric={metric} />
            ))}
          </motion.div>
        </motion.div>
        <div className="relative">
          <div className="paper-shell p-10 bg-cream">
            <HeroHalo />
            <div className="relative z-10 flex flex-col items-center text-center gap-6">
              <Image src="/emblem-mark.svg" alt="NeuroBotanica emblem" width={192} height={192} className="drop-shadow-2xl" />
              <p className="text-indigo/80 max-w-sm">
                DNA helix intertwined with living canopy roots to symbolize rigorous science rooted in ancestral stewardship.
              </p>
              <div className="flex flex-wrap justify-center gap-3 text-sm font-mono text-indigo">
                <span className="px-3 py-1 rounded-full border border-teal/20">Edge @ 285 PoPs</span>
                <span className="px-3 py-1 rounded-full border border-teal/20">HIPAA-ready</span>
                <span className="px-3 py-1 rounded-full border border-teal/20">Stripe Nonprofit</span>
              </div>
            </div>
          </div>
        </div>
      </section>

      <section className="mt-10 rounded-3xl border border-teal/30 bg-gradient-to-r from-white via-cream to-teal/10 px-6 py-6 shadow-parchment flex flex-col gap-4 sm:flex-row sm:items-center sm:justify-between">
        <div className="space-y-2">
          <p className="text-xs font-mono uppercase tracking-[0.4em] text-indigo/60">Nevada pilot intake · Q1</p>
          <h3 className="text-2xl font-display text-charcoal">Two dispensary + lab slots remain for February.</h3>
          <p className="text-sm text-indigo/80 max-w-xl">
            Includes budtender enablement, Workers edge deployment, and HIPAA-ready evidence packets delivered in under four weeks.
          </p>
        </div>
        <div className="flex flex-col gap-3 sm:items-end">
          <p className="text-sm font-mono text-indigo/70">Average onboarding timeline · <span className="font-semibold text-charcoal">4 weeks</span></p>
          <Link href="#contact-form" className="button-primary text-center whitespace-nowrap">
            Reserve pilot briefing
          </Link>
        </div>
      </section>

      <section className="mt-24 bg-white/70 rounded-3xl border border-white/60 shadow-parchment">
        <div className="overflow-hidden">
          <div className="flex gap-12 whitespace-nowrap text-sm font-mono text-indigo/70 animate-[scroll_30s_linear_infinite] px-8 py-5">
            {partnerSeals.concat(partnerSeals).map((seal, idx) => (
              <span key={`${seal}-${idx}`}>{seal}</span>
            ))}
          </div>
        </div>
      </section>

      <section className="mt-24 grid gap-8 lg:grid-cols-2">
        {platformHighlights.map((card) => (
          <motion.article
            key={card.title}
            whileHover={{ y: -8 }}
            className="relative rounded-3xl bg-cream p-8 border border-sand shadow-parchment overflow-hidden"
          >
            <div className="absolute inset-0 opacity-10 bg-[radial-gradient(circle_at_top,_rgba(67,181,138,0.8),_transparent_60%)]" />
            <div className="relative space-y-4">
              <h3 className="text-2xl font-display text-charcoal">{card.title}</h3>
              <p className="text-indigo/80">{card.description}</p>
              <ul className="space-y-2 text-sm text-charcoal/80">
                {card.bullets.map((bullet) => (
                  <li key={bullet} className="flex items-start gap-2">
                    <span className="mt-1 h-2 w-2 rounded-full bg-canopy"></span>
                    {bullet}
                  </li>
                ))}
              </ul>
              <button className="button-ghost">{card.cta}</button>
            </div>
          </motion.article>
        ))}
      </section>

      <section className="mt-24 grid gap-10 lg:grid-cols-[0.7fr_1.3fr]">
        <div className="bg-indigo text-cream rounded-3xl p-8 shadow-copper">
          <p className="uppercase text-sm tracking-[0.2em] text-cream/80">Validation Journey</p>
          <h3 className="section-heading mt-3">From manuscript to deployment in four luminous steps.</h3>
          <p className="mt-4 text-cream/90">
            Scrollytelling timeline pairs copper root systems with living data; ideal for Cloudflare Pages cinematic reveals.
          </p>
        </div>
        <div className="space-y-6">
          {journeySteps.map((step) => (
            <JourneyCard key={step.stage} step={step} />
          ))}
        </div>
      </section>

      <section className="mt-24 rounded-[40px] bg-gradient-to-br from-indigo via-teal to-canopy text-white p-10 shadow-copper overflow-hidden">
        <div className="grid gap-10 lg:grid-cols-[1.1fr_0.9fr] items-center">
          <div>
            <p className="uppercase text-sm tracking-[0.4em] text-white/70">Cannabis Intelligence</p>
            <h3 className="section-heading mt-4 text-white">Nevada-ready terpene insight cockpit.</h3>
            <p className="mt-4 text-white/80 max-w-xl">
              Heatmaps, terpene gauges, and compliance radar unify dispensary POS exhaust with NeuroBotanica inference. Motion cues convey living data without overwhelming the user.
            </p>
            <div className="mt-8 grid gap-4 sm:grid-cols-3">
              {terpeneMetrics.map((metric) => (
                <div key={metric.title} className="bg-white/15 rounded-2xl p-4 backdrop-blur border border-white/20">
                  <div className="text-sm uppercase tracking-[0.2em] text-white/70">{metric.title}</div>
                  <div className="text-2xl font-display mt-2">{metric.delta}</div>
                  <p className="text-sm text-white/80">{metric.description}</p>
                </div>
              ))}
            </div>
          </div>
          <TerpeneDial />
        </div>
      </section>

      <section id="budtender-suite" className="mt-24 rounded-[36px] bg-cream border border-sand shadow-parchment px-8 py-10">
        <div className="grid gap-8 lg:grid-cols-[1.1fr_0.9fr] items-center">
          <div>
            <p className="text-sm uppercase tracking-[0.4em] text-indigo/60">Budtender Suite</p>
            <h3 className="section-heading mt-3">Tablet-first cannabis intelligence for Nevada dispensaries.</h3>
            <p className="mt-4 text-indigo/80 max-w-xl">
              The Budtender Suite is a dedicated, shared-device interface for floor staff. It surfaces strain guidance, terpene
              explanations, and patient-friendly talking points without exposing backend controls or PHI.
            </p>
            <p className="mt-3 text-sm text-indigo/70">
              In Phase 1, this section anchors pilots and training; the live tablet UI is canonically served from the
              <span className="font-mono"> budtender.neuro-botanica.com</span> subdomain.
            </p>
          </div>
          <div className="space-y-4">
            <div className="rounded-3xl bg-white/80 border border-sand px-5 py-4 text-sm text-indigo/80">
              <p className="font-semibold text-charcoal mb-1">Floor-ready workflows</p>
              <p>
                Rapid strain lookup, terpene-centric guidance, and simple patient-facing narratives tailored to Nevada&apos;s regulatory
                context.
              </p>
            </div>
            <a
              href="https://budtender.neuro-botanica.com"
              className="inline-flex items-center justify-center rounded-full bg-canopy px-5 py-2.5 text-sm font-medium text-white shadow-copper hover:bg-canopy/90"
            >
              Open Budtender Suite (tablet)
            </a>
            <p className="text-xs text-indigo/70 font-mono">
              Note: This link reflects the Phase 1 domain plan; DNS and deployment can be staged separately.
            </p>
          </div>
        </div>
      </section>

      <section className="mt-24">
        <div className="flex flex-col gap-4 md:flex-row md:items-end md:justify-between">
          <div>
            <p className="text-sm uppercase tracking-[0.4em] text-indigo/60">Knowledge Library</p>
            <h3 className="section-heading">Case briefs with clipped parchment corners.</h3>
          </div>
          <button className="button-ghost">Browse full archive</button>
        </div>
        <div className="mt-10 grid gap-6 md:grid-cols-2">
          {caseBriefs.map((brief) => (
            <CaseBriefCard key={brief.title} brief={brief} />
          ))}
        </div>
      </section>

      <section id="contact-form" className="mt-24 rounded-[36px] border border-sand bg-cream px-8 py-12 shadow-parchment relative overflow-hidden">
        <div className="absolute inset-0 opacity-30 pointer-events-none">
          <HeroHalo />
        </div>
        <div className="relative grid gap-10 lg:grid-cols-[1.1fr_0.9fr] items-center">
          <div>
            <p className="text-sm uppercase tracking-[0.5em] text-indigo/60">Join the Stewardship Lab</p>
            <h3 className="section-heading mt-4">Book a TradTech session or request a regulated pilot.</h3>
            <p className="mt-4 text-indigo/80 max-w-xl">
              Nonprofit-friendly retainers, Stripe discounts, and cultural equity clauses baked into every engagement.
            </p>
            <div className="mt-6 flex flex-wrap gap-3 text-sm font-mono text-indigo">
              <span className="px-3 py-1 rounded-full bg-sand">GCP + Firebase secured</span>
              <span className="px-3 py-1 rounded-full bg-sand">Cloudflare Workers Edge</span>
              <span className="px-3 py-1 rounded-full bg-sand">Community Revenue Share</span>
            </div>
          </div>
          <form className="space-y-4">
            <div>
              <label className="text-sm font-semibold text-indigo/70">Work Email</label>
              <input type="email" placeholder="you@research.org" className="mt-1 w-full rounded-2xl border border-sand/80 bg-white/70 px-4 py-3 focus:outline-none focus:ring-2 focus:ring-teal" />
            </div>
            <div>
              <label className="text-sm font-semibold text-indigo/70">Organization Type</label>
              <select className="mt-1 w-full rounded-2xl border border-sand/80 bg-white/70 px-4 py-3 focus:outline-none focus:ring-2 focus:ring-teal">
                <option>Nutraceutical / Pharma</option>
                <option>Dispensary</option>
                <option>Indigenous Council</option>
                <option>Regulator / Policy</option>
              </select>
            </div>
            <div>
              <label className="text-sm font-semibold text-indigo/70">Desired Pilot Start</label>
              <input type="date" className="mt-1 w-full rounded-2xl border border-sand/80 bg-white/70 px-4 py-3 focus:outline-none focus:ring-2 focus:ring-teal" />
            </div>
            <button type="submit" className="button-primary w-full text-center">Request Strategy Call</button>
          </form>
        </div>
      </section>

      <footer className="mt-20 border-t border-sand pt-10 pb-6 text-sm text-indigo/70">
        <div className="flex flex-col gap-4 lg:flex-row lg:items-center lg:justify-between">
          <p>Cloak and Quill Research 501(c)(3) · Henderson, Nevada</p>
          <div className="flex gap-4 flex-wrap">
            <a href="#">Patent Claims</a>
            <a href="#">Community Governance</a>
            <a href="#">HIPAA Framework</a>
            <a href="#">Contact</a>
          </div>
        </div>
      </footer>
    </div>
  );
}
