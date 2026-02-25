"use client";

import Image from 'next/image';
import { motion } from 'framer-motion';
import HeroHalo from './HeroHalo';
import DataChip from './DataChip';
import JourneyCard from './JourneyStep';
import CaseBriefCard from './CaseBriefCard';
import TerpeneDial from './TerpeneDial';
import FeatureSection from './FeatureSection';
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
      <section className="grid gap-12 lg:gap-16 lg:grid-cols-[1.25fr_0.75fr] pt-16 lg:pt-28 items-start">
        <motion.div variants={stagger} initial="hidden" animate="show" className="space-y-10 max-w-3xl">
          <motion.div variants={fadeUp} className="inline-flex items-center gap-2 rounded-full border border-teal/30 px-4 py-1.5 text-sm text-indigo">
            <span className="h-2 w-2 rounded-full bg-canopy"></span>
            Cross-kingdom OS × Bias-aware evidence engine
          </motion.div>
          <motion.h1 variants={fadeUp} className="text-4xl sm:text-5xl lg:text-6xl font-display leading-tight text-charcoal">
            {heroHeadline.map((line) => (
              <span key={line} className="block">
                {line}
              </span>
            ))}
          </motion.h1>
          <motion.p variants={fadeUp} className="text-lg text-indigo/80 max-w-2xl">
            NeuroBotanica is the world&apos;s first predictive systems biology platform for botanical drug discovery. We don&apos;t just catalog
            molecules; we map their path from ancestral wisdom to FDA compliance across cannabis, marine, fungal, and plant polysaccharide domains—the operating system for cross-kingdom therapeutics.
          </motion.p>
          <motion.div variants={fadeUp} className="flex flex-wrap gap-4">
            <a href="#neurobotanica-video" className="button-ghost">
              See NeuroBotanica in 60s
            </a>
          </motion.div>
          <motion.div variants={fadeUp} id="neurobotanica-video" className="w-full">
            <div className="relative aspect-video w-full rounded-3xl border border-indigo/30 bg-gradient-to-br from-charcoal via-indigo/80 to-canopy/30 overflow-hidden">
              <div className="absolute inset-0 opacity-30 bg-[radial-gradient(circle_at_top,_rgba(67,181,138,0.5),_transparent_65%)]" />
              <div className="relative z-10 flex h-full w-full items-center justify-center text-center px-6">
                <p className="text-white/80 text-base sm:text-lg">
                  60-second overview video placeholder. Drop the final Vimeo or YouTube embed here once the cut is ready and this anchor will auto-play it.
                </p>
              </div>
            </div>
          </motion.div>
          <motion.div variants={fadeUp} className="grid gap-4 sm:grid-cols-3">
            {heroMetrics.map((metric) => (
              <DataChip key={metric.label} metric={metric} />
            ))}
          </motion.div>
        </motion.div>
        <div className="flex flex-col items-center justify-center py-8">
          <Image src="/neurobotanica_logo_1b.png" alt="NeuroBotanica logo" width={220} height={220} className="object-contain mb-2" />
        </div>
      </section>

      <section className="mt-10 rounded-3xl border border-teal/30 bg-gradient-to-r from-white via-cream to-teal/10 px-6 py-6 shadow-parchment flex flex-col gap-4 sm:flex-row sm:items-center sm:justify-between">
        <div className="space-y-2">
          <p className="text-xs font-mono uppercase tracking-[0.4em] text-indigo/60">Nevada pilot intake · Q1</p>
          <h3 className="text-2xl font-display text-charcoal">Two dispensary + lab slots remain for February.</h3>
          <p className="text-sm text-indigo/80 max-w-xl">
            Includes budtender enablement, Workers edge deployment, and evidence packets delivered in under four weeks.
          </p>
        </div>
        <div className="flex flex-col gap-3 sm:items-end">
          <p className="text-sm font-mono text-indigo/70">Average onboarding timeline · <span className="font-semibold text-charcoal">4 weeks</span></p>
          <a
            href="mailto:contessapetrini@cloakquillresearch.org?subject=Pilot%20Briefing%20Request"
            className="button-primary text-center whitespace-nowrap"
          >
            Email for pilot briefing
          </a>
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
              <div className="inline-flex items-center gap-2 text-xs font-mono uppercase tracking-[0.25em] text-indigo/60">
                <span className="h-1.5 w-1.5 rounded-full bg-canopy"></span>
                {card.engines}
              </div>
            </div>
          </motion.article>
        ))}
      </section>

      <section className="mt-24">
        <FeatureSection />
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
              In Phase 1, this section anchors pilots and training; the live tablet UI remains gated while we finalize video
              walkthrough assets for secure distribution.
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
            <div className="rounded-3xl border border-dashed border-canopy/40 bg-white/70 px-5 py-6 text-indigo/80">
              <div className="flex items-center gap-3 text-sm font-semibold text-charcoal mb-4">
                <span className="px-3 py-1 rounded-full bg-canopy/10 text-canopy text-xs tracking-[0.3em] uppercase">Coming Soon</span>
                <p>Workflow demo</p>
              </div>
              <div className="rounded-2xl bg-gradient-to-r from-indigo/10 via-white to-canopy/10 border border-sand px-5 py-4 shadow-inner">
                <p className="text-sm text-indigo/70">
                  Secure training footage will stream here once the build is ready. Request a guided preview to receive the private link.
                </p>
              </div>
            </div>
          </div>
        </div>
      </section>

      <section className="mt-24">
        <div className="flex flex-col gap-2">
          <p className="text-sm uppercase tracking-[0.4em] text-indigo/60">Knowledge Library</p>
          <h3 className="section-heading">Evidence-Based Intelligence.</h3>
          <p className="text-indigo/80 max-w-3xl">
            We don&apos;t just aggregate data; we validate it through peer-reviewed briefs and TK-aligned clinical evidence. Full archive access is provided on request to verified partners.
          </p>
        </div>
        <div className="mt-10 grid gap-6 md:grid-cols-2">
          {caseBriefs.map((brief) => (
            <CaseBriefCard key={brief.title} brief={brief} />
          ))}
        </div>
        {/* Archive access is curated — no public download. */}
      </section>

      {/* Contact form and 'Join the Mission' removed per site guidance. */}

      <footer className="mt-20 border-t border-sand pt-10 pb-6 text-sm text-indigo/70">
        <div className="max-w-6xl mx-auto flex flex-col gap-6 lg:flex-row lg:items-center lg:justify-between">
          <div className="text-center lg:text-left">
            <h4 className="font-semibold">Science in the Public Trust</h4>
            <p className="mt-2 max-w-xl">NeuroBotanica is an initiative of Cloak & Quill Research, a 501(c)(3) nonprofit organization. We are dedicated to advancing evidence-based botanical science without the pressure of shareholder returns.</p>
          </div>
          <div className="flex gap-3 items-center justify-center lg:justify-end">
            <span className="px-3 py-1 rounded-full bg-sand">GuideStar · Platinum</span>
            <span className="px-3 py-1 rounded-full bg-sand">Cloudflare Impact</span>
            <span className="px-3 py-1 rounded-full bg-sand">Stripe Nonprofit</span>
          </div>
        </div>
      </footer>
    </div>
  );
}
