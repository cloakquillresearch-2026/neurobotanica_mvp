"use client";

import { motion } from 'framer-motion';
import type { JourneyStep as JourneyType } from '@/lib/content';

type Props = {
  step: JourneyType;
};

export default function JourneyCard({ step }: Props) {
  return (
    <motion.article
      whileInView={{ opacity: [0, 1], y: [20, 0] }}
      viewport={{ once: true, amount: 0.4 }}
      transition={{ duration: 0.7, ease: 'easeOut' }}
      className="relative overflow-hidden rounded-3xl border border-sand bg-white/80 p-6"
    >
      <div className="absolute left-0 top-0 h-full w-1 bg-gradient-to-b from-copper to-canopy" />
      <div className="ml-6 space-y-2">
        <p className="font-mono text-sm text-indigo/60">{step.stage}</p>
        <h4 className="text-2xl font-display text-charcoal">{step.title}</h4>
        <p className="text-indigo/80">{step.description}</p>
        <p className="text-sm font-mono text-indigo/70">{step.kpi}</p>
      </div>
    </motion.article>
  );
}
