"use client";

import { motion } from 'framer-motion';

export default function TerpeneDial() {
  return (
    <motion.div
      initial={{ opacity: 0, y: 30 }}
      whileInView={{ opacity: 1, y: 0 }}
      viewport={{ once: true, amount: 0.4 }}
      transition={{ duration: 0.8, ease: 'easeOut' }}
      className="relative flex items-center justify-center"
    >
      <div className="relative h-72 w-72 rounded-full bg-white/10 backdrop-blur border border-white/30 p-6">
        <div
          className="absolute inset-4 rounded-full"
          style={{
            background: 'conic-gradient(from 180deg, rgba(255,255,255,0.15), rgba(255,255,255,0.85) 240deg, rgba(255,255,255,0.15))'
          }}
        />
        <div className="absolute inset-8 rounded-full border border-white/20"></div>
        <div className="absolute inset-12 rounded-full bg-white/10"></div>
        <div className="absolute left-1/2 top-6 h-8 w-1 bg-white rounded-full origin-bottom -translate-x-1/2" />
        <div className="absolute inset-0 flex items-center justify-center">
          <div className="text-center">
            <p className="text-sm uppercase tracking-[0.4em] text-white/60">Optimization</p>
            <p className="text-4xl font-display">84%</p>
            <p className="text-white/80">Projected terpene uplift</p>
          </div>
        </div>
      </div>
    </motion.div>
  );
}
