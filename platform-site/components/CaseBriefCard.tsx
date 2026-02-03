import type { CaseBrief } from '@/lib/content';

export default function CaseBriefCard({ brief }: { brief: CaseBrief }) {
  return (
    <article className="relative rounded-3xl border border-sand bg-white/80 p-6 shadow-parchment overflow-hidden">
      <div className="absolute right-0 top-0 h-24 w-24 -translate-y-12 translate-x-6 rotate-12 bg-copper/30" />
      <div className="relative space-y-3">
        <p className="text-xs font-mono uppercase tracking-[0.3em] text-indigo/60">Knowledge Library</p>
        <h4 className="text-2xl font-display text-charcoal">{brief.title}</h4>
        <p className="text-sm text-indigo/70">{brief.subtitle}</p>
        <div className="space-y-2 text-sm text-indigo/80">
          <p><span className="font-semibold text-charcoal">Topic:</span> {brief.topic}</p>
          <p><span className="font-semibold text-charcoal">Key insight:</span> {brief.keyInsight}</p>
          <p><span className="font-semibold text-charcoal">Our solution:</span> {brief.solution}</p>
        </div>
        <p className="text-xs font-mono text-indigo/60">{brief.citation}</p>
      </div>
    </article>
  );
}
