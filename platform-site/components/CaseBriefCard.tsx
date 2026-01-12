import type { CaseBrief } from '@/lib/content';

export default function CaseBriefCard({ brief }: { brief: CaseBrief }) {
  return (
    <article className="relative rounded-3xl border border-sand bg-white/80 p-6 shadow-parchment overflow-hidden">
      <div className="absolute right-0 top-0 h-24 w-24 -translate-y-12 translate-x-6 rotate-12 bg-copper/30" />
      <div className="relative space-y-3">
        <div className="flex flex-wrap gap-2 text-xs font-mono text-indigo/70">
          {brief.tags.map((tag) => (
            <span key={tag} className="rounded-full bg-sand px-3 py-1">
              {tag}
            </span>
          ))}
        </div>
        <h4 className="text-2xl font-display text-charcoal">{brief.title}</h4>
        <p className="font-semibold text-indigo/80">{brief.outcome}</p>
        <p className="text-sm text-indigo/70">{brief.detail}</p>
        <button className="text-sm font-semibold text-teal hover:underline">Download brief →</button>
      </div>
    </article>
  );
}
