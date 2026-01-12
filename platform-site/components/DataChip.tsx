import type { HeroMetric } from '@/lib/content';

type Props = {
  metric: HeroMetric;
};

export default function DataChip({ metric }: Props) {
  return (
    <div className="rounded-2xl border border-sand bg-white/70 p-4 shadow-inner">
      <p className="text-sm uppercase tracking-[0.3em] text-indigo/60">{metric.label}</p>
      <p className="mt-2 text-2xl font-display text-charcoal">{metric.value}</p>
      <p className="text-sm text-indigo/80">{metric.detail}</p>
    </div>
  );
}
