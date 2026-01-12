export default function HeroHalo() {
  return (
    <svg
      className="absolute inset-0 mx-auto my-auto h-full w-full opacity-70 animate-slow-rotate"
      viewBox="0 0 400 400"
      aria-hidden
    >
      <defs>
        <linearGradient id="halo" x1="0%" y1="0%" x2="100%" y2="100%">
          <stop offset="0%" stopColor="#007C89" stopOpacity={0.2} />
          <stop offset="50%" stopColor="#43B58A" stopOpacity={0.4} />
          <stop offset="100%" stopColor="#A3D86F" stopOpacity={0.2} />
        </linearGradient>
      </defs>
      <circle cx="200" cy="200" r="180" stroke="url(#halo)" strokeWidth="1.5" fill="none" />
      <circle cx="200" cy="200" r="140" stroke="url(#halo)" strokeWidth="1" fill="none" strokeDasharray="6 10" />
      <path
        d="M200 40 C260 140 140 260 200 360"
        stroke="#43B58A"
        strokeWidth="1"
        strokeDasharray="5 12"
        fill="none"
        opacity="0.4"
      />
      <path
        d="M200 360 C140 260 260 140 200 40"
        stroke="#007C89"
        strokeWidth="1"
        strokeDasharray="5 12"
        fill="none"
        opacity="0.4"
      />
    </svg>
  );
}
