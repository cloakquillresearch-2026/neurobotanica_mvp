function NeuroBotanicaLogo({ className = "" }: { className?: string }) {
  return (
    <div className={`inline-flex flex-col items-center ${className}`}>
      <svg width="120" height="120" viewBox="0 0 120 120" fill="none" xmlns="http://www.w3.org/2000/svg">
        {/* Botanical wreath circle */}
        <circle cx="60" cy="60" r="55" stroke="url(#gradient)" strokeWidth="3" fill="none" />

        {/* Leaf elements around circle */}
        <path d="M60 10 Q50 15 50 25 Q50 20 60 15 Z" fill="#2D9B8B" opacity="0.7" />
        <path d="M90 25 Q85 30 80 35 Q88 32 92 28 Z" fill="#4CAF50" opacity="0.7" />
        <path d="M110 60 Q105 65 100 65 Q108 63 112 60 Z" fill="#66BB6A" opacity="0.7" />
        <path d="M90 95 Q85 90 80 85 Q88 88 92 92 Z" fill="#2D9B8B" opacity="0.7" />
        <path d="M60 110 Q50 105 50 95 Q50 100 60 105 Z" fill="#4CAF50" opacity="0.7" />
        <path d="M30 95 Q35 90 40 85 Q32 88 28 92 Z" fill="#66BB6A" opacity="0.7" />
        <path d="M10 60 Q15 65 20 65 Q12 63 8 60 Z" fill="#2D9B8B" opacity="0.7" />
        <path d="M30 25 Q35 30 40 35 Q32 32 28 28 Z" fill="#4CAF50" opacity="0.7" />

        {/* Letter N in center */}
        <text x="60" y="75" fontSize="48" fontWeight="bold" fill="url(#gradient)" textAnchor="middle" fontFamily="serif">
          N
        </text>

        {/* Gradient definition */}
        <defs>
          <linearGradient id="gradient" x1="0%" y1="0%" x2="100%" y2="100%">
            <stop offset="0%" stopColor="#2D9B8B" />
            <stop offset="50%" stopColor="#4CAF50" />
            <stop offset="100%" stopColor="#66BB6A" />
          </linearGradient>
        </defs>
      </svg>

      <span className="text-xl font-serif text-slate-700 mt-2 tracking-wider">
        NEUROBOTANICA
      </span>
    </div>
  );
}

export { NeuroBotanicaLogo }