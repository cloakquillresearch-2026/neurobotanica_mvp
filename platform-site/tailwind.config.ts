import type { Config } from 'tailwindcss';

const config: Config = {
  content: [
    './app/**/*.{ts,tsx}',
    './components/**/*.{ts,tsx}',
    './lib/**/*.{ts,tsx}'
  ],
  theme: {
    extend: {
      colors: {
        canopy: 'var(--nb-canopy)',
        teal: 'var(--nb-deep-teal)',
        indigo: 'var(--nb-indigo)',
        sand: 'var(--nb-sand)',
        copper: 'var(--nb-copper)',
        charcoal: 'var(--nb-charcoal)',
        cream: 'var(--nb-cream)'
      },
      fontFamily: {
        display: ['var(--font-display)', 'serif'],
        body: ['var(--font-body)', 'sans-serif'],
        mono: ['var(--font-mono)', 'monospace']
      },
      backgroundImage: {
        'hero-gradient': 'radial-gradient(circle at 10% 20%, rgba(0,124,137,0.3), transparent 60%)',
        'orb-gradient': 'linear-gradient(135deg, #007C89 0%, #43B58A 55%, #A3D86F 100%)'
      },
      animation: {
        'slow-rotate': 'spin 40s linear infinite',
        'float': 'float 8s ease-in-out infinite',
        'fade-up': 'fadeUp 0.8s ease forwards'
      },
      keyframes: {
        float: {
          '0%, 100%': { transform: 'translateY(0px)' },
          '50%': { transform: 'translateY(-12px)' }
        },
        fadeUp: {
          '0%': { opacity: '0', transform: 'translateY(20px)' },
          '100%': { opacity: '1', transform: 'translateY(0)' }
        }
      },
      boxShadow: {
        parchment: '0 25px 60px rgba(0,0,0,0.08)',
        copper: '0 10px 40px rgba(162, 98, 42, 0.35)'
      }
    }
  },
  plugins: []
};

export default config;
