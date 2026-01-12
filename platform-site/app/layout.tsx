import type { Metadata } from 'next';
import './globals.css';
import { Fraunces, Space_Grotesk, IBM_Plex_Mono } from 'next/font/google';
import { AuthProvider } from '@/components/AuthProvider';

const fraunces = Fraunces({ subsets: ['latin'], weight: ['400', '500', '600', '700'], variable: '--font-display' });
const spaceGrotesk = Space_Grotesk({ subsets: ['latin'], weight: ['400', '500'], variable: '--font-body' });
const plexMono = IBM_Plex_Mono({ subsets: ['latin'], weight: ['500'], variable: '--font-mono' });

export const metadata: Metadata = {
  title: 'NeuroBotanica Platform',
  description: 'Validating traditional knowledge and optimizing terpene intelligence with bias-aware, evidence-grade rigor.',
  metadataBase: new URL('https://neurobotanica-mvp.pages.dev'),
  openGraph: {
    title: 'NeuroBotanica Platform',
    description: 'TradTech innovation for botanical therapeutics.',
    url: 'https://neurobotanica-mvp.pages.dev',
    siteName: 'NeuroBotanica',
    images: [
      {
        url: '/emblem-mark.svg',
        width: 1200,
        height: 630,
        alt: 'NeuroBotanica emblem hovering above an illuminated herbarium.'
      }
    ]
  }
};

export default function RootLayout({ children }: { children: React.ReactNode }) {
  return (
    <html lang="en" className={`${fraunces.variable} ${spaceGrotesk.variable} ${plexMono.variable}`}>
      <body className="antialiased bg-sand">
        <AuthProvider>{children}</AuthProvider>
      </body>
    </html>
  );
}
