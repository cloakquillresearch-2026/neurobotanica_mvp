/** @type {import('next').NextConfig} */
module.exports = {
  output: 'export',
  trailingSlash: true,
  images: {
    unoptimized: true,
  },
  env: {
    API_URL: process.env.NODE_ENV === 'production'
      ? 'https://api.neurobotanica.com'
      : 'http://localhost:8000',
  },
  // Enable PWA for offline functionality
  experimental: {
    // appDir: true, // Not needed for pages directory
  },
}