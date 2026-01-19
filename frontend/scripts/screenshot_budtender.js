const { chromium } = require('playwright');

(async () => {
  const browser = await chromium.launch({ headless: true });
  const page = await browser.newPage();
  await page.goto('http://127.0.0.1:3000', { waitUntil: 'networkidle' });

  // Click the Quick Start 'Anxiety' button
  const anxietyButton = page.getByRole('button', { name: /Anxiety/i });
  await anxietyButton.click();

  // Wait for Evidence Summary or Key Citations to appear
  await page.waitForSelector('text=Evidence Summary', { timeout: 5000 }).catch(() => {});
  await page.waitForTimeout(1000);

  // Take full page screenshot
  await page.screenshot({ path: 'runs/budtender_anxiety.png', fullPage: true });

  console.log('Screenshot saved to runs/budtender_anxiety.png');
  await browser.close();
})();