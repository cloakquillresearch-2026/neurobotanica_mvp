import { test, expect } from '@playwright/test'

const biomarkerEntries = [
  { label: 'TNF-α', value: '12.5' },
  { label: 'IL-6', value: '8.2' },
  { label: 'CRP', value: '3.4' },
  { label: 'IL-1β', value: '5.1' },
]

test.describe('TS-PS-001 synergy panel', () => {
  test('updates when biomarkers and conditions change', async ({ page }) => {
    await page.goto('/')

    // Start a consultation with an inflammatory focus
    await page.getByRole('button', { name: /Inflammation/i }).first().click()
    await page.getByRole('heading', { name: 'Customer Needs' }).waitFor()

    // Populate biomarker panel to unlock TS-PS-001 personalization
    for (const entry of biomarkerEntries) {
      await page.getByLabel(entry.label).fill(entry.value)
    }

    const saveButton = page.getByRole('button', { name: /NeuroBotanica Recommends/i })
    await saveButton.click()

    const synergyPanel = page.getByTestId('synergy-panel')
    await expect(synergyPanel).toBeVisible()

    const synergyScore = synergyPanel.getByTestId('synergy-score')
    await expect(synergyScore).toHaveText(/%/)
    const baselineScore = (await synergyScore.textContent())?.trim() || ''

    // Change biomarkers and conditions to trigger a fresh TS-PS-001 computation
    await page.getByLabel('CRP').fill('6.0')
    await page.getByRole('button', { name: /Anxiety/i }).click()
    await saveButton.click()

    // Wait for the panel to show its loading state and resolve with new insights
    await page.waitForSelector('[role="status"]:has-text("Computing")', { state: 'visible' })
    await page.waitForSelector('[role="status"]:has-text("Computing")', { state: 'detached' })

    await page.waitForFunction((previousScore: string) => {
      const node = document.querySelector('[data-testid="synergy-score"]')
      return Boolean(node && node.textContent && node.textContent.trim() !== previousScore)
    }, baselineScore)

    await expect(synergyPanel.getByText('Primary Kingdom')).toBeVisible()
    await expect(synergyPanel.getByText('Key Compounds')).toBeVisible()
  })
})
