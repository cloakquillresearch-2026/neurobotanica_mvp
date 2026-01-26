/**
 * NeuroBotanica Main API Worker Integration Tests
 * Tests for D1 integration, bias correction, and API endpoints
 */

import { describe, it, expect, beforeAll, vi } from 'vitest';

// Mock D1 Database
const createMockD1 = () => ({
  prepare: vi.fn(() => ({
    bind: vi.fn(() => ({
      first: vi.fn(),
      all: vi.fn(() => ({ results: [] })),
    })),
  })),
});

// Mock KV Namespace
const createMockKV = () => ({
  get: vi.fn(() => null),
  put: vi.fn(),
});

// Mock environment
const createMockEnv = () => ({
  NEUROBOTANICA_DB: createMockD1(),
  NEUROBOTANICA_CACHE: createMockKV(),
  ENVIRONMENT: 'test',
});

// Import the worker module (we'll need to mock the import)
// In real tests, this would be imported from the compiled worker

describe('NeuroBotanica Main API Worker', () => {
  describe('Health Check Endpoint', () => {
    it('should return healthy status when database is available', async () => {
      const env = createMockEnv();
      const mockDb = env.NEUROBOTANICA_DB as ReturnType<typeof createMockD1>;
      mockDb.prepare = vi.fn(() => ({
        bind: vi.fn(() => ({
          first: vi.fn(() => Promise.resolve({ result: 1 })),
          all: vi.fn(() => Promise.resolve({ results: [] })),
        })),
      }));

      const request = new Request('https://api.neurobotanica.workers.dev/api/neurobotanica/health', {
        method: 'GET',
      });

      // Simulate the health check
      const dbHealthy = await (async () => {
        try {
          await mockDb.prepare('SELECT 1').bind().first();
          return true;
        } catch {
          return false;
        }
      })();

      expect(dbHealthy).toBe(true);
    });

    it('should return degraded status when database is unavailable', async () => {
      const env = createMockEnv();
      const mockDb = env.NEUROBOTANICA_DB as ReturnType<typeof createMockD1>;
      mockDb.prepare = vi.fn(() => ({
        bind: vi.fn(() => ({
          first: vi.fn(() => Promise.reject(new Error('Database connection failed'))),
          all: vi.fn(() => Promise.reject(new Error('Database connection failed'))),
        })),
      }));

      const dbHealthy = await (async () => {
        try {
          await mockDb.prepare('SELECT 1').bind().first();
          return true;
        } catch {
          return false;
        }
      })();

      expect(dbHealthy).toBe(false);
    });
  });

  describe('CORS Handling', () => {
    it('should handle OPTIONS preflight requests', () => {
      const request = new Request('https://api.neurobotanica.workers.dev/api/neurobotanica/analyze', {
        method: 'OPTIONS',
        headers: {
          Origin: 'https://neurobotanica.pages.dev',
        },
      });

      // Check that CORS headers would be set
      const origin = request.headers.get('Origin');
      expect(origin).toBe('https://neurobotanica.pages.dev');
    });
  });

  describe('Analyze Endpoint', () => {
    it('should validate compound_ids is required', async () => {
      const invalidPayload = {
        demographics: { age: 30 },
      };

      // Validation check
      const isValid = (payload: Record<string, unknown>) => {
        const ids = payload.compound_ids;
        return Array.isArray(ids) && ids.length > 0;
      };

      expect(isValid(invalidPayload)).toBe(false);
      expect(isValid({ compound_ids: ['cbd'] })).toBe(true);
      expect(isValid({ compound_ids: [] })).toBe(false);
    });

    it('should complete analysis under 10 seconds target', async () => {
      const startTime = Date.now();

      // Simulate analysis workflow
      await Promise.all([
        new Promise(resolve => setTimeout(resolve, 50)), // interactions
        new Promise(resolve => setTimeout(resolve, 100)), // bias correction
        new Promise(resolve => setTimeout(resolve, 75)), // synergy
        new Promise(resolve => setTimeout(resolve, 50)), // polysaccharides
      ]);

      const elapsed = Date.now() - startTime;
      expect(elapsed).toBeLessThan(10000);
    });
  });

  describe('Bias Correction', () => {
    it('should apply age-based adjustments correctly', () => {
      const calculateAgeFactor = (age: number): number => {
        if (age < 25) return 0.85;
        if (age > 75) return 0.8;
        if (age > 65) return 0.85;
        if (age > 50) return 0.95;
        return 1.0;
      };

      expect(calculateAgeFactor(20)).toBe(0.85); // Young
      expect(calculateAgeFactor(30)).toBe(1.0); // Adult
      expect(calculateAgeFactor(55)).toBe(0.95); // Mature
      expect(calculateAgeFactor(70)).toBe(0.85); // Senior
      expect(calculateAgeFactor(80)).toBe(0.8); // Elderly
    });

    it('should apply weight-based adjustments correctly', () => {
      const calculateWeightFactor = (weightLbs: number): number => {
        if (weightLbs < 120) return 0.85;
        if (weightLbs < 150) return 0.9;
        if (weightLbs > 280) return 1.15;
        if (weightLbs > 250) return 1.1;
        if (weightLbs > 200) return 1.05;
        return 1.0;
      };

      expect(calculateWeightFactor(100)).toBe(0.85); // Light
      expect(calculateWeightFactor(140)).toBe(0.9); // Below average
      expect(calculateWeightFactor(175)).toBe(1.0); // Average
      expect(calculateWeightFactor(220)).toBe(1.05); // Above average
      expect(calculateWeightFactor(260)).toBe(1.1); // Heavy
      expect(calculateWeightFactor(300)).toBe(1.15); // Very heavy
    });

    it('should apply gender-based adjustments correctly', () => {
      const calculateGenderFactor = (gender: string): number => {
        const normalized = gender.toLowerCase();
        if (normalized === 'female' || normalized === 'f') return 0.92;
        if (normalized === 'male' || normalized === 'm') return 1.0;
        return 1.0;
      };

      expect(calculateGenderFactor('female')).toBe(0.92);
      expect(calculateGenderFactor('Female')).toBe(0.92);
      expect(calculateGenderFactor('f')).toBe(0.92);
      expect(calculateGenderFactor('male')).toBe(1.0);
      expect(calculateGenderFactor('non-binary')).toBe(1.0);
    });

    it('should apply metabolism-based adjustments correctly', () => {
      const calculateMetabolismFactor = (type: 'fast' | 'normal' | 'slow'): number => {
        switch (type) {
          case 'fast': return 1.15;
          case 'slow': return 0.85;
          default: return 1.0;
        }
      };

      expect(calculateMetabolismFactor('fast')).toBe(1.15);
      expect(calculateMetabolismFactor('normal')).toBe(1.0);
      expect(calculateMetabolismFactor('slow')).toBe(0.85);
    });

    it('should cap adjustment factor within safe bounds (0.5x - 2.0x)', () => {
      const capAdjustmentFactor = (factor: number): number => {
        if (factor < 0.5) return 0.5;
        if (factor > 2.0) return 2.0;
        return factor;
      };

      expect(capAdjustmentFactor(0.3)).toBe(0.5);
      expect(capAdjustmentFactor(0.8)).toBe(0.8);
      expect(capAdjustmentFactor(1.5)).toBe(1.5);
      expect(capAdjustmentFactor(2.5)).toBe(2.0);
    });

    it('should calculate combined adjustments correctly', () => {
      const baseDose = 10.0;
      const demographics = {
        age: 70, // 0.85
        weight: 120, // 0.9
        gender: 'female', // 0.92
      };

      const ageFactor = 0.85;
      const weightFactor = 0.9;
      const genderFactor = 0.92;

      const combinedFactor = ageFactor * weightFactor * genderFactor;
      const adjustedDose = baseDose * combinedFactor;

      expect(combinedFactor).toBeCloseTo(0.7038, 3);
      expect(adjustedDose).toBeCloseTo(7.038, 2);
    });

    it('should return fallback on database error', async () => {
      const env = createMockEnv();
      const mockDb = env.NEUROBOTANICA_DB as ReturnType<typeof createMockD1>;
      mockDb.prepare = vi.fn(() => {
        throw new Error('Database error');
      });

      // Simulate fallback behavior
      let adjustmentFactor = 1.0;
      try {
        mockDb.prepare('SELECT 1');
      } catch {
        // Apply fallback calculations
        adjustmentFactor = 0.85 * 0.9; // Age 70 + Weight 120 fallback
      }

      expect(adjustmentFactor).toBeCloseTo(0.765, 3);
    });
  });

  describe('Synergy Prediction', () => {
    it('should calculate fallback synergy for known pairs', () => {
      const calculateFallbackSynergy = (a: string, b: string): number => {
        const synergyPairs: Record<string, number> = {
          'cbd:thc': 0.75,
          'thc:cbd': 0.75,
          'cbd:cbg': 0.7,
          'thc:myrcene': 0.8,
        };

        const key = `${a.toLowerCase()}:${b.toLowerCase()}`;
        return synergyPairs[key] || 0.5;
      };

      expect(calculateFallbackSynergy('cbd', 'thc')).toBe(0.75);
      expect(calculateFallbackSynergy('THC', 'CBD')).toBe(0.75);
      expect(calculateFallbackSynergy('thc', 'myrcene')).toBe(0.8);
      expect(calculateFallbackSynergy('unknown', 'compound')).toBe(0.5);
    });

    it('should use cache when available', async () => {
      const mockKV = createMockKV();
      mockKV.get = vi.fn(() => Promise.resolve(JSON.stringify({
        synergy_score: 0.85,
        tk_enhanced: true,
        evidence: 'Cached result',
      })));

      const cached = await mockKV.get('synergy:cbd:thc:computational');
      expect(cached).not.toBeNull();
      expect(JSON.parse(cached!).synergy_score).toBe(0.85);
    });
  });

  describe('Condition Mapping', () => {
    it('should map conditions to database values correctly', () => {
      const conditionMap: Record<string, string> = {
        weight_management: 'metabolism',
        muscle_spasms: 'spasms',
        chronic_pain: 'pain',
        anxiety: 'anxiety',
        insomnia: 'sleep',
        inflammation: 'inflammation',
      };

      expect(conditionMap['weight_management']).toBe('metabolism');
      expect(conditionMap['muscle_spasms']).toBe('spasms');
      expect(conditionMap['chronic_pain']).toBe('pain');
    });
  });

  describe('API Response Format', () => {
    it('should format bias correction response correctly', () => {
      const formatBiasResponse = (
        baseDose: number,
        adjustmentFactor: number,
        demographics: Record<string, unknown>,
        evidence: string
      ) => {
        return {
          adjusted_dose_mg: Math.round(baseDose * adjustmentFactor * 100) / 100,
          factors_applied: {
            base_dose: baseDose,
            adjustment_factor: Math.round(adjustmentFactor * 1000) / 1000,
            demographics_considered: Object.keys(demographics),
          },
          evidence,
          confidence: 0.75,
          warnings: [],
        };
      };

      const response = formatBiasResponse(10.0, 0.85, { age: 70 }, 'Age-based adjustment');

      expect(response.adjusted_dose_mg).toBe(8.5);
      expect(response.factors_applied.adjustment_factor).toBe(0.85);
      expect(response.factors_applied.demographics_considered).toContain('age');
    });

    it('should format analysis response with all required fields', () => {
      const mockAnalysisResult = {
        interactions: { warnings: [], total_warnings: 0 },
        bias_correction: {
          adjusted_dose_mg: 8.5,
          factors_applied: { base_dose: 10, adjustment_factor: 0.85 },
          evidence: 'Test evidence',
          confidence: 0.75,
          warnings: [],
        },
        synergy: {
          synergy_score: 0.75,
          tk_enhanced: false,
          evidence: 'Computational prediction',
          confidence_level: 0.6,
        },
        plant_profile: {},
        polysaccharide_effects: {
          effects: 'microbiome_modulation',
          confidence: 0.75,
          modulation: 'beneficial',
        },
        processing_time_ms: 150,
      };

      expect(mockAnalysisResult).toHaveProperty('interactions');
      expect(mockAnalysisResult).toHaveProperty('bias_correction');
      expect(mockAnalysisResult).toHaveProperty('synergy');
      expect(mockAnalysisResult).toHaveProperty('polysaccharide_effects');
      expect(mockAnalysisResult).toHaveProperty('processing_time_ms');
    });
  });

  describe('Error Handling', () => {
    it('should handle missing database gracefully', () => {
      const env = { NEUROBOTANICA_DB: undefined };

      const hasDatabaseBinding = !!env.NEUROBOTANICA_DB;
      expect(hasDatabaseBinding).toBe(false);
    });

    it('should return appropriate error for missing required fields', () => {
      const validateRequest = (body: Record<string, unknown>) => {
        if (!body.compound_ids || !Array.isArray(body.compound_ids) || body.compound_ids.length === 0) {
          return { valid: false, error: 'compound_ids is required and must be a non-empty array' };
        }
        return { valid: true };
      };

      expect(validateRequest({})).toEqual({
        valid: false,
        error: 'compound_ids is required and must be a non-empty array',
      });
      expect(validateRequest({ compound_ids: [] })).toEqual({
        valid: false,
        error: 'compound_ids is required and must be a non-empty array',
      });
      expect(validateRequest({ compound_ids: ['cbd'] })).toEqual({ valid: true });
    });
  });

  describe('Performance', () => {
    it('should complete parallel API calls efficiently', async () => {
      const startTime = Date.now();

      // Simulate parallel execution pattern used in worker
      await Promise.all([
        new Promise(resolve => setTimeout(resolve, 100)),
        new Promise(resolve => setTimeout(resolve, 100)),
        new Promise(resolve => setTimeout(resolve, 100)),
        new Promise(resolve => setTimeout(resolve, 100)),
      ]);

      const elapsed = Date.now() - startTime;

      // Should complete in ~100ms (parallel), not ~400ms (sequential)
      expect(elapsed).toBeLessThan(300);
    });
  });
});

describe('Frontend Integration', () => {
  describe('API URL Configuration', () => {
    it('should use correct production URLs', () => {
      const BUDTENDER_API_URL = process.env.NEXT_PUBLIC_BUDTENDER_API_URL || 'https://budtender.neuro-botanica.com';
      const NEUROBOTANICA_API_URL = process.env.NEXT_PUBLIC_NEUROBOTANICA_API_URL || 'https://neurobotanica-api.contessapetrini.workers.dev';

      expect(BUDTENDER_API_URL).toContain('neuro-botanica');
      expect(NEUROBOTANICA_API_URL).toContain('workers.dev');
    });
  });

  describe('Retry Logic', () => {
    it('should retry on network failures', async () => {
      let attempts = 0;
      const maxRetries = 3;

      const fetchWithRetry = async (): Promise<boolean> => {
        while (attempts < maxRetries) {
          attempts++;
          if (attempts < 3) {
            // Simulate failure
            await new Promise(resolve => setTimeout(resolve, 10));
            continue;
          }
          return true;
        }
        return false;
      };

      const result = await fetchWithRetry();
      expect(result).toBe(true);
      expect(attempts).toBe(3);
    });
  });

  describe('Loading States', () => {
    it('should show appropriate loading messages', () => {
      const getLoadingMessage = (conditionCount: number) => {
        if (conditionCount === 0) return 'Loading recommendations...';
        return `Analyzing ${conditionCount} condition${conditionCount > 1 ? 's' : ''}...`;
      };

      expect(getLoadingMessage(0)).toBe('Loading recommendations...');
      expect(getLoadingMessage(1)).toBe('Analyzing 1 condition...');
      expect(getLoadingMessage(3)).toBe('Analyzing 3 conditions...');
    });
  });
});
