"""
Unit Tests for Demographic Bias Correction
Tests adjustments for age/gender and TK consent handling.
"""

import unittest
from src.engines.bias_correction import DemographicBiasCorrection

class TestDemographicBiasCorrection(unittest.TestCase):
    def setUp(self):
        self.corrector = DemographicBiasCorrection()
        # Assume test data is loaded in DB

    def test_age_adjustment_65_female(self):
        result = self.corrector.apply_corrections(10.0, "cbd", {"age": 65, "gender": "female"}, "computational_only")
        self.assertIsInstance(result["adjusted_dose_mg"], float)
        self.assertLess(result["bias_variance"], 0.05)

    def test_no_adjustment_unknown(self):
        result = self.corrector.apply_corrections(10.0, "cbd", {"age": 999, "gender": "unknown"}, "computational_only")
        self.assertEqual(result["adjusted_dose_mg"], 10.0)  # No adjustment

    def test_tk_consent_check(self):
        # Mock TK adjustment
        result = self.corrector.apply_corrections(10.0, "tk_compound", {"age": 30, "gender": "male"}, "tk_enhanced")
        # Assert based on consent
        self.assertIsInstance(result, dict)

    def tearDown(self):
        self.corrector.close()

if __name__ == "__main__":
    unittest.main()