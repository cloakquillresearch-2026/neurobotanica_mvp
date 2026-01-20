"""
Unit Tests for Synergy Prediction System
Tests ML predictions, TK enhancements, and consent handling.
"""

import unittest
from src.engines.synergy import SynergyPredictionSystem

class TestSynergyPredictionSystem(unittest.TestCase):
    def setUp(self):
        self.predictor = SynergyPredictionSystem()
        # Assume test data is loaded in DB

    def test_basic_synergy_prediction(self):
        result = self.predictor.predict_synergy("cbd", "thc", "computational_only")
        self.assertIsInstance(result["synergy_score"], float)
        self.assertGreaterEqual(result["synergy_score"], 0.0)
        self.assertLessEqual(result["synergy_score"], 1.0)
        self.assertIn("mechanisms", result)

    def test_tk_enhanced_prediction(self):
        result = self.predictor.predict_synergy("cbd", "thc", "tk_enhanced")
        self.assertIsInstance(result, dict)
        # Assert TK enhancement if consent verified
        self.assertIn("tk_enhanced", result)

    def test_confidence_calculation(self):
        confidence = self.predictor._calculate_confidence("cbd", "thc")
        self.assertIsInstance(confidence, float)
        self.assertGreaterEqual(confidence, 0.0)
        self.assertLessEqual(confidence, 1.0)

    def tearDown(self):
        self.predictor.close()

if __name__ == "__main__":
    unittest.main()