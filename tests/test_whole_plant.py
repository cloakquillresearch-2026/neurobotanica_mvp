"""
Unit Tests for Whole-Plant Analysis Engine
Tests profile aggregation, ratio computation, TK integration.
"""

import unittest
from src.engines.whole_plant import WholePlantAnalysisEngine

class TestWholePlantAnalysisEngine(unittest.TestCase):
    def setUp(self):
        self.analyzer = WholePlantAnalysisEngine()

    def test_basic_plant_analysis(self):
        result = self.analyzer.analyze_whole_plant("cbd", "anxiety", "computational_only")
        self.assertIsInstance(result, dict)
        self.assertIn("cannabinoid_ratios", result)
        self.assertFalse(result["tk_enhanced"])

    def test_tk_enhanced_analysis(self):
        result = self.analyzer.analyze_whole_plant("cbd", "anxiety", "tk_enhanced")
        self.assertIn("tk_enhanced", result)

    def test_ratio_computation(self):
        compounds = [{"compound_id": "cbd"}, {"compound_id": "thc"}]
        ratios = self.analyzer._compute_cannabinoid_ratios(compounds)
        self.assertEqual(len(ratios), 2)
        self.assertAlmostEqual(sum(ratios.values()), 1.0)

    def tearDown(self):
        self.analyzer.close()

if __name__ == "__main__":
    unittest.main()