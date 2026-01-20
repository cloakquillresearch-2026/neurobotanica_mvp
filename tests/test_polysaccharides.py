"""
Unit Tests for Cross-Kingdom Polysaccharide Integration
Tests SCFA modeling, microbiome interactions, indigenous sourcing.
"""

import unittest
from src.engines.polysaccharides import PolysaccharideIntegrationEngine

class TestPolysaccharideIntegrationEngine(unittest.TestCase):
    def setUp(self):
        self.engine = PolysaccharideIntegrationEngine()

    def test_basic_polysaccharide_prediction(self):
        microbiome = {"bifidobacteria": 100}
        result = self.engine.predict_cross_kingdom_effects("beta_glucan_001", microbiome, "computational_only")
        self.assertIsInstance(result, dict)
        self.assertIn("scfa_production", result)
        self.assertIn("microbiome_effects", result)

    def test_tk_enhanced_prediction(self):
        microbiome = {"bifidobacteria": 100}
        result = self.engine.predict_cross_kingdom_effects("beta_glucan_001", microbiome, "tk_enhanced")
        self.assertIn("tk_enhanced", result)

    def test_scfa_modeling(self):
        poly_data = {"fermentation_rate": "fast"}
        microbiome = {}
        scfa = self.engine._model_scfa_production(poly_data, microbiome)
        self.assertIn("acetate", scfa)
        self.assertGreater(scfa["acetate"], 0)

    def tearDown(self):
        self.engine.close()

if __name__ == "__main__":
    unittest.main()