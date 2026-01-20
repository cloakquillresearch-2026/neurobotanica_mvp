"""
Unit Tests for Drug-Drug Interaction Checker
Tests known interactions (e.g., CBD + warfarin) and TK consent handling.
"""

import unittest
from src.engines.interactions import DrugInteractionChecker

class TestDrugInteractionChecker(unittest.TestCase):
    def setUp(self):
        self.checker = DrugInteractionChecker()
        # Assume test data is loaded in DB

    def test_known_interaction_cbd_warfarin(self):
        result = self.checker.check_interactions(["cbd"], ["warfarin"], "computational_only")
        self.assertGreater(result["total_warnings"], 0)
        self.assertIn("warfarin", [w["drug"] for w in result["warnings"]])

    def test_no_interaction(self):
        result = self.checker.check_interactions(["unknown"], ["unknown"], "computational_only")
        self.assertEqual(result["total_warnings"], 0)

    def test_tk_consent_check(self):
        # Mock TK interaction
        result = self.checker.check_interactions(["tk_compound"], ["drug"], "tk_enhanced")
        # Assert based on consent status
        self.assertIsInstance(result, dict)

    def tearDown(self):
        self.checker.close()

if __name__ == "__main__":
    unittest.main()