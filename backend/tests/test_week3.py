"""
Week 3 Component Tests
Tests for Clinical Evidence API, Receptor Affinity Provenance, and AssayAnalyzer

Reference: NeuroBotanica MVP Development Plan - Week 3
"""
import pytest
import math
from datetime import datetime
from unittest.mock import MagicMock, patch


# ==================== AssayAnalyzer Tests ====================

class TestAssayAnalyzer:
    """Test suite for AssayAnalyzer service."""
    
    @pytest.fixture
    def analyzer(self):
        """Create analyzer instance."""
        from backend.services.assay_analyzer import AssayAnalyzer
        return AssayAnalyzer(strict_mode=False)
    
    @pytest.fixture
    def analyzer_strict(self):
        """Create strict mode analyzer."""
        from backend.services.assay_analyzer import AssayAnalyzer
        return AssayAnalyzer(strict_mode=True)
    
    @pytest.fixture
    def homogeneous_measurements(self):
        """Homogeneous measurement set."""
        return [
            {"id": 1, "value": 10.0, "unit": "nM", "assay_type": "radioligand_binding", 
             "source_type": "peer_reviewed", "confidence": 0.8, "organism": "Homo sapiens"},
            {"id": 2, "value": 12.0, "unit": "nM", "assay_type": "radioligand_binding",
             "source_type": "peer_reviewed", "confidence": 0.85, "organism": "Homo sapiens"},
            {"id": 3, "value": 9.5, "unit": "nM", "assay_type": "radioligand_binding",
             "source_type": "database_curated", "confidence": 0.75, "organism": "Homo sapiens"},
        ]
    
    @pytest.fixture
    def heterogeneous_measurements(self):
        """Heterogeneous measurement set with multiple issues."""
        return [
            {"id": 1, "value": 10.0, "unit": "nM", "assay_type": "radioligand_binding",
             "source_type": "peer_reviewed", "confidence": 0.8, "organism": "Homo sapiens"},
            {"id": 2, "value": 0.5, "unit": "µM", "assay_type": "functional_cAMP",
             "source_type": "vendor_data", "confidence": 0.4, "organism": "Rattus norvegicus"},
            {"id": 3, "value": 5000.0, "unit": "nM", "assay_type": "beta_arrestin",
             "source_type": "preprint", "confidence": 0.5, "organism": "Homo sapiens"},
        ]
    
    def test_analyze_homogeneous_data(self, analyzer, homogeneous_measurements):
        """Test analysis of homogeneous data."""
        report = analyzer.analyze(
            homogeneous_measurements,
            compound_name="THC",
            receptor="CB1"
        )
        
        assert report.compound_name == "THC"
        assert report.receptor == "CB1"
        assert report.total_measurements == 3
        assert report.heterogeneity_score < 0.3  # Low score for homogeneous data
        assert report.can_aggregate is True
        assert len(report.issues) <= 1  # Minimal issues
    
    def test_analyze_heterogeneous_data(self, analyzer, heterogeneous_measurements):
        """Test analysis of heterogeneous data."""
        report = analyzer.analyze(
            heterogeneous_measurements,
            compound_name="CBD",
            receptor="CB2"
        )
        
        assert report.total_measurements == 3
        assert report.heterogeneity_score > 0.3  # Higher score for heterogeneous data
        assert len(report.issues) > 0  # Should detect issues
        
        # Check for expected issue types
        issue_types = [i.issue_type for i in report.issues]
        assert any("unit" in t or "assay" in t or "species" in t for t in issue_types)
    
    def test_detect_unit_mismatch(self, analyzer):
        """Test unit mismatch detection."""
        measurements = [
            {"id": 1, "value": 10.0, "unit": "nM"},
            {"id": 2, "value": 0.5, "unit": "µM"},
            {"id": 3, "value": 15.0, "unit": "nM"},
        ]
        
        report = analyzer.analyze(measurements, "Test", "CB1")
        
        # Should have unit-related issue
        unit_issues = [i for i in report.issues if "unit" in i.issue_type]
        assert len(unit_issues) > 0
    
    def test_detect_assay_incompatibility(self, analyzer):
        """Test assay type incompatibility detection."""
        measurements = [
            {"id": 1, "value": 10.0, "unit": "nM", "assay_type": "radioligand_binding"},
            {"id": 2, "value": 15.0, "unit": "nM", "assay_type": "functional_calcium"},
            {"id": 3, "value": 12.0, "unit": "nM", "assay_type": "BRET"},
        ]
        
        report = analyzer.analyze(measurements, "Test", "CB1")
        
        # Should detect assay incompatibility
        assay_issues = [i for i in report.issues if "assay" in i.issue_type]
        assert len(assay_issues) > 0
    
    def test_detect_high_variance(self, analyzer):
        """Test high variance detection."""
        measurements = [
            {"id": 1, "value": 1.0, "unit": "nM"},
            {"id": 2, "value": 100.0, "unit": "nM"},
            {"id": 3, "value": 10000.0, "unit": "nM"},
        ]
        
        report = analyzer.analyze(measurements, "Test", "CB1")
        
        # Should detect high variance
        variance_issues = [i for i in report.issues if "variance" in i.issue_type or "cv" in i.issue_type]
        assert len(variance_issues) > 0
    
    def test_detect_outliers(self, analyzer):
        """Test outlier detection."""
        measurements = [
            {"id": 1, "value": 10.0, "unit": "nM"},
            {"id": 2, "value": 12.0, "unit": "nM"},
            {"id": 3, "value": 11.0, "unit": "nM"},
            {"id": 4, "value": 9.0, "unit": "nM"},
            {"id": 5, "value": 1000.0, "unit": "nM"},  # Outlier
        ]
        
        report = analyzer.analyze(measurements, "Test", "CB1")
        
        outlier_issues = [i for i in report.issues if i.issue_type == "outliers"]
        assert len(outlier_issues) > 0
        assert 5 in outlier_issues[0].affected_measurements
    
    def test_detect_species_variance(self, analyzer):
        """Test species variance detection."""
        measurements = [
            {"id": 1, "value": 10.0, "unit": "nM", "organism": "Homo sapiens"},
            {"id": 2, "value": 15.0, "unit": "nM", "organism": "Rattus norvegicus"},
            {"id": 3, "value": 12.0, "unit": "nM", "organism": "Mus musculus"},
        ]
        
        report = analyzer.analyze(measurements, "Test", "CB1")
        
        species_issues = [i for i in report.issues if "species" in i.issue_type]
        assert len(species_issues) > 0
    
    def test_empty_measurements(self, analyzer):
        """Test handling of empty measurement list."""
        report = analyzer.analyze([], "Test", "CB1")
        
        assert report.total_measurements == 0
        assert report.can_aggregate is False
        assert len(report.aggregation_warnings) > 0
    
    def test_single_measurement(self, analyzer):
        """Test handling of single measurement."""
        measurements = [{"id": 1, "value": 10.0, "unit": "nM"}]
        
        report = analyzer.analyze(measurements, "Test", "CB1")
        
        assert report.total_measurements == 1
        assert report.can_aggregate is True  # Single measurement can be "aggregated"
    
    def test_strict_mode_lower_thresholds(self, analyzer, analyzer_strict):
        """Test that strict mode uses lower thresholds."""
        measurements = [
            {"id": 1, "value": 10.0, "unit": "nM"},
            {"id": 2, "value": 20.0, "unit": "nM"},
            {"id": 3, "value": 15.0, "unit": "nM"},
        ]
        
        normal_report = analyzer.analyze(measurements, "Test", "CB1")
        strict_report = analyzer_strict.analyze(measurements, "Test", "CB1")
        
        # Strict mode may detect more issues
        assert strict_report.heterogeneity_score >= normal_report.heterogeneity_score
    
    def test_aggregate_with_confidence_weighting(self, analyzer, homogeneous_measurements):
        """Test confidence-weighted aggregation."""
        result = analyzer.aggregate_with_weights(
            homogeneous_measurements,
            weighting="confidence"
        )
        
        assert "weighted_mean" in result
        assert "geometric_mean" in result
        assert "median" in result
        assert result["n_measurements"] == 3
        assert result["unit"] == "nM"
        assert result["weighting_method"] == "confidence"
    
    def test_aggregate_with_source_weighting(self, analyzer, homogeneous_measurements):
        """Test source-weighted aggregation."""
        result = analyzer.aggregate_with_weights(
            homogeneous_measurements,
            weighting="source"
        )
        
        assert "weighted_mean" in result
        assert result["weighting_method"] == "source"
    
    def test_aggregate_with_equal_weighting(self, analyzer, homogeneous_measurements):
        """Test equal-weighted aggregation."""
        result = analyzer.aggregate_with_weights(
            homogeneous_measurements,
            weighting="equal"
        )
        
        assert result["weighting_method"] == "equal"
        # Equal weighting should give arithmetic mean
        expected_mean = (10.0 + 12.0 + 9.5) / 3
        assert abs(result["weighted_mean"] - expected_mean) < 0.01
    
    def test_aggregate_normalizes_units(self, analyzer):
        """Test that aggregation normalizes units."""
        measurements = [
            {"id": 1, "value": 10.0, "unit": "nM", "confidence": 0.8},
            {"id": 2, "value": 0.015, "unit": "µM", "confidence": 0.8},  # 15 nM
        ]
        
        result = analyzer.aggregate_with_weights(measurements, "equal")
        
        # Should be around 12.5 nM after normalization
        assert 10 < result["weighted_mean"] < 15
        assert result["unit"] == "nM"


# ==================== Receptor Affinity Provenance Tests ====================

class TestReceptorAffinityModels:
    """Test receptor affinity models and helper functions."""
    
    def test_calculate_confidence_score(self):
        """Test confidence score calculation."""
        from backend.models.receptor_affinity import calculate_confidence_score
        
        # High quality peer-reviewed
        score_high = calculate_confidence_score(
            source_type="peer_reviewed",
            has_pubmed=True,
            has_chembl=True,
            data_quality="peer_reviewed"
        )
        
        # Low quality vendor data
        score_low = calculate_confidence_score(
            source_type="vendor_data",
            has_pubmed=False,
            has_chembl=False,
            data_quality="uncurated"
        )
        
        assert score_high > score_low
        assert 0 <= score_high <= 1
        assert 0 <= score_low <= 1
        assert score_high >= 0.7  # High quality should score well
    
    def test_normalize_affinity_unit(self):
        """Test unit normalization."""
        from backend.models.receptor_affinity import normalize_affinity_unit
        
        # pM to nM
        assert normalize_affinity_unit(1000, "pM", "nM") == 1.0
        
        # µM to nM
        assert normalize_affinity_unit(1, "µM", "nM") == 1000.0
        
        # nM to nM
        assert normalize_affinity_unit(10, "nM", "nM") == 10.0
        
        # nM to µM
        assert normalize_affinity_unit(1000, "nM", "µM") == 1.0
    
    def test_detect_assay_heterogeneity_helper(self):
        """Test heterogeneity detection helper."""
        from backend.models.receptor_affinity import detect_assay_heterogeneity
        
        # Heterogeneous data
        measurements = [
            {"affinity": {"value": 10, "unit": "nM"}, "assay": {"type": "binding"}},
            {"affinity": {"value": 0.5, "unit": "µM"}, "assay": {"type": "functional"}},
        ]
        
        result = detect_assay_heterogeneity(measurements)
        
        assert result["detected"] is True
        assert len(result["issues"]) > 0
    
    def test_detect_assay_heterogeneity_homogeneous(self):
        """Test homogeneous data returns no issues."""
        from backend.models.receptor_affinity import detect_assay_heterogeneity
        
        measurements = [
            {"affinity": {"value": 10, "unit": "nM"}, "assay": {"type": "binding"}},
            {"affinity": {"value": 12, "unit": "nM"}, "assay": {"type": "binding"}},
        ]
        
        result = detect_assay_heterogeneity(measurements)
        
        assert result["detected"] is False or len(result["issues"]) == 0


# ==================== Evidence API Tests ====================

class TestEvidenceHelpers:
    """Test evidence API helper functions."""
    
    def test_calculate_evidence_strength(self):
        """Test evidence strength calculation."""
        from backend.api.evidence import calculate_evidence_strength
        
        # Strong evidence
        strength = calculate_evidence_strength(
            total_studies=10,
            favorable_ratio=0.8,
            average_confidence=0.9
        )
        assert strength == "strong"
        
        # Insufficient evidence
        strength = calculate_evidence_strength(
            total_studies=0,
            favorable_ratio=0,
            average_confidence=0
        )
        assert strength == "insufficient"
        
        # Limited evidence
        strength = calculate_evidence_strength(
            total_studies=1,
            favorable_ratio=0.5,
            average_confidence=0.5
        )
        assert strength == "limited"
    
    def test_calculate_evidence_grade(self):
        """Test evidence grade calculation."""
        from backend.api.evidence import calculate_evidence_grade
        
        # Grade A
        grade = calculate_evidence_grade(
            weighted_score=0.85,
            num_studies=10,
            has_rct=True
        )
        assert grade == "A"
        
        # Grade B
        grade = calculate_evidence_grade(
            weighted_score=0.65,
            num_studies=5,
            has_rct=False
        )
        assert grade == "B"
        
        # Grade D
        grade = calculate_evidence_grade(
            weighted_score=0.2,
            num_studies=1,
            has_rct=False
        )
        assert grade == "D"
    
    def test_parse_outcome_direction(self):
        """Test outcome direction parsing."""
        from backend.api.evidence import parse_outcome_direction
        
        # Favorable
        assert parse_outcome_direction("Significant improvement observed") == "favorable"
        assert parse_outcome_direction("Treatment was effective") == "favorable"
        assert parse_outcome_direction("Reduced symptoms significantly") == "favorable"
        
        # Negative
        assert parse_outcome_direction("No significant effect observed") == "negative"
        assert parse_outcome_direction("Treatment was ineffective") == "negative"
        
        # Neutral
        assert parse_outcome_direction("Mixed results observed") == "neutral"
        assert parse_outcome_direction("") == "neutral"
        assert parse_outcome_direction(None) == "neutral"
    
    def test_generate_recommendation(self):
        """Test recommendation generation."""
        from backend.api.evidence import generate_recommendation
        
        # Grade A recommendation
        rec = generate_recommendation("A", "chronic pain", "CBD", 0.85)
        assert "Strong evidence" in rec
        assert "CBD" in rec
        
        # Grade D recommendation
        rec = generate_recommendation("D", "unknown condition", None, 0.2)
        assert "Insufficient evidence" in rec


# ==================== Integration Tests ====================

class TestWeek3Integration:
    """Integration tests for Week 3 components."""
    
    def test_analyzer_to_report_workflow(self):
        """Test full analyzer workflow produces valid report."""
        from backend.services.assay_analyzer import AssayAnalyzer, HeterogeneitReport
        
        analyzer = AssayAnalyzer()
        measurements = [
            {"id": 1, "value": 10.0, "unit": "nM", "assay_type": "binding",
             "source_type": "peer_reviewed", "confidence": 0.8, "organism": "Homo sapiens"},
            {"id": 2, "value": 15.0, "unit": "nM", "assay_type": "binding",
             "source_type": "database_curated", "confidence": 0.7, "organism": "Homo sapiens"},
        ]
        
        report = analyzer.analyze(measurements, "THC", "CB1")
        
        # Verify report structure
        assert isinstance(report, HeterogeneitReport)
        assert report.compound_name == "THC"
        assert report.receptor == "CB1"
        assert isinstance(report.heterogeneity_score, float)
        assert isinstance(report.can_aggregate, bool)
        assert isinstance(report.issues, list)
        assert isinstance(report.recommended_actions, list)
    
    def test_aggregation_follows_analysis(self):
        """Test that aggregation uses analysis results."""
        from backend.services.assay_analyzer import AssayAnalyzer
        
        analyzer = AssayAnalyzer()
        
        # Homogeneous data
        homogeneous = [
            {"id": 1, "value": 10.0, "unit": "nM", "confidence": 0.8},
            {"id": 2, "value": 12.0, "unit": "nM", "confidence": 0.7},
        ]
        
        report = analyzer.analyze(homogeneous, "Test", "CB1")
        assert report.can_aggregate is True
        
        result = analyzer.aggregate_with_weights(homogeneous, "confidence")
        assert "error" not in result
        assert result["n_measurements"] == 2
    
    def test_heterogeneity_score_range(self):
        """Test heterogeneity score is always 0-1."""
        from backend.services.assay_analyzer import AssayAnalyzer
        
        analyzer = AssayAnalyzer()
        
        test_cases = [
            [],  # Empty
            [{"id": 1, "value": 10.0, "unit": "nM"}],  # Single
            [{"id": i, "value": 10.0 * i, "unit": "nM"} for i in range(1, 20)],  # Many
        ]
        
        for measurements in test_cases:
            report = analyzer.analyze(measurements, "Test", "CB1")
            assert 0 <= report.heterogeneity_score <= 1


# ==================== Performance Tests ====================

class TestPerformance:
    """Performance tests for Week 3 components."""
    
    def test_analyzer_performance_large_dataset(self):
        """Test analyzer performance with large dataset."""
        import time
        from backend.services.assay_analyzer import AssayAnalyzer
        
        analyzer = AssayAnalyzer()
        
        # Create 1000 measurements
        measurements = [
            {
                "id": i,
                "value": 10.0 + (i % 50),
                "unit": "nM",
                "assay_type": "radioligand_binding",
                "source_type": "peer_reviewed",
                "confidence": 0.5 + (i % 5) * 0.1,
                "organism": "Homo sapiens"
            }
            for i in range(1000)
        ]
        
        start = time.time()
        report = analyzer.analyze(measurements, "Test", "CB1")
        elapsed = time.time() - start
        
        # Should complete in under 1 second
        assert elapsed < 1.0, f"Analysis took {elapsed:.2f}s"
        assert report.total_measurements == 1000
    
    def test_aggregation_performance(self):
        """Test aggregation performance with large dataset."""
        import time
        from backend.services.assay_analyzer import AssayAnalyzer
        
        analyzer = AssayAnalyzer()
        
        measurements = [
            {"id": i, "value": 10.0 + i, "unit": "nM", "confidence": 0.8}
            for i in range(1000)
        ]
        
        start = time.time()
        result = analyzer.aggregate_with_weights(measurements, "confidence")
        elapsed = time.time() - start
        
        assert elapsed < 0.5, f"Aggregation took {elapsed:.2f}s"
        assert result["n_measurements"] == 1000


# ==================== Run Tests ====================

if __name__ == "__main__":
    pytest.main([__file__, "-v"])
