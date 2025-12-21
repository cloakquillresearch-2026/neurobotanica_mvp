"""
Week 4 Tests: Dimeric Prediction + Triangulation Framework
Test suite for dimer prediction, conformer generation, and triangulation scoring

Reference: NeuroBotanica MVP Development Plan - Week 4
"""
import pytest
import numpy as np
import sys
from pathlib import Path
from unittest.mock import Mock, patch, MagicMock
from datetime import datetime

# Add backend to path
sys.path.insert(0, str(Path(__file__).parent.parent))

# Import services
from services.dimer_predictor import (
    DimericPredictor,
    DimerPrediction,
    LinkageType,
    DimerType,
    ReactiveSite
)
from services.dimer_conformer_generator import (
    DimericConformerGenerator,
    DimericConformerResult,
    DimerGeometryAnalysis
)
from services.triangulation_scorer import (
    TriangulationScorer,
    TriangulationResult,
    ValidationEvidence,
    ValidationSource,
    ExperimentalStatus
)


# =============================================================================
# Test Data
# =============================================================================

# Known cannabinoid SMILES
THC_SMILES = "CCCCCc1cc(O)c2C3C=C(C)CCC3C(C)(C)Oc2c1"
CBD_SMILES = "CCCCCc1cc(O)c(C2C=C(C)CCC2C(C)=C)c(O)c1"
CBG_SMILES = "CCCCCc1cc(O)c(CC=C(C)CCC=C(C)C)c(O)c1"
CBN_SMILES = "CCCCCc1cc2c(c(c1)O)C1C=C(C)CCC1C(C)(C)O2"

# Simple test molecule
SIMPLE_PHENOL = "c1ccc(O)cc1"


# =============================================================================
# DimericPredictor Tests
# =============================================================================

class TestDimericPredictor:
    """Tests for DimericPredictor class."""
    
    def test_predictor_initialization(self):
        """Test predictor initializes correctly."""
        predictor = DimericPredictor()
        
        assert predictor is not None
        # Check for attributes that exist
        assert hasattr(predictor, 'known_dimers')
        assert hasattr(predictor, 'predict_homodimer')
    
    def test_identify_reactive_sites_phenol(self):
        """Test identification of phenolic reactive sites."""
        predictor = DimericPredictor()
        
        sites = predictor.identify_reactive_sites(CBG_SMILES)
        
        # CBG has phenolic OH groups
        assert len(sites) > 0
    
    def test_predict_homodimer_creates_valid_structure(self):
        """Test homodimer prediction creates valid result."""
        predictor = DimericPredictor()
        
        prediction = predictor.predict_homodimer(
            monomer_name="CBD",
            monomer_smiles=CBD_SMILES,
            linkage_type=LinkageType.METHYLENE
        )
        
        assert prediction is not None
        assert prediction.dimer_type == DimerType.HOMODIMER
        assert prediction.parent_1_name == "CBD"
        assert prediction.parent_2_name == "CBD"
        assert 0 <= prediction.formation_probability <= 1
    
    def test_predict_heterodimer_creates_valid_structure(self):
        """Test heterodimer prediction creates valid result."""
        predictor = DimericPredictor()
        
        prediction = predictor.predict_heterodimer(
            monomer1_name="THC",
            monomer1_smiles=THC_SMILES,
            monomer2_name="CBD",
            monomer2_smiles=CBD_SMILES,
            linkage_type=LinkageType.METHYLENE
        )
        
        assert prediction is not None
        assert prediction.dimer_type == DimerType.HETERODIMER
        assert prediction.parent_1_name == "THC"
        assert prediction.parent_2_name == "CBD"
    
    def test_synergy_score_for_known_pair(self):
        """Test synergy score is calculated."""
        predictor = DimericPredictor()
        
        prediction = predictor.predict_heterodimer(
            monomer1_name="THC",
            monomer1_smiles=THC_SMILES,
            monomer2_name="CBD",
            monomer2_smiles=CBD_SMILES
        )
        
        assert prediction is not None
        # Synergy should be calculated
        assert prediction.synergy_prediction >= 0
    
    def test_therapeutic_potential_predictions(self):
        """Test therapeutic potential predictions are generated."""
        predictor = DimericPredictor()
        
        prediction = predictor.predict_homodimer(
            monomer_name="CBD",
            monomer_smiles=CBD_SMILES
        )
        
        assert prediction is not None
        assert isinstance(prediction.therapeutic_potential, dict)
    
    def test_generate_all_combinations(self):
        """Test generation of all monomer combinations."""
        predictor = DimericPredictor()
        
        # Use list of tuples format
        monomers = [
            ("THC", THC_SMILES),
            ("CBD", CBD_SMILES),
        ]
        
        predictions = predictor.generate_all_combinations(
            monomers,
            include_homodimers=True
        )
        
        # Should generate: THC-THC, CBD-CBD, THC-CBD
        assert isinstance(predictions, list)
        assert len(predictions) > 0
    
    def test_rank_predictions(self):
        """Test prediction ranking by scores."""
        predictor = DimericPredictor()
        
        monomers = [
            ("THC", THC_SMILES),
            ("CBD", CBD_SMILES),
            ("CBG", CBG_SMILES),
        ]
        
        predictions = predictor.generate_all_combinations(monomers)
        
        if len(predictions) > 1:
            ranked = predictor.rank_predictions(predictions)
            assert len(ranked) == len(predictions)
    
    def test_linkage_types_enum(self):
        """Test all linkage types are available."""
        assert LinkageType.METHYLENE.value == "methylene"
        assert LinkageType.ETHER.value == "ether"
        assert LinkageType.ESTER.value == "ester"
        assert LinkageType.DIRECT.value == "direct"
        assert LinkageType.CARBON_CHAIN.value == "carbon_chain"
    
    def test_dimer_types_enum(self):
        """Test dimer type enumeration."""
        assert DimerType.HOMODIMER.value == "homodimer"
        assert DimerType.HETERODIMER.value == "heterodimer"


# =============================================================================
# TriangulationScorer Tests
# =============================================================================

class TestTriangulationScorer:
    """Tests for TriangulationScorer class."""
    
    def test_scorer_initialization(self):
        """Test scorer initializes with default weights."""
        scorer = TriangulationScorer()
        
        assert scorer is not None
        assert ValidationSource.COMPUTATIONAL in scorer.weights
        assert ValidationSource.STRUCTURAL in scorer.weights
        assert scorer.confidence_level == 0.95
    
    def test_custom_weights(self):
        """Test scorer accepts custom weights."""
        custom = {
            ValidationSource.COMPUTATIONAL: 0.5,
            ValidationSource.STRUCTURAL: 0.3,
            ValidationSource.ML_THERAPEUTIC: 0.2,
        }
        
        scorer = TriangulationScorer(custom_weights=custom)
        
        assert scorer.weights[ValidationSource.COMPUTATIONAL] == 0.5
    
    def test_triangulation_with_three_sources(self):
        """Test triangulation with standard three sources."""
        scorer = TriangulationScorer()
        
        result = scorer.calculate_triangulation_score(
            formation_probability=0.75,
            structural_similarity=0.65,
            ml_therapeutic_score=0.80
        )
        
        assert isinstance(result, TriangulationResult)
        assert 0 <= result.triangulation_score <= 1
        assert result.num_sources == 3
        assert "computational" in result.sources_used
        assert "structural" in result.sources_used
        assert "ml_therapeutic" in result.sources_used
    
    def test_triangulation_insufficient_sources(self):
        """Test triangulation with insufficient sources."""
        scorer = TriangulationScorer(require_minimum_sources=2)
        
        result = scorer.calculate_triangulation_score(
            formation_probability=0.75
        )
        
        assert result.validation_grade == "D"
        assert result.triangulation_score == 0.0
        assert "insufficient" in result.consensus_strength.lower()
    
    def test_confidence_interval_calculation(self):
        """Test confidence interval is properly bounded."""
        scorer = TriangulationScorer(confidence_level=0.95)
        
        result = scorer.calculate_triangulation_score(
            formation_probability=0.70,
            structural_similarity=0.65,
            ml_therapeutic_score=0.75
        )
        
        assert result.confidence_interval_lower <= result.triangulation_score
        assert result.confidence_interval_upper >= result.triangulation_score
        assert result.confidence_interval_lower >= 0
        assert result.confidence_interval_upper <= 1
    
    def test_experimental_confirmation_boosts_score(self):
        """Test experimental confirmation increases score."""
        scorer = TriangulationScorer()
        
        # Without experimental
        result_no_exp = scorer.calculate_triangulation_score(
            formation_probability=0.60,
            structural_similarity=0.55,
            ml_therapeutic_score=0.65
        )
        
        # With experimental confirmation
        result_exp = scorer.calculate_triangulation_score(
            formation_probability=0.60,
            structural_similarity=0.55,
            ml_therapeutic_score=0.65,
            experimental_validation="confirmed"
        )
        
        assert result_exp.triangulation_score > result_no_exp.triangulation_score
        assert result_exp.experimental_status == ExperimentalStatus.CONFIRMED
    
    def test_experimental_contradiction_reduces_score(self):
        """Test experimental contradiction reduces score."""
        scorer = TriangulationScorer()
        
        result_no_exp = scorer.calculate_triangulation_score(
            formation_probability=0.70,
            structural_similarity=0.65,
            ml_therapeutic_score=0.75
        )
        
        result_contra = scorer.calculate_triangulation_score(
            formation_probability=0.70,
            structural_similarity=0.65,
            ml_therapeutic_score=0.75,
            experimental_validation="contradicted"
        )
        
        assert result_contra.triangulation_score < result_no_exp.triangulation_score
        assert result_contra.experimental_status == ExperimentalStatus.CONTRADICTED
    
    def test_grade_assignment(self):
        """Test validation grade assignment."""
        scorer = TriangulationScorer()
        
        # High scores should get A
        result_high = scorer.calculate_triangulation_score(
            formation_probability=0.90,
            structural_similarity=0.85,
            ml_therapeutic_score=0.88
        )
        assert result_high.validation_grade in ["A", "B"]
        
        # Low scores should get C or D
        result_low = scorer.calculate_triangulation_score(
            formation_probability=0.30,
            structural_similarity=0.35,
            ml_therapeutic_score=0.40
        )
        assert result_low.validation_grade in ["C", "D"]
    
    def test_source_agreement_calculation(self):
        """Test source agreement metric."""
        scorer = TriangulationScorer()
        
        # High agreement - similar scores
        result_agree = scorer.calculate_triangulation_score(
            formation_probability=0.70,
            structural_similarity=0.72,
            ml_therapeutic_score=0.68
        )
        
        # Low agreement - varied scores
        result_disagree = scorer.calculate_triangulation_score(
            formation_probability=0.90,
            structural_similarity=0.30,
            ml_therapeutic_score=0.60
        )
        
        assert result_agree.source_agreement > result_disagree.source_agreement
    
    def test_consensus_strength_labels(self):
        """Test consensus strength is properly labeled."""
        scorer = TriangulationScorer()
        
        result = scorer.calculate_triangulation_score(
            formation_probability=0.75,
            structural_similarity=0.73,
            ml_therapeutic_score=0.77
        )
        
        assert result.consensus_strength in ["strong", "moderate", "weak", "conflicting"]
    
    def test_recommendations_generated(self):
        """Test recommendations are generated."""
        scorer = TriangulationScorer()
        
        result = scorer.calculate_triangulation_score(
            formation_probability=0.70,
            structural_similarity=0.65,
            ml_therapeutic_score=0.72
        )
        
        assert len(result.recommended_actions) > 0
        assert all(isinstance(r, str) for r in result.recommended_actions)
    
    def test_batch_triangulation(self):
        """Test batch processing of predictions."""
        scorer = TriangulationScorer()
        
        predictions = [
            {"formation_probability": 0.70, "structural_similarity": 0.65, "ml_therapeutic_score": 0.75},
            {"formation_probability": 0.80, "structural_similarity": 0.75, "ml_therapeutic_score": 0.85},
            {"formation_probability": 0.50, "structural_similarity": 0.55, "ml_therapeutic_score": 0.60},
        ]
        
        results = scorer.batch_triangulate(predictions)
        
        assert len(results) == 3
        assert all(isinstance(r, TriangulationResult) for r in results)
    
    def test_rank_by_confidence(self):
        """Test ranking by confidence."""
        scorer = TriangulationScorer()
        
        result1 = scorer.calculate_triangulation_score(
            formation_probability=0.50, structural_similarity=0.55, ml_therapeutic_score=0.60
        )
        result2 = scorer.calculate_triangulation_score(
            formation_probability=0.85, structural_similarity=0.80, ml_therapeutic_score=0.82
        )
        result3 = scorer.calculate_triangulation_score(
            formation_probability=0.65, structural_similarity=0.70, ml_therapeutic_score=0.68
        )
        
        items = [("pred1", result1), ("pred2", result2), ("pred3", result3)]
        ranked = scorer.rank_by_confidence(items)
        
        # Highest score should be first
        assert ranked[0][1].triangulation_score >= ranked[1][1].triangulation_score
        assert ranked[1][1].triangulation_score >= ranked[2][1].triangulation_score
    
    def test_to_dict_serialization(self):
        """Test TriangulationResult serializes to dict."""
        scorer = TriangulationScorer()
        
        result = scorer.calculate_triangulation_score(
            formation_probability=0.70,
            structural_similarity=0.65,
            ml_therapeutic_score=0.75
        )
        
        result_dict = result.to_dict()
        
        assert "triangulation_score" in result_dict
        assert "confidence_interval" in result_dict
        assert "validation_grade" in result_dict
        assert result_dict["triangulation_score"] == result.triangulation_score


# =============================================================================
# DimericConformerGenerator Tests
# =============================================================================

class TestDimericConformerGenerator:
    """Tests for DimericConformerGenerator class."""
    
    def test_generator_initialization(self):
        """Test generator initializes correctly."""
        generator = DimericConformerGenerator()
        
        assert generator is not None
        assert hasattr(generator, 'generate_dimer_conformers')
    
    def test_simple_molecule_conformers(self):
        """Test conformer generation for simple molecule."""
        generator = DimericConformerGenerator(num_conformers=5)
        
        # Use base class method
        result = generator.generate_conformers("CCCCCC")  # Hexane
        
        assert result is not None
        assert result.success or result.num_conformers_generated >= 1
    
    def test_dimer_conformers_generation(self):
        """Test dimer conformer generation."""
        generator = DimericConformerGenerator(num_conformers=5)
        
        # Simple dimer-like structure (biphenyl)
        result = generator.generate_dimer_conformers(
            dimer_smiles="c1ccccc1-c1ccccc1",
            monomer1_smiles="c1ccccc1",
            monomer2_smiles="c1ccccc1"
        )
        
        assert result is not None
        assert isinstance(result, DimericConformerResult)
    
    def test_conformer_generation_base(self):
        """Test base conformer generation works."""
        generator = DimericConformerGenerator(num_conformers=10)
        
        result = generator.generate_conformers("CCCCCCCCC")  # Nonane
        
        assert result is not None
        if result.success:
            assert result.num_conformers_generated >= 1


# =============================================================================
# Integration Tests
# =============================================================================

class TestDimerPipelineIntegration:
    """Integration tests for the full dimer prediction pipeline."""
    
    def test_full_prediction_to_triangulation(self):
        """Test full pipeline from prediction to triangulation."""
        predictor = DimericPredictor()
        scorer = TriangulationScorer()
        
        # Generate prediction
        prediction = predictor.predict_homodimer(
            monomer_name="CBD",
            monomer_smiles=CBD_SMILES
        )
        
        assert prediction is not None
        
        # Get therapeutic max
        therapeutic_max = max(prediction.therapeutic_potential.values()) if prediction.therapeutic_potential else 0.5
        
        # Calculate triangulation
        result = scorer.calculate_triangulation_score(
            formation_probability=prediction.formation_probability,
            structural_similarity=1 - prediction.novelty_score,
            ml_therapeutic_score=therapeutic_max
        )
        
        assert result.num_sources >= 2
        assert result.triangulation_score >= 0
    
    def test_batch_processing_pipeline(self):
        """Test batch processing of multiple dimers."""
        predictor = DimericPredictor()
        scorer = TriangulationScorer()
        
        monomers = [
            ("THC", THC_SMILES),
            ("CBD", CBD_SMILES),
        ]
        
        predictions = predictor.generate_all_combinations(monomers)
        
        # Convert to triangulation input format
        tri_inputs = []
        for pred in predictions:
            therapeutic_max = max(pred.therapeutic_potential.values()) if pred.therapeutic_potential else 0.5
            tri_inputs.append({
                "formation_probability": pred.formation_probability,
                "structural_similarity": 1 - pred.novelty_score,
                "ml_therapeutic_score": therapeutic_max
            })
        
        if tri_inputs:
            results = scorer.batch_triangulate(tri_inputs)
            assert len(results) == len(predictions)


# =============================================================================
# Validation Source Tests
# =============================================================================

class TestValidationSource:
    """Tests for ValidationSource and ValidationEvidence."""
    
    def test_validation_source_enum(self):
        """Test validation source enumeration."""
        assert ValidationSource.COMPUTATIONAL.value == "computational"
        assert ValidationSource.STRUCTURAL.value == "structural"
        assert ValidationSource.ML_THERAPEUTIC.value == "ml_therapeutic"
        assert ValidationSource.EXPERIMENTAL.value == "experimental"
    
    def test_validation_evidence_weighted_contribution(self):
        """Test weighted contribution calculation."""
        evidence = ValidationEvidence(
            source=ValidationSource.COMPUTATIONAL,
            score=0.80,
            confidence=0.90,
            weight=0.40
        )
        
        contribution = evidence.weighted_contribution()
        expected = 0.80 * 0.90 * 0.40
        
        assert abs(contribution - expected) < 0.001
    
    def test_experimental_status_enum(self):
        """Test experimental status enumeration."""
        assert ExperimentalStatus.CONFIRMED.value == "confirmed"
        assert ExperimentalStatus.TENTATIVE.value == "tentative"
        assert ExperimentalStatus.PREDICTED.value == "predicted"
        assert ExperimentalStatus.CONTRADICTED.value == "contradicted"


# =============================================================================
# Edge Cases
# =============================================================================

class TestEdgeCases:
    """Tests for edge cases and error handling."""
    
    def test_empty_monomers_list(self):
        """Test handling of empty monomers list."""
        predictor = DimericPredictor()
        
        predictions = predictor.generate_all_combinations([])
        
        assert predictions == []
    
    def test_invalid_smiles_prediction(self):
        """Test handling of invalid SMILES."""
        predictor = DimericPredictor()
        
        prediction = predictor.predict_homodimer(
            monomer_name="Invalid",
            monomer_smiles="invalid_smiles"
        )
        
        # Should return a failed prediction with low scores
        assert prediction is not None
        assert prediction.formation_probability == 0.0 or "failed" in prediction.confidence_level.lower()
    
    def test_triangulation_all_zeros(self):
        """Test triangulation with all zero scores."""
        scorer = TriangulationScorer()
        
        result = scorer.calculate_triangulation_score(
            formation_probability=0.0,
            structural_similarity=0.0,
            ml_therapeutic_score=0.0
        )
        
        assert result.triangulation_score == 0.0
        assert result.validation_grade == "D"
    
    def test_triangulation_all_ones(self):
        """Test triangulation with all perfect scores."""
        scorer = TriangulationScorer()
        
        result = scorer.calculate_triangulation_score(
            formation_probability=1.0,
            structural_similarity=1.0,
            ml_therapeutic_score=1.0
        )
        
        assert result.triangulation_score > 0.9
        assert result.validation_grade == "A"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
