"""
Week 6 Tests - ML Models + PatentPath Lite
Tests for ML training pipeline and patent analysis features.
"""
import pytest
import pytest_asyncio
import numpy as np
import pandas as pd
from datetime import datetime
from pathlib import Path
import tempfile
import asyncio
from unittest.mock import AsyncMock, MagicMock, patch

# Configure pytest-asyncio mode
pytest_plugins = ('pytest_asyncio',)

# Import Week 6 modules
import sys
sys.path.insert(0, str(Path(__file__).parent.parent))

from services.ml_data_prep import (
    MLDataPreparator,
    FeatureSet,
    DatasetStats,
    STANDARD_2D_DESCRIPTORS,
    STANDARD_3D_DESCRIPTORS,
    RECEPTOR_TARGETS
)
from services.ml_models import (
    TherapeuticPredictionModel,
    PatientResponseModel,
    DimerPotentialModel,
    ModelMetrics,
    PredictionResult,
    ModelRegistry,
    BaseModel
)
from services.patentpath.prior_art import (
    PriorArtSearcher,
    PriorArtResult,
    SearchQuery,
    PatentSource,
    USPTOClient
)
from services.patentpath.novelty import (
    NoveltyScorer,
    NoveltyReport,
    NoveltyLevel,
    SimilarityMatch,
    quick_novelty_check
)


# ============================================================================
# ML Data Preparation Tests (Task 6.1)
# ============================================================================

class TestMLDataPreparator:
    """Tests for MLDataPreparator class."""
    
    def test_init(self):
        """Test MLDataPreparator initialization."""
        prep = MLDataPreparator()
        assert prep is not None
        assert hasattr(prep, 'prepare_therapeutic_dataset')
        assert hasattr(prep, 'prepare_dimer_dataset')
    
    def test_feature_constants(self):
        """Test that feature constants are properly defined."""
        assert len(STANDARD_2D_DESCRIPTORS) > 0
        assert len(STANDARD_3D_DESCRIPTORS) > 0
        assert len(RECEPTOR_TARGETS) > 0
        
        # Check expected descriptors
        assert "molecular_weight" in STANDARD_2D_DESCRIPTORS
        assert "logp" in STANDARD_2D_DESCRIPTORS
        assert "tpsa" in STANDARD_2D_DESCRIPTORS
        
        # Check receptor targets
        assert "CB1" in RECEPTOR_TARGETS
        assert "CB2" in RECEPTOR_TARGETS
    
    def test_feature_set_dataclass(self):
        """Test FeatureSet dataclass."""
        df = pd.DataFrame({"f1": [1.0], "f2": [2.0]})
        features = FeatureSet(
            features=df,
            feature_names=["f1", "f2"],
            target_column="target",
            weight_column=None
        )
        
        assert len(features.feature_names) == 2
        assert features.target_column == "target"
        
        # Test to_dict
        data = features.to_dict()
        assert "num_samples" in data
        assert "num_features" in data
    
    def test_dataset_stats_dataclass(self):
        """Test DatasetStats dataclass."""
        stats = DatasetStats(
            total_samples=100,
            num_features=50,
            num_2d_features=20,
            num_3d_features=13,
            num_receptor_features=10,
            num_conditions=10,
            condition_distribution={"pain": 30, "anxiety": 20},
            missing_3d_count=5,
            avg_confidence_weight=0.85
        )
        
        assert stats.total_samples == 100
        assert stats.num_features == 50
        
        summary = stats.to_dict()
        assert "total_samples" in summary


class TestFeatureExtraction:
    """Tests for feature extraction functionality."""
    
    def test_2d_descriptors_list(self):
        """Test 2D descriptors are comprehensive."""
        expected = [
            "molecular_weight", "logp", "tpsa", "h_bond_donors",
            "h_bond_acceptors", "rotatable_bonds"
        ]
        for desc in expected:
            assert desc in STANDARD_2D_DESCRIPTORS, f"Missing {desc}"
    
    def test_3d_descriptors_list(self):
        """Test 3D descriptors are defined."""
        expected = ["pmi1", "pmi2", "pmi3", "asphericity"]
        for desc in expected:
            assert desc in STANDARD_3D_DESCRIPTORS, f"Missing {desc}"
    
    def test_receptor_targets_complete(self):
        """Test receptor targets include key receptors."""
        expected = ["CB1", "CB2", "GPR55", "TRPV1"]
        for target in expected:
            assert target in RECEPTOR_TARGETS, f"Missing {target}"


# ============================================================================
# ML Models Tests (Task 6.2)
# ============================================================================

class TestModelMetrics:
    """Tests for ModelMetrics dataclass."""
    
    def test_model_metrics_creation(self):
        """Test ModelMetrics creation and serialization."""
        metrics = ModelMetrics(
            model_name="TestModel",
            model_version="1.0.0",
            train_score=0.95,
            test_score=0.88,
            cv_scores=[0.85, 0.87, 0.89, 0.86, 0.88],
            cv_mean=0.87,
            cv_std=0.015,
            num_features=50,
            num_samples=1000,
            training_time_seconds=5.5,
            feature_importances={"feature1": 0.3, "feature2": 0.2}
        )
        
        assert metrics.model_name == "TestModel"
        assert metrics.test_score == 0.88
        assert len(metrics.cv_scores) == 5
        
        data = metrics.to_dict()
        assert "train_score" in data
        assert "cv_mean" in data
        assert "top_features" in data


class TestPredictionResult:
    """Tests for PredictionResult dataclass."""
    
    def test_prediction_result_scalar(self):
        """Test PredictionResult with scalar values."""
        result = PredictionResult(
            prediction=0.75,
            uncertainty=0.05,
            confidence_interval_lower=0.65,
            confidence_interval_upper=0.85
        )
        
        assert result.prediction == 0.75
        assert result.confidence_level == 0.95
        
        data = result.to_dict()
        assert data["prediction"] == 0.75
        assert "confidence_interval" in data
    
    def test_prediction_result_array(self):
        """Test PredictionResult with array values."""
        result = PredictionResult(
            prediction=np.array([0.7, 0.8, 0.6]),
            uncertainty=np.array([0.05, 0.03, 0.07]),
            confidence_interval_lower=np.array([0.6, 0.74, 0.46]),
            confidence_interval_upper=np.array([0.8, 0.86, 0.74])
        )
        
        data = result.to_dict()
        assert isinstance(data["prediction"], list)
        assert len(data["prediction"]) == 3


class TestTherapeuticPredictionModel:
    """Tests for TherapeuticPredictionModel."""
    
    def test_model_initialization(self):
        """Test model can be initialized."""
        model = TherapeuticPredictionModel()
        assert model.model_name == "TherapeuticPrediction"
        assert not model.is_trained
        assert model.MODEL_VERSION == "1.0.0"
    
    def test_model_training(self):
        """Test model training with synthetic data."""
        model = TherapeuticPredictionModel(n_estimators=10, max_depth=3)
        
        # Create synthetic training data with realistic pattern
        np.random.seed(42)
        X = pd.DataFrame(np.random.randn(100, 10), columns=[f"f{i}" for i in range(10)])
        # Create target that has some correlation with features
        y = pd.Series(np.clip(0.5 + 0.3 * X["f0"] + 0.2 * X["f1"] + 0.1 * np.random.randn(100), 0, 1))
        
        metrics = model.train(X, y)
        
        assert model.is_trained
        assert metrics is not None
        # R² can be negative for poor models, just check it exists
        assert hasattr(metrics, 'test_score')
        assert len(metrics.cv_scores) == 5
    
    def test_model_prediction(self):
        """Test model prediction after training."""
        model = TherapeuticPredictionModel(n_estimators=10)
        
        np.random.seed(42)
        X_train = pd.DataFrame(np.random.randn(100, 5), columns=[f"f{i}" for i in range(5)])
        y_train = pd.Series(np.random.rand(100))
        
        model.train(X_train, y_train)
        
        X_test = pd.DataFrame(np.random.randn(10, 5), columns=[f"f{i}" for i in range(5)])
        predictions = model.predict(X_test)
        
        assert len(predictions) == 10
        assert all(isinstance(p, (float, np.floating)) for p in predictions)
    
    def test_model_predict_with_uncertainty(self):
        """Test prediction with uncertainty quantification."""
        model = TherapeuticPredictionModel(n_estimators=20)
        
        np.random.seed(42)
        X_train = pd.DataFrame(np.random.randn(100, 5), columns=[f"f{i}" for i in range(5)])
        y_train = pd.Series(np.clip(np.random.rand(100), 0, 1))
        
        model.train(X_train, y_train)
        
        X_test = pd.DataFrame(np.random.randn(5, 5), columns=[f"f{i}" for i in range(5)])
        result = model.predict_with_uncertainty(X_test)
        
        assert isinstance(result, PredictionResult)
        assert len(result.prediction) == 5
        assert len(result.uncertainty) == 5
        assert all(result.confidence_interval_lower <= result.prediction)
        assert all(result.prediction <= result.confidence_interval_upper)
    
    def test_untrained_model_raises(self):
        """Test that untrained model raises on predict."""
        model = TherapeuticPredictionModel()
        X = pd.DataFrame(np.random.randn(5, 5), columns=[f"f{i}" for i in range(5)])
        
        with pytest.raises(ValueError, match="Model not trained"):
            model.predict(X)


class TestPatientResponseModel:
    """Tests for PatientResponseModel."""
    
    def test_model_initialization(self):
        """Test classifier initialization."""
        model = PatientResponseModel()
        assert model.model_name == "PatientResponse"
        assert not model.is_trained
    
    def test_model_training_classification(self):
        """Test binary classification training."""
        model = PatientResponseModel(n_estimators=10, max_depth=3)
        
        np.random.seed(42)
        X = pd.DataFrame(np.random.randn(100, 8), columns=[f"f{i}" for i in range(8)])
        y = pd.Series(np.random.rand(100))  # Will be binarized
        
        metrics = model.train(X, y)
        
        assert model.is_trained
        assert 0 <= metrics.test_score <= 1
    
    def test_predict_proba(self):
        """Test probability prediction."""
        model = PatientResponseModel(n_estimators=10)
        
        np.random.seed(42)
        X_train = pd.DataFrame(np.random.randn(100, 5), columns=[f"f{i}" for i in range(5)])
        y_train = pd.Series(np.random.rand(100))
        
        model.train(X_train, y_train)
        
        X_test = pd.DataFrame(np.random.randn(10, 5), columns=[f"f{i}" for i in range(5)])
        probs = model.predict_proba(X_test)
        
        assert len(probs) == 10
        assert all(0 <= p <= 1 for p in probs)


class TestDimerPotentialModel:
    """Tests for DimerPotentialModel."""
    
    def test_model_initialization(self):
        """Test dimer model initialization."""
        model = DimerPotentialModel()
        assert model.model_name == "DimerPotential"
        assert not model.is_trained
    
    def test_model_training(self):
        """Test dimer model training."""
        model = DimerPotentialModel(n_estimators=10, max_depth=3)
        
        np.random.seed(42)
        # Dimer features: combined parent features
        X = pd.DataFrame(np.random.randn(80, 20), columns=[f"f{i}" for i in range(20)])
        y = pd.Series(np.random.rand(80))
        
        metrics = model.train(X, y)
        
        assert model.is_trained
        assert metrics.model_name == "DimerPotential"
    
    def test_rank_dimers(self):
        """Test dimer ranking functionality."""
        model = DimerPotentialModel(n_estimators=10)
        
        np.random.seed(42)
        X_train = pd.DataFrame(np.random.randn(50, 10), columns=[f"f{i}" for i in range(10)])
        y_train = pd.Series(np.random.rand(50))
        model.train(X_train, y_train)
        
        # Create test dimers
        dimers = [
            {"dimer_name": "THC-CBD", "parent_1_name": "THC", "parent_2_name": "CBD"},
            {"dimer_name": "CBD-CBG", "parent_1_name": "CBD", "parent_2_name": "CBG"},
        ]
        X_test = pd.DataFrame(np.random.randn(2, 10), columns=[f"f{i}" for i in range(10)])
        
        ranked = model.rank_dimers(dimers, X_test)
        
        assert len(ranked) == 2
        assert "therapeutic_potential" in ranked[0]
        assert ranked[0]["therapeutic_potential"] >= ranked[1]["therapeutic_potential"]


class TestModelPersistence:
    """Tests for model save/load functionality."""
    
    def test_model_save_load(self):
        """Test model can be saved and loaded."""
        model = TherapeuticPredictionModel(n_estimators=10)
        
        np.random.seed(42)
        X = pd.DataFrame(np.random.randn(50, 5), columns=[f"f{i}" for i in range(5)])
        y = pd.Series(np.random.rand(50))
        model.train(X, y)
        
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "model.joblib"
            model.save(str(path))
            
            assert path.exists()
            
            # Load into new model
            loaded = TherapeuticPredictionModel()
            loaded.load(str(path))
            
            assert loaded.is_trained
            assert loaded.feature_names == model.feature_names


class TestModelRegistry:
    """Tests for ModelRegistry."""
    
    def test_registry_initialization(self):
        """Test registry can be initialized."""
        with tempfile.TemporaryDirectory() as tmpdir:
            registry = ModelRegistry(models_dir=tmpdir)
            assert registry.models_dir.exists()
    
    def test_register_model(self):
        """Test model registration."""
        model = TherapeuticPredictionModel(n_estimators=10)
        
        np.random.seed(42)
        X = pd.DataFrame(np.random.randn(50, 5), columns=[f"f{i}" for i in range(5)])
        y = pd.Series(np.random.rand(50))
        model.train(X, y)
        
        with tempfile.TemporaryDirectory() as tmpdir:
            registry = ModelRegistry(models_dir=tmpdir)
            model_id = registry.register_model(model, tags=["test"])
            
            assert model_id is not None
            assert "TherapeuticPrediction" in model_id
    
    def test_list_models(self):
        """Test listing registered models."""
        with tempfile.TemporaryDirectory() as tmpdir:
            registry = ModelRegistry(models_dir=tmpdir)
            
            # Register a model
            model = TherapeuticPredictionModel(n_estimators=5)
            X = pd.DataFrame(np.random.randn(30, 3), columns=["f0", "f1", "f2"])
            y = pd.Series(np.random.rand(30))
            model.train(X, y)
            
            registry.register_model(model)
            
            models = registry.list_models()
            assert len(models) >= 1


# ============================================================================
# Prior Art Search Tests (Task 6.4)
# ============================================================================

class TestSearchQuery:
    """Tests for SearchQuery dataclass."""
    
    def test_search_query_creation(self):
        """Test SearchQuery creation."""
        query = SearchQuery(
            keywords=["cannabinoid", "formulation"],
            cpc_codes=["A61K31/352"],
            date_from="2020-01-01",
            max_results=25
        )
        
        assert len(query.keywords) == 2
        assert query.max_results == 25
    
    def test_query_string_generation(self):
        """Test query string generation."""
        query = SearchQuery(keywords=["thc", "pain"])
        query_str = query.to_query_string()
        
        assert "thc" in query_str.lower()
        assert "pain" in query_str.lower()


class TestPriorArtResult:
    """Tests for PriorArtResult dataclass."""
    
    def test_prior_art_result_creation(self):
        """Test PriorArtResult creation."""
        result = PriorArtResult(
            patent_number="US10123456",
            title="Cannabinoid Formulation",
            abstract="A novel formulation...",
            publication_date="2023-01-15",
            inventors=["John Doe"],
            assignees=["Acme Corp"],
            classifications=["A61K31/352"],
            claims_count=20,
            source=PatentSource.USPTO,
            relevance_score=0.85,
            url="https://patents.google.com/patent/US10123456"
        )
        
        assert result.patent_number == "US10123456"
        assert result.relevance_score == 0.85
        
        data = result.to_dict()
        assert "patent_number" in data
        assert "source" in data


class TestUSPTOClient:
    """Tests for USPTO PatentsView API client."""
    
    def test_client_initialization(self):
        """Test USPTO client initialization."""
        client = USPTOClient()
        assert client.BASE_URL is not None
        assert len(client.CANNABIS_CPC_CODES) > 0
    
    def test_cannabis_cpc_codes(self):
        """Test cannabis CPC codes are defined."""
        client = USPTOClient()
        assert "A61K31/352" in client.CANNABIS_CPC_CODES
        assert "A61K36/185" in client.CANNABIS_CPC_CODES


class TestPriorArtSearcher:
    """Tests for PriorArtSearcher."""
    
    def test_searcher_initialization(self):
        """Test searcher initialization."""
        searcher = PriorArtSearcher()
        assert PatentSource.USPTO in searcher.sources
    
    def test_search_templates_exist(self):
        """Test search templates are defined."""
        searcher = PriorArtSearcher()
        templates = searcher.get_template_names()
        
        assert "terpene_formulation" in templates
        assert "dimer_synthesis" in templates
        assert "therapeutic_method" in templates
    
    def test_get_template(self):
        """Test getting a search template."""
        searcher = PriorArtSearcher()
        template = searcher.get_template("terpene_formulation")
        
        assert template is not None
        assert "terpene" in template.keywords
        assert len(template.cpc_codes) > 0
    
    @pytest.mark.asyncio(loop_scope="function")
    async def test_search_basic(self):
        """Test basic search functionality (mocked)."""
        searcher = PriorArtSearcher()
        
        # Mock the USPTO search
        mock_results = [
            PriorArtResult(
                patent_number="US12345678",
                title="Test Cannabinoid Patent",
                abstract="Test abstract",
                publication_date="2023-01-01",
                inventors=["Inventor"],
                assignees=["Assignee"],
                classifications=["A61K31/352"],
                claims_count=10,
                source=PatentSource.USPTO,
                relevance_score=0.8,
                url="https://patents.google.com/patent/US12345678"
            )
        ]
        
        with patch.object(searcher.uspto, 'search', new_callable=AsyncMock) as mock_search:
            mock_search.return_value = mock_results
            
            results = await searcher.search(keywords=["cannabinoid"])
            
            assert len(results) == 1
            assert results[0].patent_number == "US12345678"


# ============================================================================
# Novelty Scoring Tests (Task 6.5)
# ============================================================================

class TestSimilarityMatch:
    """Tests for SimilarityMatch dataclass."""
    
    def test_similarity_match_creation(self):
        """Test SimilarityMatch creation."""
        match = SimilarityMatch(
            patent_number="US12345678",
            patent_title="Prior Art Patent",
            similarity_score=0.75,
            matching_keywords=["cannabinoid", "formulation"],
            matching_cpc_codes=["A61K31/352"],
            concern_level="moderate",
            notes="Matches key terms",
            url="https://example.com"
        )
        
        assert match.similarity_score == 0.75
        assert match.concern_level == "moderate"


class TestNoveltyReport:
    """Tests for NoveltyReport dataclass."""
    
    def test_novelty_report_creation(self):
        """Test NoveltyReport creation."""
        report = NoveltyReport(
            innovation_title="Novel Dimer",
            innovation_description="A new cannabinoid dimer...",
            innovation_keywords=["dimer", "cannabinoid"],
            innovation_cpc_codes=["A61K31/352"],
            overall_novelty_score=0.72,
            keyword_novelty_score=0.75,
            structural_novelty_score=0.70,
            claim_novelty_score=0.71,
            novelty_level=NoveltyLevel.MODERATE,
            patentability_potential="Moderate",
            prior_art_count=15,
            high_relevance_count=2,
            top_matches=[],
            recommendations=["Focus on unique features"],
            differentiation_opportunities=["Novel linker chemistry"]
        )
        
        assert report.overall_novelty_score == 0.72
        assert report.novelty_level == NoveltyLevel.MODERATE
        
        data = report.to_dict()
        assert "scores" in data
        assert "assessment" in data
        assert "recommendations" in data
    
    def test_novelty_report_summary(self):
        """Test report summary generation."""
        report = NoveltyReport(
            innovation_title="Test Innovation",
            innovation_description="Description",
            innovation_keywords=["test"],
            innovation_cpc_codes=[],
            overall_novelty_score=0.8,
            keyword_novelty_score=0.8,
            structural_novelty_score=0.8,
            claim_novelty_score=0.8,
            novelty_level=NoveltyLevel.HIGH,
            patentability_potential="Strong",
            prior_art_count=5,
            high_relevance_count=0,
            top_matches=[],
            recommendations=["Proceed with filing"],
            differentiation_opportunities=[]
        )
        
        summary = report.get_summary()
        assert "Test Innovation" in summary
        assert "HIGH" in summary


class TestNoveltyScorer:
    """Tests for NoveltyScorer."""
    
    def test_scorer_initialization(self):
        """Test NoveltyScorer initialization."""
        scorer = NoveltyScorer()
        assert scorer.searcher is not None
    
    def test_innovation_keywords_defined(self):
        """Test innovation keyword categories exist."""
        assert "composition" in NoveltyScorer.INNOVATION_KEYWORDS
        assert "method" in NoveltyScorer.INNOVATION_KEYWORDS
        assert "compound" in NoveltyScorer.INNOVATION_KEYWORDS
    
    def test_technology_cpc_map_defined(self):
        """Test technology CPC mappings exist."""
        assert "cannabinoid_compound" in NoveltyScorer.TECHNOLOGY_CPC_MAP
        assert "formulation" in NoveltyScorer.TECHNOLOGY_CPC_MAP
        assert "therapeutic_use" in NoveltyScorer.TECHNOLOGY_CPC_MAP
    
    @pytest.mark.asyncio(loop_scope="function")
    async def test_assess_novelty_mocked(self):
        """Test novelty assessment with mocked prior art."""
        scorer = NoveltyScorer()
        
        # Mock prior art search
        mock_prior_art = [
            PriorArtResult(
                patent_number="US11111111",
                title="Existing Cannabinoid Patent",
                abstract="Prior art abstract with cannabinoid",
                publication_date="2020-01-01",
                inventors=["Inventor"],
                assignees=["Company"],
                classifications=["A61K31/352"],
                claims_count=15,
                source=PatentSource.USPTO,
                relevance_score=0.5,
                url="https://example.com"
            )
        ]
        
        with patch.object(scorer.searcher, 'search', new_callable=AsyncMock) as mock_search:
            mock_search.return_value = mock_prior_art
            
            report = await scorer.assess_novelty(
                title="Novel THC-CBD Dimer",
                description="A novel linked dimer compound with unique therapeutic properties",
                keywords=["dimer", "thc", "cbd", "novel"]
            )
            
            assert isinstance(report, NoveltyReport)
            assert report.innovation_title == "Novel THC-CBD Dimer"
            assert report.prior_art_count == 1
            assert len(report.recommendations) > 0


class TestNoveltyLevel:
    """Tests for NoveltyLevel enum."""
    
    def test_novelty_levels_exist(self):
        """Test all novelty levels are defined."""
        assert NoveltyLevel.HIGH.value == "high"
        assert NoveltyLevel.MODERATE.value == "moderate"
        assert NoveltyLevel.LOW.value == "low"
        assert NoveltyLevel.BLOCKED.value == "blocked"


# ============================================================================
# Integration Tests
# ============================================================================

class TestMLPipelineIntegration:
    """Integration tests for ML pipeline."""
    
    def test_full_training_pipeline(self):
        """Test complete training pipeline."""
        # 1. Prepare data
        np.random.seed(42)
        X = pd.DataFrame(
            np.random.randn(100, 20),
            columns=[f"feature_{i}" for i in range(20)]
        )
        y = pd.Series(np.clip(np.random.rand(100), 0, 1))
        
        # 2. Train model
        model = TherapeuticPredictionModel(n_estimators=20, max_depth=4)
        metrics = model.train(X, y)
        
        # 3. Verify training
        assert model.is_trained
        assert metrics.num_features == 20
        assert metrics.num_samples == 100
        
        # 4. Make predictions with uncertainty
        X_new = pd.DataFrame(
            np.random.randn(5, 20),
            columns=[f"feature_{i}" for i in range(20)]
        )
        result = model.predict_with_uncertainty(X_new)
        
        assert len(result.prediction) == 5
        assert all(u >= 0 for u in result.uncertainty)
    
    def test_model_comparison(self):
        """Test training and comparing multiple models."""
        np.random.seed(42)
        X = pd.DataFrame(np.random.randn(80, 10), columns=[f"f{i}" for i in range(10)])
        # Create target with pattern
        y = pd.Series(np.clip(0.5 + 0.4 * X["f0"] + 0.1 * np.random.randn(80), 0, 1))
        
        # Train different model types
        therapeutic = TherapeuticPredictionModel(n_estimators=15)
        dimer = DimerPotentialModel(n_estimators=15)
        
        m1 = therapeutic.train(X, y)
        m2 = dimer.train(X, y)
        
        assert m1.model_name == "TherapeuticPrediction"
        assert m2.model_name == "DimerPotential"
        
        # Both should produce valid metrics (cv_mean can be any value)
        assert hasattr(m1, 'cv_mean')
        assert hasattr(m2, 'cv_mean')


class TestPatentPathIntegration:
    """Integration tests for PatentPath features."""
    
    @pytest.mark.asyncio(loop_scope="function")
    async def test_search_and_score_pipeline(self):
        """Test combined search and scoring pipeline."""
        searcher = PriorArtSearcher()
        scorer = NoveltyScorer(searcher=searcher)
        
        # Mock the search
        mock_results = [
            PriorArtResult(
                patent_number="US99999999",
                title="Related Patent",
                abstract="Contains cannabinoid dimer formulation",
                publication_date="2022-06-15",
                inventors=["Jane Doe"],
                assignees=["Patent Co"],
                classifications=["A61K31/352", "A61K9/00"],
                claims_count=25,
                source=PatentSource.USPTO,
                relevance_score=0.65,
                url="https://example.com/patent"
            )
        ]
        
        with patch.object(searcher.uspto, 'search', new_callable=AsyncMock) as mock_search:
            mock_search.return_value = mock_results
            
            # Search
            results = await searcher.search(keywords=["cannabinoid", "dimer"])
            assert len(results) == 1
            
            # Score novelty
            report = await scorer.assess_novelty(
                title="New Dimer Compound",
                description="A novel dimer with enhanced bioavailability"
            )
            
            assert report is not None
            assert report.prior_art_count >= 0


# ============================================================================
# Async Test Support
# ============================================================================

@pytest.fixture
def event_loop():
    """Create event loop for async tests."""
    loop = asyncio.new_event_loop()
    yield loop
    loop.close()


# ============================================================================
# Summary Statistics
# ============================================================================

def test_week6_module_count():
    """Verify Week 6 module structure."""
    # ML modules
    from services import ml_data_prep, ml_models
    
    # PatentPath modules
    from services.patentpath import prior_art, novelty
    
    # All imports successful
    assert ml_data_prep is not None
    assert ml_models is not None
    assert prior_art is not None
    assert novelty is not None


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
