"""Comprehensive tests for the whole-plant cannabis analysis engine."""

from __future__ import annotations

import json
import sqlite3
import time
from typing import Dict, List

import pytest
from fastapi.testclient import TestClient
from pydantic import ValidationError

from backend.database import get_db_connection
from backend.engines.whole_plant import WholePlantAnalyzer
from backend.models.compound_profile import CompoundProfile, WholePlantAnalysisResult
from backend.services.terpene_synergy import TerpeneSynergyScorer


# ============================================================================
# TEST DATABASE SETUP
# ============================================================================


def create_test_db() -> sqlite3.Connection:
    """Provision an in-memory database with representative data."""
    conn = sqlite3.connect(":memory:")
    cursor = conn.cursor()

    cursor.execute(
        """
        CREATE TABLE neurobotanica_compounds (
            compound_id TEXT PRIMARY KEY,
            compound_name TEXT,
            therapeutic_targets TEXT,
            primary_mechanisms TEXT
        )
        """
    )
    cursor.executemany(
        "INSERT INTO neurobotanica_compounds VALUES (?, ?, ?, ?)",
        [
            ("cbd", "CBD", json.dumps(["anxiety", "epilepsy", "inflammation"]), "5ht1a"),
            ("thc", "THC", json.dumps(["chronic_pain", "appetite_loss"]), "cb1"),
            ("cbg", "CBG", json.dumps(["anxiety", "inflammation"]), "cb2"),
        ],
    )

    synthetic_rows: List[tuple[str, str, str, str]] = []
    for idx in range(1, 185):
        compound_id = f"cmp_{idx:03d}"
        compound_name = f"Compound {idx}"
        base_condition = "anxiety" if idx % 2 == 0 else "chronic_pain"
        synthetic_rows.append(
            (
                compound_id,
                compound_name,
                json.dumps([base_condition]),
                "cb1",
            )
        )
    cursor.executemany(
        "INSERT INTO neurobotanica_compounds VALUES (?, ?, ?, ?)", synthetic_rows
    )

    cursor.execute(
        """
        CREATE TABLE neurobotanica_synergy_predictions (
            synergy_id TEXT PRIMARY KEY,
            compound_a_id TEXT,
            compound_b_id TEXT,
            synergy_score REAL,
            mechanism TEXT,
            therapeutic_context TEXT,
            confidence_score REAL,
            terpene TEXT,
            condition TEXT,
            correlation_strength REAL,
            evidence_weight REAL,
            citations TEXT,
            prediction_type TEXT
        )
        """
    )
    cursor.executemany(
        """
        INSERT INTO neurobotanica_synergy_predictions
        (synergy_id, compound_a_id, compound_b_id, synergy_score, mechanism,
         therapeutic_context, confidence_score, terpene, condition,
         correlation_strength, evidence_weight, citations, prediction_type)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        [
            (
                "cbd_linalool_anxiety",
                "cbd",
                "linalool",
                0.84,
                "serotonergic modulation",
                "anxiety,sleep",
                0.78,
                "linalool",
                "anxiety",
                0.82,
                0.74,
                "PMID:12345",
                "terpene_condition",
            ),
            (
                "cbd_myrcene_sleep",
                "cbd",
                "myrcene",
                0.76,
                "sedative",
                "sleep",
                0.7,
                "myrcene",
                "sleep",
                0.9,
                0.72,
                "PMID:67890",
                "terpene_condition",
            ),
            (
                "thc_b_caryophyllene_pain",
                "thc",
                "beta_caryophyllene",
                0.88,
                "cb2 agonism",
                "chronic_pain",
                0.81,
                "beta_caryophyllene",
                "chronic_pain",
                0.9,
                0.8,
                "PMID:11111",
                "terpene_condition",
            ),
        ],
    )

    cursor.execute(
        """
        CREATE TABLE neurobotanica_clinical_studies (
            study_id TEXT PRIMARY KEY,
            study_type TEXT,
            condition TEXT,
            intervention TEXT,
            key_findings TEXT,
            citation TEXT,
            confidence_score REAL,
            sample_size INTEGER,
            publication_year INTEGER
        )
        """
    )
    cursor.executemany(
        "INSERT INTO neurobotanica_clinical_studies VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
        [
            (
                "study_anxiety_1",
                "RCT",
                "anxiety",
                "CBD 300mg",
                "Significant reduction in anxiety scores",
                "J Clin 2024",
                0.82,
                72,
                2024,
            ),
            (
                "study_sleep_1",
                "Observational",
                "sleep",
                "CBD + myrcene",
                "Improved sleep latency",
                "Sleep Med 2023",
                0.7,
                60,
                2023,
            ),
        ],
    )

    conn.commit()
    return conn


# ============================================================================
# FIXTURES
# ============================================================================


@pytest.fixture
def db_connection() -> sqlite3.Connection:
    conn = create_test_db()
    yield conn
    conn.close()


@pytest.fixture
def sample_anxiety_profile() -> CompoundProfile:
    return CompoundProfile(
        plant_id="anxiety_test_001",
        cannabinoid_profile={"thc": 8.5, "cbd": 15.0, "cbg": 1.2, "cbn": 0.3},
        terpene_profile={"myrcene": 1.0, "limonene": 0.8, "pinene": 0.4, "linalool": 1.2},
        target_condition="anxiety",
    )


@pytest.fixture
def sample_pain_profile() -> CompoundProfile:
    return CompoundProfile(
        plant_id="pain_test_001",
        cannabinoid_profile={"thc": 18.5, "cbd": 2.0, "cbg": 0.8, "cbn": 1.5},
        terpene_profile={"myrcene": 1.5, "beta_caryophyllene": 1.8, "pinene": 0.5, "linalool": 0.3},
        target_condition="chronic_pain",
    )


@pytest.fixture
def high_thc_profile() -> CompoundProfile:
    return CompoundProfile(
        plant_id="high_thc_test",
        cannabinoid_profile={"thc": 28.0, "cbd": 0.5, "cbg": 0.2, "cbn": 0.1},
        terpene_profile={"myrcene": 0.8, "limonene": 0.4, "pinene": 0.2, "linalool": 0.1},
        target_condition="anxiety",
    )


@pytest.fixture
def api_test_client(db_connection: sqlite3.Connection) -> TestClient:
    from backend.main import app

    def _override_db():
        return db_connection

    app.dependency_overrides[get_db_connection] = _override_db
    client = TestClient(app)
    yield client
    app.dependency_overrides.pop(get_db_connection, None)


# ============================================================================
# MODEL TESTS
# ============================================================================


class TestCompoundProfile:
    def test_valid_profile_creation(self, sample_anxiety_profile: CompoundProfile) -> None:
        assert sample_anxiety_profile.plant_id == "anxiety_test_001"
        assert sample_anxiety_profile.cannabinoid_profile["cbd"] == 15.0
        assert sample_anxiety_profile.target_condition == "anxiety"

    def test_profile_dict_conversion(self, sample_anxiety_profile: CompoundProfile) -> None:
        profile_dict = sample_anxiety_profile.model_dump()
        assert "plant_id" in profile_dict
        assert profile_dict["cannabinoid_profile"]["thc"] == 8.5

    def test_profile_json_serialization(self, sample_anxiety_profile: CompoundProfile) -> None:
        json_str = sample_anxiety_profile.model_dump_json()
        assert "anxiety_test_001" in json_str
        assert "cannabinoid_profile" in json_str


class TestWholePlantAnalysisResult:
    def test_confidence_range_validation(self) -> None:
        result = WholePlantAnalysisResult(
            synergy_score=0.8,
            confidence=0.75,
            recommended_conditions=["anxiety"],
            entourage_compounds=[],
            clinical_evidence=[],
        )
        assert 0.5 <= result.confidence <= 0.9

        with pytest.raises(ValidationError):
            WholePlantAnalysisResult(
                synergy_score=0.8,
                confidence=0.3,
                recommended_conditions=["anxiety"],
                entourage_compounds=[],
                clinical_evidence=[],
            )

        with pytest.raises(ValidationError):
            WholePlantAnalysisResult(
                synergy_score=0.8,
                confidence=0.95,
                recommended_conditions=["anxiety"],
                entourage_compounds=[],
                clinical_evidence=[],
            )

    def test_synergy_score_range(self) -> None:
        result = WholePlantAnalysisResult(
            synergy_score=0.85,
            confidence=0.7,
            recommended_conditions=[],
            entourage_compounds=[],
            clinical_evidence=[],
        )
        assert 0.0 <= result.synergy_score <= 1.0


# ============================================================================
# TERPENE SYNERGY TESTS
# ============================================================================


class TestTerpeneSynergyScorer:
    def test_scorer_initialization(self, db_connection: sqlite3.Connection) -> None:
        scorer = TerpeneSynergyScorer(db_connection)
        assert scorer.db is not None

    def test_anxiety_terpene_scoring(self, db_connection: sqlite3.Connection) -> None:
        scorer = TerpeneSynergyScorer(db_connection)
        terpene_profile = {"linalool": 1.2, "myrcene": 1.0, "limonene": 0.8}
        synergy_score, confidence = scorer.calculate_synergy_score(terpene_profile, "anxiety")
        assert 0.0 <= synergy_score <= 1.0
        assert 0.5 <= confidence <= 0.9
        assert synergy_score > 0.6

    def test_pain_terpene_scoring(self, db_connection: sqlite3.Connection) -> None:
        scorer = TerpeneSynergyScorer(db_connection)
        terpene_profile = {"beta_caryophyllene": 1.8, "myrcene": 1.5}
        synergy_score, _ = scorer.calculate_synergy_score(terpene_profile, "chronic_pain")
        assert synergy_score > 0.7

    def test_empty_terpene_profile(self, db_connection: sqlite3.Connection) -> None:
        scorer = TerpeneSynergyScorer(db_connection)
        synergy_score, confidence = scorer.calculate_synergy_score({}, "anxiety")
        assert synergy_score < 0.3
        assert confidence == 0.5

    def test_top_terpenes_for_condition(self, db_connection: sqlite3.Connection) -> None:
        scorer = TerpeneSynergyScorer(db_connection)
        top_terpenes = scorer.get_top_terpenes_for_condition("anxiety", top_n=3)
        assert len(top_terpenes) == 3
        assert any(entry["terpene"] == "linalool" for entry in top_terpenes)


# ============================================================================
# WHOLE-PLANT ANALYZER TESTS
# ============================================================================


class TestWholePlantAnalyzer:
    def test_analyzer_initialization(self, db_connection: sqlite3.Connection) -> None:
        analyzer = WholePlantAnalyzer(db_connection)
        assert analyzer.db is not None
        assert analyzer.terpene_scorer is not None

    def test_anxiety_profile_analysis(
        self, db_connection: sqlite3.Connection, sample_anxiety_profile: CompoundProfile
    ) -> None:
        analyzer = WholePlantAnalyzer(db_connection)
        result = analyzer.analyze(sample_anxiety_profile)
        assert isinstance(result, WholePlantAnalysisResult)
        assert 0.0 <= result.synergy_score <= 1.0
        assert 0.5 <= result.confidence <= 0.9
        assert "anxiety" in [cond.lower() for cond in result.recommended_conditions]
        assert result.synergy_score > 0.7

    def test_pain_profile_analysis(
        self, db_connection: sqlite3.Connection, sample_pain_profile: CompoundProfile
    ) -> None:
        analyzer = WholePlantAnalyzer(db_connection)
        result = analyzer.analyze(sample_pain_profile)
        assert "chronic_pain" in result.recommended_conditions
        assert result.synergy_score > 0.7

    def test_high_thc_anxiety_warning(
        self, db_connection: sqlite3.Connection, high_thc_profile: CompoundProfile
    ) -> None:
        analyzer = WholePlantAnalyzer(db_connection)
        result = analyzer.analyze(high_thc_profile)
        assert result.synergy_score < 0.6

    def test_entourage_effect_detection(
        self, db_connection: sqlite3.Connection, sample_anxiety_profile: CompoundProfile
    ) -> None:
        analyzer = WholePlantAnalyzer(db_connection)
        result = analyzer.analyze(sample_anxiety_profile)
        assert len(result.entourage_compounds) > 0
        for compound in result.entourage_compounds:
            assert "synergy_score" in compound
            assert "compound_a" in compound
            assert "compound_b" in compound

    def test_clinical_evidence_retrieval(
        self, db_connection: sqlite3.Connection, sample_anxiety_profile: CompoundProfile
    ) -> None:
        analyzer = WholePlantAnalyzer(db_connection)
        result = analyzer.analyze(sample_anxiety_profile)
        assert len(result.clinical_evidence) > 0
        for evidence in result.clinical_evidence:
            assert "study" in evidence or "study_id" in evidence
            assert "citation" in evidence
            assert "finding" in evidence or "summary" in evidence

    def test_multi_condition_recommendations(
        self, db_connection: sqlite3.Connection, sample_anxiety_profile: CompoundProfile
    ) -> None:
        analyzer = WholePlantAnalyzer(db_connection)
        result = analyzer.analyze(sample_anxiety_profile)
        assert len(result.recommended_conditions) >= 2
        overlap = set(result.recommended_conditions) & {"anxiety", "sleep", "ptsd"}
        assert len(overlap) >= 2


# ============================================================================
# PERFORMANCE TESTS
# ============================================================================


class TestPerformance:
    def test_analysis_speed(
        self, db_connection: sqlite3.Connection, sample_anxiety_profile: CompoundProfile
    ) -> None:
        analyzer = WholePlantAnalyzer(db_connection)
        start = time.time()
        analyzer.analyze(sample_anxiety_profile)
        duration = time.time() - start
        assert duration < 2.0, f"Analysis took {duration:.2f}s (target: <2s)"

    def test_batch_analysis_speed(self, db_connection: sqlite3.Connection) -> None:
        analyzer = WholePlantAnalyzer(db_connection)
        profiles = [
            CompoundProfile(
                plant_id=f"batch_{idx}",
                cannabinoid_profile={"thc": 10.0 + idx, "cbd": 5.0},
                terpene_profile={"myrcene": 1.0, "limonene": 0.5},
                target_condition="anxiety",
            )
            for idx in range(10)
        ]
        start = time.time()
        results = [analyzer.analyze(profile) for profile in profiles]
        avg_time = (time.time() - start) / len(profiles)
        assert avg_time < 2.0, f"Average analysis time: {avg_time:.2f}s (target: <2s)"
        assert len(results) == 10


# ============================================================================
# VALIDATION TESTS (vs. Budtender App)
# ============================================================================


class TestBudtenderValidation:
    @pytest.fixture
    def budtender_test_cases(self) -> List[Dict[str, object]]:
        return [
            {
                "profile": CompoundProfile(
                    plant_id="budtender_001",
                    cannabinoid_profile={"thc": 15.0, "cbd": 10.0, "cbg": 1.0, "cbn": 0.5},
                    terpene_profile={"myrcene": 1.2, "linalool": 0.8},
                    target_condition="anxiety",
                ),
                "expected_synergy": 0.58,
                "expected_confidence": 0.75,
            }
        ]

    def test_synergy_score_accuracy(
        self, db_connection: sqlite3.Connection, budtender_test_cases: List[Dict[str, object]]
    ) -> None:
        analyzer = WholePlantAnalyzer(db_connection)
        matches = 0
        for case in budtender_test_cases:
            result = analyzer.analyze(case["profile"])
            expected = case["expected_synergy"]
            variance = abs(result.synergy_score - expected) / expected
            if variance < 0.05:
                matches += 1
        accuracy = matches / len(budtender_test_cases)
        assert accuracy >= 0.95


# ============================================================================
# EDGE CASE TESTS
# ============================================================================


class TestEdgeCases:
    def test_zero_cannabinoids(self, db_connection: sqlite3.Connection) -> None:
        profile = CompoundProfile(
            plant_id="zero_test",
            cannabinoid_profile={"thc": 0.0, "cbd": 0.0},
            terpene_profile={"myrcene": 1.0},
            target_condition="anxiety",
        )
        analyzer = WholePlantAnalyzer(db_connection)
        result = analyzer.analyze(profile)
        assert isinstance(result, WholePlantAnalysisResult)
        assert result.synergy_score < 0.3

    def test_unknown_condition(self, db_connection: sqlite3.Connection) -> None:
        profile = CompoundProfile(
            plant_id="unknown_test",
            cannabinoid_profile={"thc": 10.0, "cbd": 5.0},
            terpene_profile={"myrcene": 1.0},
            target_condition="unknown_condition_xyz",
        )
        analyzer = WholePlantAnalyzer(db_connection)
        result = analyzer.analyze(profile)
        assert result.confidence >= 0.5

    def test_all_184_cannabinoids_supported(self, db_connection: sqlite3.Connection) -> None:
        analyzer = WholePlantAnalyzer(db_connection)
        cursor = db_connection.cursor()
        cursor.execute("SELECT LOWER(compound_id) FROM neurobotanica_compounds")
        compound_ids = {row[0] for row in cursor.fetchall()}
        assert len(compound_ids) >= 184
        assert compound_ids.issubset(set(analyzer.compound_metadata.keys()))

    def test_22_conditions_supported(self, db_connection: sqlite3.Connection) -> None:
        conditions = [
            "anxiety",
            "chronic_pain",
            "sleep",
            "ptsd",
            "inflammation",
            "depression",
            "nausea",
            "appetite_loss",
            "seizures",
            "glaucoma",
            "fibromyalgia",
            "migraine",
            "arthritis",
            "neuropathy",
            "epilepsy",
            "ibs",
            "crohns",
            "autism",
            "adhd",
            "menstrual_pain",
            "cancer_pain",
            "parkinsons",
        ]
        analyzer = WholePlantAnalyzer(db_connection)
        template = CompoundProfile(
            plant_id="condition_test",
            cannabinoid_profile={"thc": 10.0, "cbd": 5.0},
            terpene_profile={"myrcene": 1.0},
            target_condition="anxiety",
        )
        for condition in conditions:
            profile = template.model_copy(update={"target_condition": condition})
            result = analyzer.analyze(profile)
            assert isinstance(result, WholePlantAnalysisResult)


# ============================================================================
# INTEGRATION TESTS (API Endpoint)
# ============================================================================


class TestAPIEndpoint:
    def test_whole_plant_endpoint_success(
        self,
        api_test_client: TestClient,
        sample_anxiety_profile: CompoundProfile,
    ) -> None:
        response = api_test_client.post(
            "/api/v1/analysis/whole-plant",
            json=sample_anxiety_profile.model_dump(),
        )
        assert response.status_code == 200
        data = response.json()
        assert "synergy_score" in data
        assert "confidence" in data
        assert "recommended_conditions" in data

    def test_whole_plant_endpoint_invalid_input(self, api_test_client: TestClient) -> None:
        response = api_test_client.post(
            "/api/v1/analysis/whole-plant",
            json={"invalid": "data"},
        )
        assert response.status_code == 422