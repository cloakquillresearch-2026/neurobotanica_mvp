"""
NeuroBotanica Unit Tests - Models
=================================

Test database models: ClinicalStudy, Cannabinoid, Patient, Treatment

Run tests:
    pytest backend/tests/test_models.py -v
"""
import pytest
from datetime import datetime
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

# Import models
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

from models.database import Base
from models.study import ClinicalStudy, StudyType, EvidenceGrade
from models.compound import Cannabinoid
from models.patient import Patient, TreatmentResponse, ConsentStatus
from models.treatment import Treatment, TreatmentTemplate


# Create test database
TEST_DATABASE_URL = "sqlite:///:memory:"
test_engine = create_engine(TEST_DATABASE_URL, connect_args={"check_same_thread": False})
TestSession = sessionmaker(autocommit=False, autoflush=False, bind=test_engine)


@pytest.fixture(scope="function")
def db_session():
    """Create a fresh database session for each test."""
    Base.metadata.create_all(bind=test_engine)
    session = TestSession()
    try:
        yield session
    finally:
        session.close()
        Base.metadata.drop_all(bind=test_engine)


class TestClinicalStudyModel:
    """Tests for ClinicalStudy model."""
    
    def test_create_study(self, db_session):
        """Test creating a clinical study."""
        study = ClinicalStudy(
            study_id="NB-001",
            condition="chronic_pain",
            study_type="RCT",
            study_title="Test Study for Chronic Pain",
            cannabinoid="CBD",
            sample_size="50",
            year=2023,
            evidence_grade="Level 1",
            confidence_weight=0.85
        )
        
        db_session.add(study)
        db_session.commit()
        
        retrieved = db_session.query(ClinicalStudy).filter_by(study_id="NB-001").first()
        assert retrieved is not None
        assert retrieved.condition == "chronic_pain"
        assert retrieved.cannabinoid == "CBD"
        assert retrieved.confidence_weight == 0.85
    
    def test_study_pharmacology_package(self, db_session):
        """Test generating pharmacology package from study."""
        study = ClinicalStudy(
            study_id="NB-002",
            condition="epilepsy",
            study_type="RCT",
            cannabinoid="CBD",
            mechanism_of_action="GABA-A receptor modulation",
            receptor_targets=["GABA-A", "5-HT1A"],
            sample_size="100",
            effect_size="50% seizure reduction",
            adverse_events="Drowsiness, fatigue"
        )
        
        db_session.add(study)
        db_session.commit()
        
        package = study.to_pharmacology_package()
        
        assert package["study_id"] == "NB-002"
        assert package["cannabinoid"] == "CBD"
        assert package["mechanism_of_action"] == "GABA-A receptor modulation"
        assert "GABA-A" in package["receptor_targets"]
        assert package["safety_profile"]["adverse_events"] == "Drowsiness, fatigue"
    
    def test_study_efficacy_comparison(self, db_session):
        """Test generating efficacy comparison from study."""
        study = ClinicalStudy(
            study_id="NB-003",
            condition="anxiety",
            study_type="RCT",
            cannabinoid="CBD",
            dosage="25mg",
            effect_size="d=0.8"
        )
        
        db_session.add(study)
        db_session.commit()
        
        comparison = study.to_efficacy_comparison()
        
        assert comparison["condition"] == "anxiety"
        assert comparison["intervention"]["cannabinoid"] == "CBD"
        assert comparison["intervention"]["dosage"] == "25mg"


class TestCannabinoidModel:
    """Tests for Cannabinoid model."""
    
    def test_create_cannabinoid(self, db_session):
        """Test creating a cannabinoid compound."""
        cannabinoid = Cannabinoid(
            name="Cannabidiol",
            abbreviation="CBD",
            smiles="CCCCCC1=CC(=C(C(=C1)O)C2C=C(CCC2C(=C)C)C)O",
            molecular_weight=314.47,
            logp=5.79,
            is_natural=True,
            is_dimeric=False,
            compound_class="Phytocannabinoid"
        )
        
        db_session.add(cannabinoid)
        db_session.commit()
        
        retrieved = db_session.query(Cannabinoid).filter_by(name="Cannabidiol").first()
        assert retrieved is not None
        assert retrieved.abbreviation == "CBD"
        assert retrieved.molecular_weight == 314.47
        assert retrieved.is_natural is True
    
    def test_cannabinoid_cmc_profile(self, db_session):
        """Test generating CMC profile from cannabinoid."""
        cannabinoid = Cannabinoid(
            name="Delta-9-THC",
            abbreviation="THC",
            smiles="CCCCCC1=CC(=C2C3C=C(CCC3C(OC2=C1)(C)C)C)O",
            inchi="InChI=1S/C21H30O2/...",
            inchi_key="CYQFCXCEBYINGO-IAGOWNOFSA-N",
            molecular_weight=314.47,
            logp=7.68,
            fda_approved=True,
            fda_drug_name="Marinol",
            schedule_classification="III"
        )
        
        db_session.add(cannabinoid)
        db_session.commit()
        
        cmc = cannabinoid.to_cmc_profile()
        
        assert cmc["compound_identity"]["name"] == "Delta-9-THC"
        assert cmc["compound_identity"]["abbreviation"] == "THC"
        assert cmc["physicochemical_properties"]["molecular_weight"] == 314.47
        assert cmc["classification"]["fda_approved"] is True
        assert cmc["classification"]["schedule"] == "III"
    
    def test_cannabinoid_3d_fields(self, db_session):
        """Test 3D conformer fields on cannabinoid."""
        cannabinoid = Cannabinoid(
            name="Test Compound",
            abbreviation="TC",
            smiles="CCO",
            has_conformers=True,
            conformer_generation_method="ETKDG",
            num_conformers_generated=50,
            conformer_metadata={
                "lowest_energy_kcal_mol": -10.5,
                "energy_range_kcal_mol": 5.2
            },
            rdkit_descriptors_3d={
                "asphericity": 0.15,
                "radius_of_gyration": 2.5
            }
        )
        
        db_session.add(cannabinoid)
        db_session.commit()
        
        retrieved = db_session.query(Cannabinoid).filter_by(name="Test Compound").first()
        assert retrieved.has_conformers is True
        assert retrieved.num_conformers_generated == 50
        assert retrieved.conformer_metadata["lowest_energy_kcal_mol"] == -10.5
        assert retrieved.rdkit_descriptors_3d["asphericity"] == 0.15


class TestPatientModel:
    """Tests for Patient model."""
    
    def test_create_patient(self, db_session):
        """Test creating a patient profile."""
        patient = Patient(
            patient_hash="abc123hash",
            profile_code="NB-00001",
            age_range="35-44",
            biological_sex="female",
            state_code="NV",
            primary_conditions=["chronic_pain", "insomnia"],
            cannabinoid_experience_level="occasional",
            primary_treatment_goal="pain relief",
            omnipath_consent_status="full_consent"
        )
        
        db_session.add(patient)
        db_session.commit()
        
        retrieved = db_session.query(Patient).filter_by(profile_code="NB-00001").first()
        assert retrieved is not None
        assert retrieved.age_range == "35-44"
        assert "chronic_pain" in retrieved.primary_conditions
    
    def test_patient_recommendation_profile(self, db_session):
        """Test generating recommendation profile from patient."""
        patient = Patient(
            patient_hash="xyz456hash",
            profile_code="NB-00002",
            age_range="25-34",
            biological_sex="male",
            weight_category="normal",
            state_code="CA",
            jurisdiction_type="recreational",
            primary_conditions=["anxiety"],
            cannabinoid_experience_level="regular",
            thc_tolerance="medium",
            cbd_preference=5.0,  # 1:5 THC:CBD ratio
            preferred_delivery_methods=["sublingual", "edible"],
            primary_treatment_goal="anxiety reduction"
        )
        
        db_session.add(patient)
        db_session.commit()
        
        profile = patient.to_recommendation_profile()
        
        assert profile["patient_code"] == "NB-00002"
        assert profile["demographics"]["age_range"] == "25-34"
        assert profile["cannabinoid_profile"]["thc_tolerance"] == "medium"
        assert "sublingual" in profile["cannabinoid_profile"]["delivery_preferences"]
        assert profile["treatment_goals"]["primary"] == "anxiety reduction"
    
    def test_patient_completeness_calculation(self, db_session):
        """Test profile completeness calculation."""
        # Complete patient
        complete_patient = Patient(
            patient_hash="comp123",
            profile_code="NB-00003",
            age_range="45-54",
            biological_sex="female",
            weight_category="normal",
            state_code="NV",
            primary_conditions=["chronic_pain"],
            cannabinoid_experience_level="daily",
            primary_treatment_goal="pain management",
            current_medications=["NSAID"],
            thc_tolerance="high",
            preferred_delivery_methods=["inhalation"],
            medical_history_summary="Previous spine injury"
        )
        
        db_session.add(complete_patient)
        db_session.commit()
        
        completeness = complete_patient.calculate_completeness()
        assert completeness >= 0.9  # Should be near 100%
        
        # Incomplete patient
        incomplete_patient = Patient(
            patient_hash="incomp456",
            profile_code="NB-00004",
            age_range="25-34"
        )
        
        db_session.add(incomplete_patient)
        db_session.commit()
        
        completeness = incomplete_patient.calculate_completeness()
        assert completeness < 0.5  # Should be low


class TestTreatmentModel:
    """Tests for Treatment model."""
    
    def test_create_treatment(self, db_session):
        """Test creating a treatment recommendation."""
        # First create a patient
        patient = Patient(
            patient_hash="pat789",
            profile_code="NB-00005",
            age_range="35-44",
            biological_sex="male",
            state_code="NV",
            primary_conditions=["anxiety"]
        )
        db_session.add(patient)
        db_session.commit()
        
        # Create treatment
        treatment = Treatment(
            treatment_code="NB-TRT-00001",
            patient_id=patient.id,
            primary_condition="anxiety",
            condition_severity="moderate",
            primary_cannabinoid_name="CBD",
            primary_cannabinoid_dosage="25mg",
            cannabinoid_ratio="CBD only",
            delivery_method="sublingual",
            initial_dose="10mg",
            maintenance_dose="25mg",
            max_daily_dose="100mg",
            dosing_frequency="BID",
            duration_weeks=8,
            evidence_strength="strong",
            schedule_classification="unscheduled",
            status="proposed"
        )
        
        db_session.add(treatment)
        db_session.commit()
        
        retrieved = db_session.query(Treatment).filter_by(treatment_code="NB-TRT-00001").first()
        assert retrieved is not None
        assert retrieved.primary_cannabinoid_name == "CBD"
        assert retrieved.delivery_method == "sublingual"
    
    def test_treatment_patient_summary(self, db_session):
        """Test generating patient summary from treatment."""
        patient = Patient(
            patient_hash="sum123",
            profile_code="NB-00006",
            state_code="NV"
        )
        db_session.add(patient)
        db_session.commit()
        
        treatment = Treatment(
            treatment_code="NB-TRT-00002",
            patient_id=patient.id,
            primary_condition="chronic_pain",
            primary_cannabinoid_name="THC:CBD",
            primary_cannabinoid_dosage="5mg:20mg",
            cannabinoid_ratio="THC:CBD 1:4",
            delivery_method="oral",
            dosing_frequency="TID",
            duration_weeks=12,
            time_to_effect_days=14,
            response_probability=0.75,
            warnings=["Do not drive", "May cause drowsiness"],
            monitoring_requirements=["Monthly check-in"]
        )
        
        db_session.add(treatment)
        db_session.commit()
        
        summary = treatment.to_patient_summary()
        
        assert summary["condition"] == "chronic_pain"
        assert summary["recommendation"]["ratio"] == "THC:CBD 1:4"
        assert summary["schedule"]["duration_weeks"] == 12
        assert summary["expected_outcomes"]["response_probability"] == "75%"
        assert "Do not drive" in summary["important_warnings"]
    
    def test_treatment_fda_documentation(self, db_session):
        """Test generating FDA documentation from treatment."""
        patient = Patient(
            patient_hash="fda123",
            profile_code="NB-00007",
            state_code="NV"
        )
        db_session.add(patient)
        db_session.commit()
        
        treatment = Treatment(
            treatment_code="NB-TRT-00003",
            patient_id=patient.id,
            primary_condition="epilepsy",
            condition_severity="severe",
            primary_cannabinoid_name="CBD",
            initial_dose="2.5mg/kg",
            maintenance_dose="10mg/kg",
            max_daily_dose="20mg/kg",
            delivery_method="oral",
            dosing_frequency="BID",
            duration_weeks=52,
            schedule_classification="unscheduled",
            fda_drug_reference="Epidiolex",
            evidence_citations=["NB-EPI-001", "NB-EPI-002"],
            evidence_strength="strong"
        )
        
        db_session.add(treatment)
        db_session.commit()
        
        fda_doc = treatment.to_fda_documentation()
        
        assert fda_doc["treatment_protocol_id"] == "NB-TRT-00003"
        assert fda_doc["indication"]["primary_condition"] == "epilepsy"
        assert fda_doc["dosing_regimen"]["maintenance_dose"] == "10mg/kg"
        assert fda_doc["regulatory_classification"]["fda_reference_drug"] == "Epidiolex"


class TestTreatmentResponse:
    """Tests for TreatmentResponse model."""
    
    def test_create_response(self, db_session):
        """Test creating a treatment response."""
        patient = Patient(
            patient_hash="resp123",
            profile_code="NB-00008",
            state_code="NV"
        )
        db_session.add(patient)
        db_session.commit()
        
        treatment = Treatment(
            treatment_code="NB-TRT-00004",
            patient_id=patient.id,
            primary_condition="chronic_pain",
            primary_cannabinoid_name="CBD"
        )
        db_session.add(treatment)
        db_session.commit()
        
        response = TreatmentResponse(
            patient_id=patient.id,
            treatment_id=treatment.id,
            assessment_date=datetime.now(),
            assessment_type="follow_up",
            days_on_treatment=30,
            primary_symptom_baseline=8.0,
            primary_symptom_current=4.0,
            primary_symptom_change=-4.0,
            symptom_relief_score=7.5,
            satisfaction_score=8.0,
            would_recommend=True
        )
        
        db_session.add(response)
        db_session.commit()
        
        retrieved = db_session.query(TreatmentResponse).filter_by(
            treatment_id=treatment.id
        ).first()
        
        assert retrieved is not None
        assert retrieved.primary_symptom_change == -4.0
        assert retrieved.would_recommend is True
    
    def test_efficacy_score_calculation(self, db_session):
        """Test efficacy score calculation from response."""
        patient = Patient(
            patient_hash="eff123",
            profile_code="NB-00009",
            state_code="NV"
        )
        db_session.add(patient)
        db_session.commit()
        
        treatment = Treatment(
            treatment_code="NB-TRT-00005",
            patient_id=patient.id,
            primary_condition="anxiety",
            primary_cannabinoid_name="CBD"
        )
        db_session.add(treatment)
        db_session.commit()
        
        response = TreatmentResponse(
            patient_id=patient.id,
            treatment_id=treatment.id,
            assessment_date=datetime.now(),
            primary_symptom_change=-5.0,  # 5 point improvement
            symptom_relief_score=8.0,
            functional_improvement=7.0,
            quality_of_life_change=3.0  # +3 on -5 to +5 scale
        )
        
        db_session.add(response)
        db_session.commit()
        
        efficacy = response.calculate_efficacy_score()
        
        # Should be a good score given the positive outcomes
        assert efficacy is not None
        assert efficacy > 5.0  # Better than neutral


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
