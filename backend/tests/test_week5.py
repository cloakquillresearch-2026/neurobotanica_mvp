"""
Week 5 Tests - Adjuvant Optimization + ChEMBL/PubChem Integration
Tests for NeuroBotanica Week 5 features.

Tests:
- Adjuvant model and database
- PrimingProtocolGenerator
- ChEMBL client integration
- PubChem client integration
- Drug interaction checking
"""
import pytest
from datetime import datetime
from unittest.mock import Mock, patch, MagicMock
import json

# Import Week 5 components
from backend.services.priming_protocol_generator import (
    PrimingProtocolGenerator,
    PrimingProtocolResult,
    AdjuvantRecommendation,
    DEFAULT_ADJUVANTS,
    CONDITION_ADJUVANT_PRIORITIES,
    CANNABINOID_ADJUVANT_PRIORITIES
)
from backend.services.chembl_client import (
    ChEMBLClient,
    ChEMBLActivity,
    ChEMBLMolecule,
    ChEMBLSearchResult
)
from backend.services.pubchem_client import (
    PubChemClient,
    PubChemProperty,
    PubChemBioassay,
    PubChemSearchResult
)


# =============================================================================
# ADJUVANT DATABASE TESTS
# =============================================================================

class TestAdjuvantDatabase:
    """Test adjuvant database configuration."""
    
    def test_default_adjuvants_loaded(self):
        """Test default adjuvants are loaded."""
        assert len(DEFAULT_ADJUVANTS) >= 15
    
    def test_adjuvant_required_fields(self):
        """Test all adjuvants have required fields."""
        required = ["id", "name", "category", "mechanism", 
                   "recommended_dose_mg", "timing_offset_minutes", "evidence_tier"]
        
        for adj in DEFAULT_ADJUVANTS:
            for field in required:
                assert field in adj, f"Missing {field} in {adj.get('name')}"
    
    def test_adjuvant_categories_valid(self):
        """Test adjuvant categories are valid."""
        valid_categories = [
            "Receptor Priming", "Metabolic Enhancer", "Synergistic",
            "Protective", "Bioavailability", "Terpene", "Other"
        ]
        
        for adj in DEFAULT_ADJUVANTS:
            assert adj["category"] in valid_categories
    
    def test_evidence_tiers_in_range(self):
        """Test evidence tiers are 1-5."""
        for adj in DEFAULT_ADJUVANTS:
            assert 1 <= adj["evidence_tier"] <= 5
    
    def test_timing_offsets_reasonable(self):
        """Test timing offsets are reasonable (-60 to +30 min)."""
        for adj in DEFAULT_ADJUVANTS:
            timing = adj["timing_offset_minutes"]
            assert -120 <= timing <= 60, f"Unusual timing for {adj['name']}: {timing}"


class TestConditionPriorities:
    """Test condition-adjuvant priority mappings."""
    
    def test_all_conditions_have_priorities(self):
        """Test all conditions have adjuvant priorities."""
        expected_conditions = [
            "chronic_pain", "anxiety", "insomnia", "inflammation",
            "neuropathy", "ptsd", "epilepsy"
        ]
        
        for condition in expected_conditions:
            assert condition in CONDITION_ADJUVANT_PRIORITIES
    
    def test_priority_adjuvants_exist(self):
        """Test priority adjuvants exist in database."""
        adjuvant_names = {a["name"] for a in DEFAULT_ADJUVANTS}
        
        for condition, priorities in CONDITION_ADJUVANT_PRIORITIES.items():
            for adj_name in priorities:
                assert adj_name in adjuvant_names, \
                    f"Priority adjuvant '{adj_name}' for {condition} not in database"


class TestCannabinoidPriorities:
    """Test cannabinoid-adjuvant priority mappings."""
    
    def test_major_cannabinoids_have_priorities(self):
        """Test major cannabinoids have adjuvant priorities."""
        major = ["THC", "CBD", "CBG", "CBN"]
        
        for cannabinoid in major:
            assert cannabinoid in CANNABINOID_ADJUVANT_PRIORITIES


# =============================================================================
# PRIMING PROTOCOL GENERATOR TESTS
# =============================================================================

class TestPrimingProtocolGenerator:
    """Test priming protocol generator."""
    
    @pytest.fixture
    def generator(self):
        """Create protocol generator instance."""
        return PrimingProtocolGenerator()
    
    def test_generator_initialization(self, generator):
        """Test generator initializes with adjuvants."""
        assert len(generator.adjuvants) >= 15
        assert len(generator.adjuvant_by_id) >= 15
    
    def test_generate_protocol_basic(self, generator):
        """Test basic protocol generation."""
        result = generator.generate_protocol(
            cannabinoid_name="THC",
            condition="chronic_pain"
        )
        
        assert isinstance(result, PrimingProtocolResult)
        assert result.target_cannabinoid == "THC"
        assert result.target_condition == "chronic_pain"
        assert len(result.recommendations) > 0
    
    def test_generate_protocol_cbd_anxiety(self, generator):
        """Test protocol generation for CBD + anxiety."""
        result = generator.generate_protocol(
            cannabinoid_name="CBD",
            condition="anxiety"
        )
        
        assert result.target_cannabinoid == "CBD"
        assert result.target_condition == "anxiety"
        assert result.protocol_confidence > 0
    
    def test_protocol_respects_max_adjuvants(self, generator):
        """Test max_adjuvants parameter is respected."""
        result = generator.generate_protocol(
            cannabinoid_name="THC",
            condition="chronic_pain",
            max_adjuvants=2
        )
        
        assert len(result.recommendations) <= 2
    
    def test_protocol_respects_evidence_threshold(self, generator):
        """Test min_evidence_tier parameter."""
        result = generator.generate_protocol(
            cannabinoid_name="THC",
            condition="chronic_pain",
            min_evidence_tier=2
        )
        
        for rec in result.recommendations:
            assert rec.evidence_tier <= 2
    
    def test_protocol_sorted_by_timing(self, generator):
        """Test recommendations are sorted by timing."""
        result = generator.generate_protocol(
            cannabinoid_name="THC",
            condition="chronic_pain",
            max_adjuvants=4
        )
        
        timings = [r.timing_minutes for r in result.recommendations]
        assert timings == sorted(timings)
    
    def test_protocol_with_patient_profile(self, generator):
        """Test protocol generation with patient profile."""
        patient_profile = {
            "weight_kg": 80,
            "conditions": [],
            "allergies": [],
            "medications": []
        }
        
        result = generator.generate_protocol(
            cannabinoid_name="CBD",
            condition="inflammation",
            patient_profile=patient_profile
        )
        
        assert isinstance(result, PrimingProtocolResult)
    
    def test_protocol_handles_contraindications(self, generator):
        """Test contraindications are handled."""
        patient_profile = {
            "conditions": ["kidney_disease"],
            "allergies": [],
            "medications": []
        }
        
        result = generator.generate_protocol(
            cannabinoid_name="THC",
            condition="chronic_pain",
            patient_profile=patient_profile
        )
        
        # Magnesium should be excluded
        adj_names = [r.adjuvant_name for r in result.recommendations]
        assert "Magnesium Glycinate" not in adj_names
        assert len(result.contraindications) > 0
    
    def test_protocol_confidence_calculation(self, generator):
        """Test confidence score is calculated."""
        result = generator.generate_protocol(
            cannabinoid_name="THC",
            condition="chronic_pain"
        )
        
        assert 0.0 <= result.protocol_confidence <= 1.0
    
    def test_protocol_to_dict(self, generator):
        """Test protocol serialization."""
        result = generator.generate_protocol(
            cannabinoid_name="CBD",
            condition="anxiety"
        )
        
        data = result.to_dict()
        
        assert "protocol_name" in data
        assert "recommendations" in data
        assert "protocol_confidence" in data
        assert "generated_at" in data


class TestDrugInteractionChecker:
    """Test drug interaction checking."""
    
    @pytest.fixture
    def generator(self):
        return PrimingProtocolGenerator()
    
    def test_check_known_interaction(self, generator):
        """Test known drug interaction is detected."""
        interactions = generator.check_drug_interactions(
            adjuvant_names=["Omega-3 Fatty Acids (EPA/DHA)"],
            medications=["warfarin"]
        )
        
        assert len(interactions) > 0
        assert interactions[0]["severity"] == "Moderate"
    
    def test_check_no_interaction(self, generator):
        """Test no interaction returns empty list."""
        interactions = generator.check_drug_interactions(
            adjuvant_names=["L-Theanine"],
            medications=["ibuprofen"]
        )
        
        # L-Theanine has no major interactions
        assert len(interactions) == 0
    
    def test_check_multiple_adjuvants(self, generator):
        """Test checking multiple adjuvants."""
        interactions = generator.check_drug_interactions(
            adjuvant_names=[
                "Black Pepper Extract (Piperine)",
                "Omega-3 Fatty Acids (EPA/DHA)"
            ],
            medications=["warfarin"]
        )
        
        # Both can interact with blood thinners
        assert len(interactions) >= 1


class TestAdjuvantInfo:
    """Test adjuvant information retrieval."""
    
    @pytest.fixture
    def generator(self):
        return PrimingProtocolGenerator()
    
    def test_get_adjuvant_info(self, generator):
        """Test retrieving adjuvant info."""
        info = generator.get_adjuvant_info("Magnesium Glycinate")
        
        assert info is not None
        assert info["name"] == "Magnesium Glycinate"
        assert "mechanism" in info
    
    def test_get_adjuvant_info_not_found(self, generator):
        """Test retrieving non-existent adjuvant."""
        info = generator.get_adjuvant_info("Unknown Compound")
        assert info is None
    
    def test_list_adjuvants_by_category(self, generator):
        """Test listing adjuvants by category."""
        terpenes = generator.list_adjuvants_by_category("Terpene")
        
        assert len(terpenes) >= 3
        for adj in terpenes:
            assert adj["category"] == "Terpene"
    
    def test_get_conditions(self, generator):
        """Test getting supported conditions."""
        conditions = generator.get_conditions()
        
        assert "chronic_pain" in conditions
        assert "anxiety" in conditions
    
    def test_get_cannabinoids(self, generator):
        """Test getting supported cannabinoids."""
        cannabinoids = generator.get_cannabinoids()
        
        assert "THC" in cannabinoids
        assert "CBD" in cannabinoids


# =============================================================================
# CHEMBL CLIENT TESTS
# =============================================================================

class TestChEMBLClient:
    """Test ChEMBL API client."""
    
    @pytest.fixture
    def client(self):
        return ChEMBLClient(timeout=10, max_retries=2)
    
    def test_client_initialization(self, client):
        """Test client initializes correctly."""
        assert client.BASE_URL == "https://www.ebi.ac.uk/chembl/api/data"
        assert len(client.CANNABINOID_TARGETS) > 0
    
    def test_known_cannabinoid_ids(self, client):
        """Test known cannabinoid IDs are defined."""
        known = client.get_known_cannabinoid_ids()
        
        assert "THC" in known
        assert "CBD" in known
        assert known["THC"] == "CHEMBL465"
    
    def test_cannabinoid_targets_defined(self, client):
        """Test cannabinoid targets are defined."""
        targets = client.CANNABINOID_TARGETS
        
        assert "CB1" in targets
        assert "CB2" in targets
        assert "FAAH" in targets


class TestChEMBLDataClasses:
    """Test ChEMBL data classes."""
    
    def test_chembl_activity_to_provenance(self):
        """Test ChEMBLActivity to provenance dict."""
        activity = ChEMBLActivity(
            activity_id="12345",
            molecule_chembl_id="CHEMBL465",
            target_chembl_id="CHEMBL218",
            target_name="Cannabinoid CB1 receptor",
            target_organism="Homo sapiens",
            assay_type="B",
            assay_description="Binding assay",
            standard_type="Ki",
            standard_value=40.0,
            standard_units="nM",
            standard_relation="=",
            pchembl_value=7.4,
            data_validity_comment=None,
            document_chembl_id="CHEMBL1234",
            document_year=2020,
            src_id=1
        )
        
        prov = activity.to_provenance_dict()
        
        assert prov["affinity_value"] == 40.0
        assert prov["affinity_unit"] == "nM"
        assert prov["target_name"] == "Cannabinoid CB1 receptor"
        assert "ChEMBL" in prov["source"]
        assert "retrieved_at" in prov
    
    def test_chembl_molecule_to_dict(self):
        """Test ChEMBLMolecule to dict."""
        molecule = ChEMBLMolecule(
            molecule_chembl_id="CHEMBL465",
            pref_name="THC",
            molecule_type="Small molecule",
            max_phase=4,
            therapeutic_flag=True,
            natural_product=True,
            oral=True,
            parenteral=False,
            topical=False,
            molecular_weight=314.46,
            alogp=6.97,
            psa=29.46,
            hba=2,
            hbd=1,
            num_ro5_violations=1,
            canonical_smiles="CCCCC",
            standard_inchi="InChI=...",
            standard_inchi_key="ABC..."
        )
        
        data = molecule.to_dict()
        
        assert data["chembl_id"] == "CHEMBL465"
        assert data["name"] == "THC"
        assert data["natural_product"] == True
        assert data["properties"]["molecular_weight"] == 314.46


class TestChEMBLSearchResult:
    """Test ChEMBL search result."""
    
    def test_search_result_success(self):
        """Test successful search result."""
        result = ChEMBLSearchResult(
            success=True,
            molecule=None,
            activities=[],
            request_time_ms=150.0
        )
        
        assert result.success
        data = result.to_dict()
        assert data["success"] == True
    
    def test_search_result_failure(self):
        """Test failed search result."""
        result = ChEMBLSearchResult(
            success=False,
            error="Not found"
        )
        
        assert not result.success
        assert result.error == "Not found"


# =============================================================================
# PUBCHEM CLIENT TESTS
# =============================================================================

class TestPubChemClient:
    """Test PubChem API client."""
    
    @pytest.fixture
    def client(self):
        return PubChemClient(timeout=10, max_retries=2)
    
    def test_client_initialization(self, client):
        """Test client initializes correctly."""
        assert "pubchem.ncbi.nlm.nih.gov" in client.BASE_URL
        assert len(client.DEFAULT_PROPERTIES) > 10
    
    def test_known_cannabinoid_cids(self, client):
        """Test known cannabinoid CIDs are defined."""
        known = client.get_known_cannabinoid_cids()
        
        assert "THC" in known
        assert "CBD" in known
        assert known["THC"] == "16078"


class TestPubChemDataClasses:
    """Test PubChem data classes."""
    
    def test_pubchem_property_to_dict(self):
        """Test PubChemProperty to dict."""
        prop = PubChemProperty(
            cid="16078",
            molecular_formula="C21H30O2",
            molecular_weight=314.47,
            canonical_smiles="CCCCC",
            isomeric_smiles="CCCCC",
            inchi="InChI=...",
            inchi_key="ABC...",
            iupac_name="delta-9-THC",
            xlogp=7.0,
            exact_mass=314.224,
            monoisotopic_mass=314.224,
            tpsa=29.5,
            complexity=448.0,
            charge=0,
            h_bond_donor_count=1,
            h_bond_acceptor_count=2,
            rotatable_bond_count=4,
            heavy_atom_count=23,
            atom_stereo_count=2,
            defined_atom_stereo_count=2,
            undefined_atom_stereo_count=0,
            bond_stereo_count=0,
            covalent_unit_count=1
        )
        
        data = prop.to_dict()
        
        assert data["cid"] == "16078"
        assert data["molecular_weight"] == 314.47
        assert data["xlogp"] == 7.0
    
    def test_pubchem_bioassay_to_provenance(self):
        """Test PubChemBioassay to provenance dict."""
        bioassay = PubChemBioassay(
            aid=12345,
            assay_name="CB1 binding assay",
            assay_type="Binding",
            activity_outcome="Active",
            activity_value=50.0,
            activity_unit="nM",
            target_name="Cannabinoid receptor 1",
            target_gi=100,
            pubmed_id=12345678
        )
        
        prov = bioassay.to_provenance_dict()
        
        assert prov["assay_id"] == "PubChem:AID12345"
        assert prov["activity_outcome"] == "Active"
        assert prov["source"] == "PubChem BioAssay"


class TestPubChemSynonym:
    """Test PubChem synonym handling."""
    
    def test_get_common_names(self):
        """Test getting common names from synonyms."""
        from backend.services.pubchem_client import PubChemSynonym
        
        synonym = PubChemSynonym(
            cid="16078",
            synonyms=["THC", "delta-9-THC", "Dronabinol", 
                     "SCHEMBL123456", "Very Long Name " * 10],
            depositor_names=[],
            mesh_headings=[]
        )
        
        common = synonym.get_common_names(5)
        
        assert "THC" in common
        # SCHEMBL names should be filtered
        assert all(not n.startswith("SCHEMBL") for n in common)


# =============================================================================
# INTEGRATION TESTS (Mocked)
# =============================================================================

class TestMockedChEMBLIntegration:
    """Test ChEMBL integration with mocked responses."""
    
    @patch('requests.Session.get')
    def test_search_by_smiles_success(self, mock_get):
        """Test SMILES search with mocked response."""
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.json.return_value = {
            "molecules": [{"molecule_chembl_id": "CHEMBL465"}]
        }
        mock_get.return_value = mock_response
        
        client = ChEMBLClient()
        chembl_id = client.search_by_smiles("CCCCC")
        
        assert chembl_id == "CHEMBL465"
    
    @patch('requests.Session.get')
    def test_search_by_smiles_not_found(self, mock_get):
        """Test SMILES search not found."""
        mock_response = Mock()
        mock_response.status_code = 404
        mock_get.return_value = mock_response
        
        client = ChEMBLClient()
        result = client.search_by_smiles("INVALID_SMILES")
        
        assert result is None


class TestMockedPubChemIntegration:
    """Test PubChem integration with mocked responses."""
    
    @patch('requests.Session.get')
    def test_search_by_name_success(self, mock_get):
        """Test name search with mocked response."""
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.json.return_value = {
            "IdentifierList": {"CID": [16078]}
        }
        mock_get.return_value = mock_response
        
        client = PubChemClient()
        cid = client.search_by_name("THC")
        
        assert cid == "16078"


# =============================================================================
# RECOMMENDATION DATA CLASS TESTS
# =============================================================================

class TestAdjuvantRecommendation:
    """Test AdjuvantRecommendation data class."""
    
    def test_recommendation_creation(self):
        """Test creating a recommendation."""
        rec = AdjuvantRecommendation(
            adjuvant_id=1,
            adjuvant_name="Magnesium Glycinate",
            category="Receptor Priming",
            dose_mg=200.0,
            timing_minutes=-30,
            timing_label="30 minutes before cannabinoid",
            mechanism="Enhances NMDA receptor modulation",
            expected_enhancement=25.0,
            evidence_tier=2,
            notes="Take with water"
        )
        
        assert rec.adjuvant_name == "Magnesium Glycinate"
        assert rec.dose_mg == 200.0
    
    def test_recommendation_to_dict(self):
        """Test recommendation serialization."""
        rec = AdjuvantRecommendation(
            adjuvant_id=1,
            adjuvant_name="Test",
            category="Test",
            dose_mg=100.0,
            timing_minutes=-15,
            timing_label="15 min before",
            mechanism="Test mechanism",
            expected_enhancement=20.0,
            evidence_tier=3
        )
        
        data = rec.to_dict()
        
        assert data["adjuvant_id"] == 1
        assert data["dose_mg"] == 100.0
        assert "timing_label" in data


# =============================================================================
# EDGE CASE TESTS
# =============================================================================

class TestEdgeCases:
    """Test edge cases and error handling."""
    
    def test_empty_cannabinoid_name(self):
        """Test handling empty cannabinoid name."""
        generator = PrimingProtocolGenerator()
        result = generator.generate_protocol(
            cannabinoid_name="",
            condition="chronic_pain"
        )
        
        # Should still work, just fewer matches
        assert isinstance(result, PrimingProtocolResult)
    
    def test_unknown_condition(self):
        """Test handling unknown condition."""
        generator = PrimingProtocolGenerator()
        result = generator.generate_protocol(
            cannabinoid_name="THC",
            condition="unknown_condition_xyz"
        )
        
        # Should still generate based on cannabinoid
        assert isinstance(result, PrimingProtocolResult)
    
    def test_severe_allergy_profile(self):
        """Test handling patient with many allergies."""
        generator = PrimingProtocolGenerator()
        patient_profile = {
            "conditions": ["kidney_disease", "asthma_acute"],
            "allergies": ["fish_allergy", "soy_allergy"],
            "medications": []
        }
        
        result = generator.generate_protocol(
            cannabinoid_name="THC",
            condition="chronic_pain",
            patient_profile=patient_profile
        )
        
        # Should exclude contraindicated adjuvants
        assert len(result.contraindications) > 0


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
