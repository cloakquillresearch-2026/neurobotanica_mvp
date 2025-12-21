"""
Week 9 Tests - ChemPath Chemical Characterization Engine

Tests for:
- ChemPath Analyzer (structure validation, descriptor computation)
- COA Quality Control validation
- Report Generator (Markdown, HTML, JSON)
- API endpoints
"""
import pytest
from unittest.mock import patch, MagicMock
from fastapi.testclient import TestClient

from backend.main import app

client = TestClient(app)


# ============================================================================
# ChemPath Analyzer Tests
# ============================================================================

class TestChemPathAnalyzer:
    """Test ChemPath Analyzer functionality."""
    
    def test_analyzer_initialization(self):
        """Test analyzer initializes correctly."""
        from backend.services.chempath.analyzer import ChemPathAnalyzer
        
        analyzer = ChemPathAnalyzer()
        
        assert hasattr(analyzer, 'QC_CODES')
        assert hasattr(analyzer, 'KNOWN_CANNABINOIDS')
        assert "STRUCTURE_INVALID" in analyzer.QC_CODES
        assert "COA_TOTAL_EXCEEDS_100" in analyzer.QC_CODES
    
    def test_qc_severity_enum(self):
        """Test QC severity enum values."""
        from backend.services.chempath.analyzer import QCSeverity
        
        assert QCSeverity.ERROR.value == "error"
        assert QCSeverity.WARNING.value == "warning"
        assert QCSeverity.INFO.value == "info"
    
    def test_qc_flag_dataclass(self):
        """Test QCFlag dataclass."""
        from backend.services.chempath.analyzer import QCFlag, QCSeverity
        
        flag = QCFlag(
            code="TEST_CODE",
            severity=QCSeverity.WARNING,
            message="Test message",
            field="test_field"
        )
        
        flag_dict = flag.to_dict()
        assert flag_dict["code"] == "TEST_CODE"
        assert flag_dict["severity"] == "warning"
        assert flag_dict["message"] == "Test message"
        assert flag_dict["field"] == "test_field"
    
    def test_compound_input_dataclass(self):
        """Test CompoundInput dataclass."""
        from backend.services.chempath.analyzer import CompoundInput
        
        compound = CompoundInput(
            name="Test Compound",
            smiles="CCO",
            source="test"
        )
        
        compound_dict = compound.to_dict()
        assert compound_dict["name"] == "Test Compound"
        assert compound_dict["smiles"] == "CCO"
        assert compound_dict["source"] == "test"
    
    def test_coa_input_dataclass(self):
        """Test COAInput dataclass."""
        from backend.services.chempath.analyzer import COAInput
        
        coa = COAInput(
            lab_name="Test Lab",
            sample_id="SAMPLE001",
            units="%",
            cannabinoids=[{"name": "CBD", "value": 15.0}]
        )
        
        coa_dict = coa.to_dict()
        assert coa_dict["lab_name"] == "Test Lab"
        assert coa_dict["units"] == "%"
        assert len(coa_dict["cannabinoids"]) == 1
    
    def test_analyze_valid_smiles(self):
        """Test analysis with valid SMILES."""
        from backend.services.chempath.analyzer import (
            ChemPathAnalyzer, CompoundInput, ChemPathRequest
        )
        
        analyzer = ChemPathAnalyzer()
        request = ChemPathRequest(
            compound_input=CompoundInput(
                name="Ethanol",
                smiles="CCO"
            )
        )
        
        response = analyzer.analyze(request)
        
        assert response.status == "succeeded"
        assert response.compound_name == "Ethanol"
        assert response.normalized_structure.is_valid is True
        assert response.data_completeness_score > 0
    
    def test_analyze_invalid_smiles(self):
        """Test analysis with invalid SMILES."""
        from backend.services.chempath.analyzer import (
            ChemPathAnalyzer, CompoundInput, ChemPathRequest
        )
        
        analyzer = ChemPathAnalyzer()
        request = ChemPathRequest(
            compound_input=CompoundInput(
                name="Invalid",
                smiles="XXXINVALID((("
            )
        )
        
        response = analyzer.analyze(request)
        
        # Should have structure error
        error_codes = [f.code for f in response.structure_qc_flags]
        assert "STRUCTURE_INVALID" in error_codes
    
    def test_analyze_cbd_smiles(self):
        """Test analysis with CBD SMILES."""
        from backend.services.chempath.analyzer import (
            ChemPathAnalyzer, CompoundInput, ChemPathRequest
        )
        
        analyzer = ChemPathAnalyzer()
        cbd_smiles = "CCCCCC1=CC(=C(C(=C1)O)C2C=C(CCC2C(=C)C)C)O"
        
        request = ChemPathRequest(
            compound_input=CompoundInput(
                name="CBD",
                smiles=cbd_smiles
            )
        )
        
        response = analyzer.analyze(request)
        
        assert response.status == "succeeded"
        assert response.normalized_structure.is_valid is True
        
        # Check descriptors computed
        desc = response.computed_descriptors.descriptors_2d
        assert "molecular_weight" in desc
        assert "logp" in desc
    
    def test_analyze_with_coa(self):
        """Test analysis with COA data."""
        from backend.services.chempath.analyzer import (
            ChemPathAnalyzer, CompoundInput, COAInput, ChemPathRequest
        )
        
        analyzer = ChemPathAnalyzer()
        request = ChemPathRequest(
            compound_input=CompoundInput(
                name="Test Compound",
                smiles="CCO"
            ),
            coa_input=COAInput(
                lab_name="Test Lab",
                sample_id="SAMPLE001",
                units="%",
                cannabinoids=[
                    {"name": "CBD", "value": 15.0},
                    {"name": "THC", "value": 0.3}
                ],
                total_cbd=15.0,
                total_thc=0.3
            )
        )
        
        response = analyzer.analyze(request)
        
        assert response.status == "succeeded"
        # Completeness should be higher with COA
        assert response.data_completeness_score > 40
    
    def test_coa_validation_exceeds_100(self):
        """Test COA validation catches >100% total."""
        from backend.services.chempath.analyzer import (
            ChemPathAnalyzer, CompoundInput, COAInput, ChemPathRequest
        )
        
        analyzer = ChemPathAnalyzer()
        request = ChemPathRequest(
            compound_input=CompoundInput(
                name="Test",
                smiles="CCO"
            ),
            coa_input=COAInput(
                lab_name="Test Lab",
                sample_id="SAMPLE001",
                units="%",
                cannabinoids=[
                    {"name": "CBD", "value": 80.0},
                    {"name": "THC", "value": 30.0}
                ]
            )
        )
        
        response = analyzer.analyze(request)
        
        # Should have error flag
        coa_codes = [f.code for f in response.coa_qc_flags]
        assert "COA_TOTAL_EXCEEDS_100" in coa_codes
    
    def test_coa_validation_negative_value(self):
        """Test COA validation catches negative values."""
        from backend.services.chempath.analyzer import (
            ChemPathAnalyzer, CompoundInput, COAInput, ChemPathRequest
        )
        
        analyzer = ChemPathAnalyzer()
        request = ChemPathRequest(
            compound_input=CompoundInput(
                name="Test",
                smiles="CCO"
            ),
            coa_input=COAInput(
                lab_name="Test Lab",
                sample_id="SAMPLE001",
                units="%",
                cannabinoids=[
                    {"name": "CBD", "value": -5.0}
                ]
            )
        )
        
        response = analyzer.analyze(request)
        
        coa_codes = [f.code for f in response.coa_qc_flags]
        assert "COA_NEGATIVE_VALUE" in coa_codes
    
    def test_coa_validation_missing_lod(self):
        """Test COA validation flags missing LOD."""
        from backend.services.chempath.analyzer import (
            ChemPathAnalyzer, CompoundInput, COAInput, ChemPathRequest
        )
        
        analyzer = ChemPathAnalyzer()
        request = ChemPathRequest(
            compound_input=CompoundInput(
                name="Test",
                smiles="CCO"
            ),
            coa_input=COAInput(
                lab_name="Test Lab",
                sample_id="SAMPLE001",
                units="%",
                cannabinoids=[{"name": "CBD", "value": 15.0}]
                # No lod_loq
            )
        )
        
        response = analyzer.analyze(request)
        
        coa_codes = [f.code for f in response.coa_qc_flags]
        assert "COA_MISSING_LOD" in coa_codes
    
    def test_validate_smiles_method(self):
        """Test quick SMILES validation method."""
        from backend.services.chempath.analyzer import ChemPathAnalyzer
        
        analyzer = ChemPathAnalyzer()
        
        # Valid SMILES
        result = analyzer.validate_smiles("CCO")
        assert result["is_valid"] is True
        
        # Invalid SMILES
        result = analyzer.validate_smiles("INVALID(((")
        assert result["is_valid"] is False
    
    def test_get_qc_codes(self):
        """Test getting QC codes."""
        from backend.services.chempath.analyzer import ChemPathAnalyzer
        
        analyzer = ChemPathAnalyzer()
        codes = analyzer.get_qc_codes()
        
        assert "STRUCTURE_INVALID" in codes
        assert "COA_TOTAL_EXCEEDS_100" in codes
        assert len(codes) > 10
    
    def test_get_known_cannabinoids(self):
        """Test getting known cannabinoid list."""
        from backend.services.chempath.analyzer import ChemPathAnalyzer
        
        analyzer = ChemPathAnalyzer()
        cannabinoids = analyzer.get_known_cannabinoids()
        
        assert "CBD" in cannabinoids
        assert "THC" in cannabinoids
    
    def test_response_to_dict(self):
        """Test ChemPathResponse serialization."""
        from backend.services.chempath.analyzer import (
            ChemPathAnalyzer, CompoundInput, ChemPathRequest
        )
        
        analyzer = ChemPathAnalyzer()
        request = ChemPathRequest(
            compound_input=CompoundInput(name="Test", smiles="CCO")
        )
        
        response = analyzer.analyze(request)
        response_dict = response.to_dict()
        
        assert "chempath_job_id" in response_dict
        assert "compound_name" in response_dict
        assert "normalized_structure" in response_dict
        assert "computed_descriptors" in response_dict
        assert "data_completeness_score" in response_dict


# ============================================================================
# Report Generator Tests
# ============================================================================

class TestChemPathReportGenerator:
    """Test ChemPath Report Generator."""
    
    def test_report_format_enum(self):
        """Test report format enum."""
        from backend.services.chempath.report_generator import ReportFormat
        
        assert ReportFormat.MARKDOWN.value == "markdown"
        assert ReportFormat.HTML.value == "html"
        assert ReportFormat.JSON.value == "json"
    
    def test_generate_markdown_report(self):
        """Test Markdown report generation."""
        from backend.services.chempath.analyzer import (
            ChemPathAnalyzer, CompoundInput, ChemPathRequest
        )
        from backend.services.chempath.report_generator import (
            ChemPathReportGenerator, ReportFormat
        )
        
        analyzer = ChemPathAnalyzer()
        generator = ChemPathReportGenerator()
        
        request = ChemPathRequest(
            compound_input=CompoundInput(name="Test Compound", smiles="CCO")
        )
        response = analyzer.analyze(request)
        
        report = generator.generate_report(response, format=ReportFormat.MARKDOWN)
        
        assert "# Molecular Characterization Report" in report
        assert "Test Compound" in report
        assert "## 1. Identity & Structure" in report
        assert "## 2. Physicochemical Descriptors" in report
    
    def test_generate_html_report(self):
        """Test HTML report generation."""
        from backend.services.chempath.analyzer import (
            ChemPathAnalyzer, CompoundInput, ChemPathRequest
        )
        from backend.services.chempath.report_generator import (
            ChemPathReportGenerator, ReportFormat
        )
        
        analyzer = ChemPathAnalyzer()
        generator = ChemPathReportGenerator()
        
        request = ChemPathRequest(
            compound_input=CompoundInput(name="Test", smiles="CCO")
        )
        response = analyzer.analyze(request)
        
        report = generator.generate_report(response, format=ReportFormat.HTML)
        
        assert "<!DOCTYPE html>" in report
        assert "<html>" in report
        assert "ChemPath Report" in report
    
    def test_generate_json_report(self):
        """Test JSON report generation."""
        import json
        from backend.services.chempath.analyzer import (
            ChemPathAnalyzer, CompoundInput, ChemPathRequest
        )
        from backend.services.chempath.report_generator import (
            ChemPathReportGenerator, ReportFormat
        )
        
        analyzer = ChemPathAnalyzer()
        generator = ChemPathReportGenerator()
        
        request = ChemPathRequest(
            compound_input=CompoundInput(name="Test", smiles="CCO")
        )
        response = analyzer.analyze(request)
        
        report = generator.generate_report(response, format=ReportFormat.JSON)
        
        # Should be valid JSON
        data = json.loads(report)
        assert "report_type" in data
        assert data["report_type"] == "molecular_characterization"
    
    def test_generate_summary(self):
        """Test summary generation."""
        from backend.services.chempath.analyzer import (
            ChemPathAnalyzer, CompoundInput, ChemPathRequest
        )
        from backend.services.chempath.report_generator import ChemPathReportGenerator
        
        analyzer = ChemPathAnalyzer()
        generator = ChemPathReportGenerator()
        
        request = ChemPathRequest(
            compound_input=CompoundInput(name="Test", smiles="CCO")
        )
        response = analyzer.analyze(request)
        
        summary = generator.generate_summary(response)
        
        assert "compound_name" in summary
        assert "job_id" in summary
        assert "structure_valid" in summary
        assert "completeness_score" in summary
        assert "ready_for_toxpath" in summary
    
    def test_report_includes_qc_flags(self):
        """Test report includes QC flag section."""
        from backend.services.chempath.analyzer import (
            ChemPathAnalyzer, CompoundInput, ChemPathRequest
        )
        from backend.services.chempath.report_generator import (
            ChemPathReportGenerator, ReportFormat
        )
        
        analyzer = ChemPathAnalyzer()
        generator = ChemPathReportGenerator()
        
        request = ChemPathRequest(
            compound_input=CompoundInput(name="Test", smiles="CCO")
        )
        response = analyzer.analyze(request)
        
        report = generator.generate_report(response, format=ReportFormat.MARKDOWN)
        
        assert "## 4. Structure QC Flags" in report or "## 5. Data Completeness" in report
    
    def test_report_completeness_visualization(self):
        """Test completeness score visualization in report."""
        from backend.services.chempath.analyzer import (
            ChemPathAnalyzer, CompoundInput, ChemPathRequest
        )
        from backend.services.chempath.report_generator import (
            ChemPathReportGenerator, ReportFormat
        )
        
        analyzer = ChemPathAnalyzer()
        generator = ChemPathReportGenerator()
        
        request = ChemPathRequest(
            compound_input=CompoundInput(name="Test", smiles="CCO")
        )
        response = analyzer.analyze(request)
        
        report = generator.generate_report(response, format=ReportFormat.MARKDOWN)
        
        # Should have progress bar
        assert "█" in report or "░" in report


# ============================================================================
# ChemPath API Endpoint Tests
# ============================================================================

class TestChemPathAPIEndpoints:
    """Test ChemPath API endpoints."""
    
    def test_analyze_endpoint(self):
        """Test /analyze endpoint."""
        response = client.post(
            "/api/v1/chempath/analyze",
            json={
                "compound": {
                    "name": "Ethanol",
                    "smiles": "CCO"
                }
            }
        )
        
        assert response.status_code == 200
        data = response.json()
        assert data["status"] == "success"
        assert "analysis" in data
        assert data["analysis"]["compound_name"] == "Ethanol"
    
    def test_analyze_with_coa(self):
        """Test /analyze with COA data."""
        response = client.post(
            "/api/v1/chempath/analyze",
            json={
                "compound": {
                    "name": "CBD Extract",
                    "smiles": "CCO"
                },
                "coa": {
                    "lab_name": "Test Lab",
                    "sample_id": "SAMPLE001",
                    "units": "%",
                    "cannabinoids": [
                        {"name": "CBD", "value": 15.0}
                    ]
                }
            }
        )
        
        assert response.status_code == 200
        data = response.json()
        assert data["status"] == "success"
    
    def test_validate_smiles_endpoint(self):
        """Test /validate-smiles endpoint."""
        response = client.post(
            "/api/v1/chempath/validate-smiles",
            json={"smiles": "CCO"}
        )
        
        assert response.status_code == 200
        data = response.json()
        assert data["status"] == "success"
        assert data["validation"]["is_valid"] is True
    
    def test_validate_invalid_smiles(self):
        """Test validation of invalid SMILES."""
        response = client.post(
            "/api/v1/chempath/validate-smiles",
            json={"smiles": "INVALID((("}
        )
        
        assert response.status_code == 200
        data = response.json()
        assert data["validation"]["is_valid"] is False
    
    def test_report_endpoint_markdown(self):
        """Test /report endpoint with markdown format."""
        response = client.post(
            "/api/v1/chempath/report?format=markdown",
            json={
                "compound": {
                    "name": "Test",
                    "smiles": "CCO"
                }
            }
        )
        
        assert response.status_code == 200
        data = response.json()
        assert data["status"] == "success"
        assert data["format"] == "markdown"
        assert "# Molecular Characterization Report" in data["report"]
    
    def test_report_endpoint_json(self):
        """Test /report endpoint with JSON format."""
        response = client.post(
            "/api/v1/chempath/report?format=json",
            json={
                "compound": {
                    "name": "Test",
                    "smiles": "CCO"
                }
            }
        )
        
        assert response.status_code == 200
        data = response.json()
        assert data["format"] == "json"
    
    def test_summary_endpoint(self):
        """Test /summary endpoint."""
        response = client.post(
            "/api/v1/chempath/summary",
            json={
                "compound": {
                    "name": "Test",
                    "smiles": "CCO"
                }
            }
        )
        
        assert response.status_code == 200
        data = response.json()
        assert data["status"] == "success"
        assert "summary" in data
        assert "compound_name" in data["summary"]
    
    def test_qc_codes_endpoint(self):
        """Test /qc-codes endpoint."""
        response = client.get("/api/v1/chempath/qc-codes")
        
        assert response.status_code == 200
        data = response.json()
        assert data["status"] == "success"
        assert "qc_codes" in data
        assert "STRUCTURE_INVALID" in data["qc_codes"]
    
    def test_known_cannabinoids_endpoint(self):
        """Test /known-cannabinoids endpoint."""
        response = client.get("/api/v1/chempath/known-cannabinoids")
        
        assert response.status_code == 200
        data = response.json()
        assert data["status"] == "success"
        assert "cannabinoids" in data
        assert "CBD" in data["cannabinoids"]
    
    def test_report_formats_endpoint(self):
        """Test /report-formats endpoint."""
        response = client.get("/api/v1/chempath/report-formats")
        
        assert response.status_code == 200
        data = response.json()
        assert data["status"] == "success"
        assert len(data["formats"]) == 3
    
    def test_descriptor_list_endpoint(self):
        """Test /descriptor-list endpoint."""
        response = client.get("/api/v1/chempath/descriptor-list")
        
        assert response.status_code == 200
        data = response.json()
        assert data["status"] == "success"
        assert "descriptors_2d" in data
        assert "descriptors_3d" in data
        assert "molecular_weight" in data["descriptors_2d"]
    
    def test_health_endpoint(self):
        """Test /health endpoint."""
        response = client.get("/api/v1/chempath/health")
        
        assert response.status_code == 200
        data = response.json()
        assert data["status"] == "healthy"
        assert data["service"] == "chempath"


# ============================================================================
# Integration Tests
# ============================================================================

class TestWeek9Integration:
    """Integration tests for Week 9 ChemPath features."""
    
    def test_full_analysis_pipeline(self):
        """Test full ChemPath analysis pipeline."""
        from backend.services.chempath import (
            ChemPathAnalyzer,
            ChemPathReportGenerator,
            CompoundInput,
            COAInput,
            ChemPathRequest,
            ReportFormat
        )
        
        # Initialize
        analyzer = ChemPathAnalyzer()
        generator = ChemPathReportGenerator()
        
        # Create request with COA
        request = ChemPathRequest(
            compound_input=CompoundInput(
                name="Test CBD Product",
                smiles="CCCCCC1=CC(=C(C(=C1)O)C2C=C(CCC2C(=C)C)C)O"
            ),
            coa_input=COAInput(
                lab_name="Nevada Testing Lab",
                sample_id="NV-2025-001",
                units="%",
                cannabinoids=[
                    {"name": "CBD", "value": 18.5},
                    {"name": "THC", "value": 0.2}
                ],
                terpenes=[
                    {"name": "Myrcene", "value": 1.2},
                    {"name": "Limonene", "value": 0.8}
                ],
                total_cbd=18.5,
                total_thc=0.2,
                lod_loq={"CBD": 0.01, "THC": 0.01}
            )
        )
        
        # Analyze
        response = analyzer.analyze(request)
        
        # Generate report
        report = generator.generate_report(response, format=ReportFormat.MARKDOWN)
        
        # Generate summary
        summary = generator.generate_summary(response)
        
        # Assertions
        assert response.status == "succeeded"
        assert response.data_completeness_score > 60
        assert "Test CBD Product" in report
        assert summary["structure_valid"] is True
    
    def test_chempath_module_imports(self):
        """Test all ChemPath exports are accessible."""
        from backend.services.chempath import (
            ChemPathAnalyzer,
            CompoundInput,
            COAInput,
            ChemPathRequest,
            ChemPathResponse,
            QCFlag,
            QCSeverity,
            NormalizedStructure,
            DescriptorSet,
            ChemPathReportGenerator,
            ReportFormat
        )
        
        # All should be importable
        assert ChemPathAnalyzer is not None
        assert CompoundInput is not None
        assert COAInput is not None
        assert ChemPathRequest is not None
        assert ChemPathResponse is not None
        assert QCFlag is not None
        assert QCSeverity is not None
        assert NormalizedStructure is not None
        assert DescriptorSet is not None
        assert ChemPathReportGenerator is not None
        assert ReportFormat is not None


# ============================================================================
# Validation Tests
# ============================================================================

class TestWeek9Validation:
    """Validation tests for Week 9 requirements."""
    
    def test_descriptor_computation_accuracy(self):
        """Test that computed descriptors are reasonable."""
        from backend.services.chempath.analyzer import (
            ChemPathAnalyzer, CompoundInput, ChemPathRequest
        )
        
        analyzer = ChemPathAnalyzer()
        
        # Test with ethanol (known properties)
        request = ChemPathRequest(
            compound_input=CompoundInput(name="Ethanol", smiles="CCO")
        )
        response = analyzer.analyze(request)
        
        desc = response.computed_descriptors.descriptors_2d
        
        # Ethanol MW should be ~46
        if "molecular_weight" in desc:
            assert 44 < desc["molecular_weight"] < 48
        
        # Ethanol should be hydrophilic (negative or low logP)
        if "logp" in desc:
            assert desc["logp"] < 1
    
    def test_coa_contaminant_detection(self):
        """Test COA contaminant detection."""
        from backend.services.chempath.analyzer import (
            ChemPathAnalyzer, CompoundInput, COAInput, ChemPathRequest
        )
        
        analyzer = ChemPathAnalyzer()
        
        # COA with contaminant failure
        request = ChemPathRequest(
            compound_input=CompoundInput(name="Test", smiles="CCO"),
            coa_input=COAInput(
                lab_name="Test Lab",
                sample_id="SAMPLE001",
                units="%",
                cannabinoids=[{"name": "CBD", "value": 15.0}],
                heavy_metals=[
                    {"name": "Lead", "value": 1.5, "result": "fail"}
                ]
            )
        )
        
        response = analyzer.analyze(request)
        
        coa_codes = [f.code for f in response.coa_qc_flags]
        assert "COA_HEAVY_METAL_FAIL" in coa_codes
    
    def test_data_completeness_scoring(self):
        """Test data completeness scoring logic."""
        from backend.services.chempath.analyzer import (
            ChemPathAnalyzer, CompoundInput, COAInput, ChemPathRequest
        )
        
        analyzer = ChemPathAnalyzer()
        
        # Minimal request
        request_minimal = ChemPathRequest(
            compound_input=CompoundInput(name="Test", smiles="CCO")
        )
        response_minimal = analyzer.analyze(request_minimal)
        
        # Full request
        request_full = ChemPathRequest(
            compound_input=CompoundInput(name="Test", smiles="CCO"),
            coa_input=COAInput(
                lab_name="Test Lab",
                sample_id="SAMPLE001",
                units="%",
                cannabinoids=[{"name": "CBD", "value": 15.0}],
                terpenes=[{"name": "Myrcene", "value": 1.0}],
                total_cbd=15.0,
                heavy_metals=[{"name": "Lead", "result": "pass"}],
                pesticides=[{"name": "Pesticide1", "result": "pass"}],
                solvents=[{"name": "Ethanol", "result": "pass"}],
                lod_loq={"CBD": 0.01}
            )
        )
        response_full = analyzer.analyze(request_full)
        
        # Full should have higher completeness
        assert response_full.data_completeness_score > response_minimal.data_completeness_score
    
    def test_salt_normalization(self):
        """Test salt form detection and normalization."""
        from backend.services.chempath.analyzer import (
            ChemPathAnalyzer, CompoundInput, ChemPathRequest
        )
        
        analyzer = ChemPathAnalyzer()
        
        # Multi-fragment SMILES (more complex salt form)
        request = ChemPathRequest(
            compound_input=CompoundInput(
                name="Salt Form",
                smiles="CC(=O)O.[Na+]"  # Acetic acid with sodium - common salt
            )
        )
        
        response = analyzer.analyze(request)
        
        # Should have processed successfully or flagged the fragment
        # Either SALT_NORMALIZED or FRAGMENT_DETECTED should be present
        codes = [f.code for f in response.structure_qc_flags]
        has_salt_flag = any("SALT" in c or "FRAGMENT" in c for c in codes)
        
        # If RDKit handled it, structure should be valid with notes
        if response.normalized_structure.is_valid:
            assert len(response.normalized_structure.normalization_notes) >= 1
        else:
            # Or it flagged as fragment
            assert has_salt_flag
