"""
Simple test to verify NeuroBotanica installation
================================================

Run this test to confirm everything is set up correctly:
    pytest tests/test_installation.py -v
"""

import json
from pathlib import Path
import pytest


def test_python_version():
    """Verify Python version is 3.11+"""
    import sys
    version = sys.version_info
    assert version.major == 3
    assert version.minor >= 11, f"Python 3.11+ required, found {version.major}.{version.minor}"


def test_required_packages():
    """Verify all required packages are installed"""
    try:
        import numpy
        import pandas
        import sklearn
        import matplotlib
        import pytest
        assert True
    except ImportError as e:
        pytest.fail(f"Required package missing: {e}")


def test_rdkit_available():
    """Verify RDKit is installed (most common installation issue)"""
    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors
        
        # Test basic functionality
        mol = Chem.MolFromSmiles("CCO")  # Ethanol
        assert mol is not None
        
        mw = Descriptors.MolWt(mol)
        assert 45 < mw < 47  # Molecular weight should be ~46
        
    except ImportError:
        pytest.fail("RDKit not installed. Try: pip install rdkit OR conda install -c conda-forge rdkit")


def test_project_structure():
    """Verify project directory structure exists"""
    project_root = Path(__file__).parent.parent
    
    required_dirs = [
        project_root / "src",
        project_root / "src/analysis",
        project_root / "src/ml_models",
        project_root / "data",
        project_root / "data/training",
        project_root / "tests"
    ]
    
    for directory in required_dirs:
        assert directory.exists(), f"Missing directory: {directory}"


def test_training_data_exists():
    """Verify training data files are present"""
    project_root = Path(__file__).parent.parent
    
    required_files = [
        project_root / "data/training/neurobotanica_complete_dataset_63compounds.json",
        project_root / "data/training/neurobotanica_dimeric_predictions.json"
    ]
    
    for filepath in required_files:
        assert filepath.exists(), f"Missing data file: {filepath}"


def test_dataset_structure():
    """Verify dataset has correct structure"""
    project_root = Path(__file__).parent.parent
    dataset_path = project_root / "data/training/neurobotanica_complete_dataset_63compounds.json"
    
    with open(dataset_path, 'r') as f:
        data = json.load(f)
    
    # Verify main keys exist
    assert 'compounds' in data, "Dataset missing 'compounds' key"
    assert 'clinical_studies' in data, "Dataset missing 'clinical_studies' key"
    
    # Verify we have cannabinoids
    assert len(data['compounds']) > 0, "No compounds in dataset"
    
    # Verify structure of first cannabinoid
    first_cannabinoid = data['compounds'][0]
    required_fields = ['compound_name', 'smiles', 'rdkit_descriptors']
    for field in required_fields:
        assert field in first_cannabinoid, f"Cannabinoid missing field: {field}"
    
    # Verify RDKit descriptors exist
    assert 'logP' in first_cannabinoid['rdkit_descriptors']
    assert 'tpsa' in first_cannabinoid['rdkit_descriptors']


def test_can_import_project_modules():
    """Verify project modules can be imported"""
    import sys
    from pathlib import Path
    
    project_root = Path(__file__).parent.parent
    sys.path.insert(0, str(project_root))
    
    # These should not raise ImportError
    # (Even if modules are empty, they should exist)
    try:
        from src.analysis import cannabinoid_analyzer
        assert True
    except ImportError:
        pytest.fail("Cannot import project modules. Check __init__.py files exist.")


if __name__ == "__main__":
    # Run tests when executed directly
    pytest.main([__file__, "-v"])
