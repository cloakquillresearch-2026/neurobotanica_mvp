"""
NeuroBotanica Unit Tests - Conformer Generation
===============================================

Test 3D conformer generation service with ETKDG method

Run tests:
    pytest backend/tests/test_conformers.py -v
    
Note: These tests require RDKit to be installed.
"""
import pytest
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

# Check for RDKit availability
try:
    from rdkit import Chem
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False

# Import conformer generator if RDKit is available
if RDKIT_AVAILABLE:
    from services.conformer_generator import (
        ConformerGenerator,
        BatchConformerGenerator,
        ConformerResult,
        generate_conformers_for_smiles,
        calculate_3d_descriptors
    )


# Skip all tests if RDKit not available
pytestmark = pytest.mark.skipif(
    not RDKIT_AVAILABLE,
    reason="RDKit not installed - conformer tests skipped"
)


class TestConformerGenerator:
    """Tests for ConformerGenerator class."""
    
    @pytest.fixture
    def generator(self):
        """Create conformer generator with default settings."""
        return ConformerGenerator(num_conformers=10, random_seed=42)
    
    def test_generator_initialization(self, generator):
        """Test generator initializes with correct parameters."""
        assert generator.num_conformers == 10
        assert generator.random_seed == 42
        assert generator.force_field == "MMFF94"
    
    def test_simple_molecule(self, generator):
        """Test conformer generation for simple molecule (ethanol)."""
        result = generator.generate_conformers("CCO")
        
        assert result.success is True
        assert result.num_conformers_generated > 0
        assert result.num_conformers_generated <= 10
        assert result.lowest_energy is not None
        assert result.descriptors_3d is not None
    
    def test_cbd_molecule(self, generator):
        """Test conformer generation for CBD."""
        # CBD SMILES
        cbd_smiles = "CCCCCC1=CC(=C(C(=C1)O)C2C=C(CCC2C(=C)C)C)O"
        result = generator.generate_conformers(cbd_smiles)
        
        assert result.success is True
        assert result.num_conformers_generated > 0
        assert result.conformer_metadata is not None
        assert "lowest_energy_kcal_mol" in result.conformer_metadata
    
    def test_thc_molecule(self, generator):
        """Test conformer generation for THC."""
        # THC SMILES
        thc_smiles = "CCCCCC1=CC(=C2C3C=C(CCC3C(OC2=C1)(C)C)C)O"
        result = generator.generate_conformers(thc_smiles)
        
        assert result.success is True
        assert result.num_conformers_generated > 0
    
    def test_invalid_smiles(self, generator):
        """Test handling of invalid SMILES."""
        result = generator.generate_conformers("invalid_smiles_xyz")
        
        assert result.success is False
        assert result.error is not None
        assert "Invalid SMILES" in result.error
    
    def test_empty_smiles(self, generator):
        """Test handling of empty SMILES."""
        result = generator.generate_conformers("")
        
        assert result.success is False
    
    def test_3d_descriptors_calculated(self, generator):
        """Test that 3D descriptors are calculated correctly."""
        result = generator.generate_conformers("CCCC")  # Butane
        
        assert result.success is True
        desc = result.descriptors_3d
        
        # Check key descriptors are present
        assert "pmi1" in desc
        assert "pmi2" in desc
        assert "pmi3" in desc
        assert "asphericity" in desc
        assert "radius_of_gyration" in desc
        
        # Check values are reasonable (non-negative)
        assert desc["pmi1"] >= 0
        assert desc["radius_of_gyration"] > 0
    
    def test_energy_ranking(self, generator):
        """Test conformers are ranked by energy."""
        result = generator.generate_conformers("CCCCCC")  # Hexane
        
        assert result.success is True
        assert result.lowest_energy <= result.highest_energy
        assert result.energy_range >= 0
    
    def test_rmsd_clustering(self, generator):
        """Test RMSD-based clustering."""
        # Use a larger conformer count for meaningful clustering
        large_generator = ConformerGenerator(num_conformers=30, random_seed=42)
        result = large_generator.generate_conformers("CCCCCCCC")  # Octane
        
        assert result.success is True
        assert result.num_clusters_rmsd_2A >= 1
        assert result.num_clusters_rmsd_2A <= result.num_conformers_generated
    
    def test_reproducibility(self, generator):
        """Test that same seed produces same results."""
        gen1 = ConformerGenerator(num_conformers=5, random_seed=42)
        gen2 = ConformerGenerator(num_conformers=5, random_seed=42)
        
        result1 = gen1.generate_conformers("CCC")
        result2 = gen2.generate_conformers("CCC")
        
        assert result1.num_conformers_generated == result2.num_conformers_generated
        assert result1.lowest_energy == result2.lowest_energy
    
    def test_single_conformer_generation(self, generator):
        """Test single conformer generation for quick analysis."""
        desc = generator.generate_single_conformer("CCO")
        
        assert desc is not None
        assert "asphericity" in desc
        assert "radius_of_gyration" in desc


class TestBatchConformerGenerator:
    """Tests for BatchConformerGenerator class."""
    
    def test_batch_processing(self):
        """Test processing multiple compounds."""
        generator = ConformerGenerator(num_conformers=5)
        batch = BatchConformerGenerator(generator)
        
        compounds = [
            {"name": "Ethanol", "smiles": "CCO"},
            {"name": "Propanol", "smiles": "CCCO"},
            {"name": "Invalid", "smiles": "not_valid"}
        ]
        
        result = batch.process_batch(compounds)
        
        assert result["total_compounds"] == 3
        assert result["success_count"] == 2
        assert result["failed_count"] == 1
    
    def test_get_failed_compounds(self):
        """Test retrieving failed compounds."""
        generator = ConformerGenerator(num_conformers=5)
        batch = BatchConformerGenerator(generator)
        
        compounds = [
            {"name": "Good", "smiles": "CCO"},
            {"name": "Bad", "smiles": "invalid"}
        ]
        
        batch.process_batch(compounds)
        failed = batch.get_failed_compounds()
        
        assert len(failed) == 1
        assert failed[0]["name"] == "Bad"


class TestConvenienceFunctions:
    """Tests for convenience functions."""
    
    def test_generate_conformers_for_smiles(self):
        """Test quick conformer generation function."""
        result = generate_conformers_for_smiles("CCO", num_conformers=5)
        
        assert result.success is True
        assert result.num_conformers_generated > 0
    
    def test_calculate_3d_descriptors(self):
        """Test quick 3D descriptor calculation."""
        desc = calculate_3d_descriptors("CCO")
        
        assert desc is not None
        assert "asphericity" in desc


class TestConformerResult:
    """Tests for ConformerResult dataclass."""
    
    def test_result_creation(self):
        """Test creating ConformerResult."""
        result = ConformerResult(
            success=True,
            num_conformers_generated=10,
            num_conformers_requested=10,
            lowest_energy=-50.5,
            highest_energy=-40.2,
            energy_range=10.3,
            num_clusters_rmsd_2A=3
        )
        
        assert result.success is True
        assert result.num_conformers_generated == 10
        assert result.energy_range == 10.3
    
    def test_failed_result(self):
        """Test creating failed ConformerResult."""
        result = ConformerResult(
            success=False,
            error="Generation failed"
        )
        
        assert result.success is False
        assert result.error == "Generation failed"
        assert result.num_conformers_generated == 0


class TestEdgeCases:
    """Tests for edge cases and unusual molecules."""
    
    def test_charged_molecule(self):
        """Test conformer generation for charged molecule."""
        generator = ConformerGenerator(num_conformers=5)
        # Acetate ion
        result = generator.generate_conformers("CC(=O)[O-]")
        # May succeed or fail depending on force field
        assert isinstance(result.success, bool)
    
    def test_aromatic_molecule(self):
        """Test conformer generation for aromatic molecule."""
        generator = ConformerGenerator(num_conformers=5)
        # Benzene
        result = generator.generate_conformers("c1ccccc1")
        
        assert result.success is True
    
    def test_large_flexible_molecule(self):
        """Test conformer generation for large flexible molecule."""
        generator = ConformerGenerator(num_conformers=20)
        # Long alkane chain
        result = generator.generate_conformers("CCCCCCCCCCCCCC")  # Tetradecane
        
        assert result.success is True
        # Should have multiple RMSD clusters due to flexibility
        assert result.num_clusters_rmsd_2A >= 1
    
    def test_ring_system(self):
        """Test conformer generation for ring system."""
        generator = ConformerGenerator(num_conformers=10)
        # Cyclohexane
        result = generator.generate_conformers("C1CCCCC1")
        
        assert result.success is True
        assert result.num_conformers_generated > 0


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
