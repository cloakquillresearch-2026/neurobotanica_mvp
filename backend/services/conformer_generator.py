"""
3D Conformer Generator Service
ETKDG-based conformer generation with MMFF94 optimization

Patent Claims Support:
- 3D molecular descriptor calculation (PMI, asphericity, etc.)
- Conformer ensemble generation for receptor docking preparation
- Energy-ranked conformer clustering

Reference: NeuroBotanica MVP Development Plan - Week 1 Days 4-5
"""
from typing import Dict, List, Tuple, Optional, Any
from dataclasses import dataclass
from datetime import datetime
import logging
import json

# RDKit imports
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, Descriptors3D, rdMolAlign
    from rdkit.Chem.rdMolDescriptors import (
        CalcPBF,  # Plane of Best Fit
    )
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    logging.warning("RDKit not available. Conformer generation will be disabled.")

import numpy as np

logger = logging.getLogger(__name__)


@dataclass
class ConformerResult:
    """Result container for conformer generation."""
    success: bool
    num_conformers_generated: int = 0
    num_conformers_requested: int = 0
    lowest_energy: Optional[float] = None
    highest_energy: Optional[float] = None
    energy_range: Optional[float] = None
    num_clusters_rmsd_2A: int = 0
    descriptors_3d: Optional[Dict] = None
    conformer_metadata: Optional[Dict] = None
    error: Optional[str] = None
    generation_time_seconds: float = 0.0


class ConformerGenerator:
    """Generate and analyze 3D conformers using ETKDG method.
    
    This service implements the conformer generation pipeline from
    NeuroBotanica MVP Development Plan Week 1 Days 4-5.
    
    Features:
    - ETKDG (Experimental-Torsion Knowledge Distance Geometry) conformer generation
    - MMFF94 force field optimization
    - Energy-ranked conformer selection
    - RMSD-based clustering
    - 3D descriptor calculation (PMI, asphericity, eccentricity, etc.)
    
    Usage:
        generator = ConformerGenerator(num_conformers=50)
        result = generator.generate_conformers("CCO")  # Ethanol example
        
        if result.success:
            print(f"Generated {result.num_conformers_generated} conformers")
            print(f"3D descriptors: {result.descriptors_3d}")
    """
    
    def __init__(
        self,
        num_conformers: int = 50,
        max_attempts: int = 100,
        random_seed: int = 42,
        rmsd_threshold: float = 2.0,
        prune_rms_thresh: float = 0.5,
        use_random_coords: bool = False,
        force_field: str = "MMFF94"
    ):
        """Initialize conformer generator.
        
        Args:
            num_conformers: Number of conformers to generate (default: 50)
            max_attempts: Maximum embedding attempts (default: 100)
            random_seed: Random seed for reproducibility (default: 42)
            rmsd_threshold: RMSD threshold for clustering in Angstroms (default: 2.0)
            prune_rms_thresh: RMSD threshold for pruning similar conformers (default: 0.5)
            use_random_coords: Use random coordinates for embedding (default: False)
            force_field: Force field for optimization ("MMFF94" or "UFF")
        """
        if not RDKIT_AVAILABLE:
            raise RuntimeError("RDKit is required for conformer generation. Install with: pip install rdkit")
        
        self.num_conformers = num_conformers
        self.max_attempts = max_attempts
        self.random_seed = random_seed
        self.rmsd_threshold = rmsd_threshold
        self.prune_rms_thresh = prune_rms_thresh
        self.use_random_coords = use_random_coords
        self.force_field = force_field
        
        logger.info(f"ConformerGenerator initialized: {num_conformers} conformers, "
                   f"seed={random_seed}, force_field={force_field}")
    
    def generate_conformers(self, smiles: str) -> ConformerResult:
        """Generate conformer ensemble for a molecule.
        
        Args:
            smiles: SMILES string of the molecule
            
        Returns:
            ConformerResult with conformer data, energies, and 3D descriptors
        """
        import time
        start_time = time.time()
        
        # Validate SMILES
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return ConformerResult(
                success=False,
                error=f"Invalid SMILES string: {smiles}"
            )
        
        # Add hydrogens for proper 3D structure
        mol = Chem.AddHs(mol)
        
        try:
            # Set up ETKDG parameters
            # Use try/except for version compatibility across RDKit versions
            try:
                params = AllChem.ETKDGv3()
            except AttributeError:
                params = AllChem.ETKDG()
            
            params.randomSeed = self.random_seed
            
            # Set optional parameters if available (version-dependent)
            try:
                params.useRandomCoords = self.use_random_coords
            except AttributeError:
                pass
            
            try:
                params.pruneRmsThresh = self.prune_rms_thresh
            except AttributeError:
                pass
            
            # Generate conformers using ETKDG
            conf_ids = AllChem.EmbedMultipleConfs(
                mol,
                numConfs=self.num_conformers,
                params=params
            )
            
            if len(conf_ids) == 0:
                return ConformerResult(
                    success=False,
                    error="Conformer generation failed - no conformers embedded",
                    num_conformers_requested=self.num_conformers
                )
            
            logger.info(f"Embedded {len(conf_ids)} conformers for SMILES: {smiles[:50]}...")
            
            # Optimize with force field and calculate energies
            energies = self._optimize_conformers(mol, conf_ids)
            
            # Filter failed optimizations
            valid_confs = [(i, e) for i, e in enumerate(energies) if e is not None]
            
            if len(valid_confs) == 0:
                return ConformerResult(
                    success=False,
                    error="All force field optimizations failed",
                    num_conformers_requested=self.num_conformers
                )
            
            # Sort by energy (lowest first)
            valid_confs.sort(key=lambda x: x[1])
            lowest_conf_id = valid_confs[0][0]
            
            # Calculate 3D descriptors on lowest energy conformer
            descriptors_3d = self._calculate_3d_descriptors(mol, lowest_conf_id)
            
            # Cluster conformers by RMSD
            clusters = self._cluster_conformers_by_rmsd(
                mol, 
                [c[0] for c in valid_confs],
                threshold=self.rmsd_threshold
            )
            
            # Calculate energy statistics
            energies_valid = [e for _, e in valid_confs]
            lowest_energy = energies_valid[0]
            highest_energy = energies_valid[-1]
            energy_range = highest_energy - lowest_energy
            
            # Build conformer metadata
            conformer_metadata = {
                "lowest_energy_kcal_mol": round(lowest_energy, 3),
                "highest_energy_kcal_mol": round(highest_energy, 3),
                "energy_range_kcal_mol": round(energy_range, 3),
                "energy_mean_kcal_mol": round(np.mean(energies_valid), 3),
                "energy_std_kcal_mol": round(np.std(energies_valid), 3),
                "num_clusters_rmsd_2A": len(clusters),
                "cluster_sizes": [len(c) for c in clusters],
                "generation_method": "ETKDG",
                "optimization_method": self.force_field,
                "random_seed": self.random_seed
            }
            
            generation_time = time.time() - start_time
            
            return ConformerResult(
                success=True,
                num_conformers_generated=len(valid_confs),
                num_conformers_requested=self.num_conformers,
                lowest_energy=round(lowest_energy, 3),
                highest_energy=round(highest_energy, 3),
                energy_range=round(energy_range, 3),
                num_clusters_rmsd_2A=len(clusters),
                descriptors_3d=descriptors_3d,
                conformer_metadata=conformer_metadata,
                generation_time_seconds=round(generation_time, 2)
            )
            
        except Exception as e:
            logger.error(f"Conformer generation error: {str(e)}")
            return ConformerResult(
                success=False,
                error=f"Conformer generation error: {str(e)}",
                num_conformers_requested=self.num_conformers
            )
    
    def _optimize_conformers(self, mol, conf_ids: List[int]) -> List[Optional[float]]:
        """Optimize conformers with force field and return energies.
        
        Args:
            mol: RDKit molecule with conformers
            conf_ids: List of conformer IDs
            
        Returns:
            List of energies (None for failed optimizations)
        """
        energies = []
        
        for conf_id in conf_ids:
            try:
                if self.force_field == "MMFF94":
                    # Get MMFF properties
                    mmff_props = AllChem.MMFFGetMoleculeProperties(mol)
                    if mmff_props is None:
                        # Fall back to UFF if MMFF fails
                        converged = AllChem.UFFOptimizeMolecule(mol, confId=conf_id, maxIters=500)
                        if converged == 0:
                            ff = AllChem.UFFGetMoleculeForceField(mol, confId=conf_id)
                            if ff is not None:
                                energies.append(ff.CalcEnergy())
                            else:
                                energies.append(None)
                        else:
                            energies.append(None)
                    else:
                        # Optimize with MMFF94
                        converged = AllChem.MMFFOptimizeMolecule(
                            mol, 
                            mmffVariant="MMFF94",
                            confId=conf_id,
                            maxIters=500
                        )
                        if converged == 0:
                            ff = AllChem.MMFFGetMoleculeForceField(
                                mol, 
                                mmff_props, 
                                confId=conf_id
                            )
                            if ff is not None:
                                energies.append(ff.CalcEnergy())
                            else:
                                energies.append(None)
                        else:
                            energies.append(None)
                else:  # UFF
                    converged = AllChem.UFFOptimizeMolecule(mol, confId=conf_id, maxIters=500)
                    if converged == 0:
                        ff = AllChem.UFFGetMoleculeForceField(mol, confId=conf_id)
                        if ff is not None:
                            energies.append(ff.CalcEnergy())
                        else:
                            energies.append(None)
                    else:
                        energies.append(None)
                        
            except Exception as e:
                logger.warning(f"Optimization failed for conformer {conf_id}: {e}")
                energies.append(None)
        
        return energies
    
    def _calculate_3d_descriptors(self, mol, conf_id: int) -> Dict[str, float]:
        """Calculate 3D molecular descriptors.
        
        Args:
            mol: RDKit molecule with conformers
            conf_id: Conformer ID to calculate descriptors for
            
        Returns:
            Dictionary of 3D descriptor names and values
        """
        descriptors = {}
        
        try:
            # Principal Moments of Inertia (PMI)
            descriptors["pmi1"] = round(Descriptors3D.PMI1(mol, confId=conf_id), 4)
            descriptors["pmi2"] = round(Descriptors3D.PMI2(mol, confId=conf_id), 4)
            descriptors["pmi3"] = round(Descriptors3D.PMI3(mol, confId=conf_id), 4)
            
            # Normalized PMI ratios
            descriptors["npr1"] = round(Descriptors3D.NPR1(mol, confId=conf_id), 4)
            descriptors["npr2"] = round(Descriptors3D.NPR2(mol, confId=conf_id), 4)
            
            # Shape descriptors
            descriptors["radius_of_gyration"] = round(
                Descriptors3D.RadiusOfGyration(mol, confId=conf_id), 4
            )
            descriptors["inertial_shape_factor"] = round(
                Descriptors3D.InertialShapeFactor(mol, confId=conf_id), 4
            )
            descriptors["asphericity"] = round(
                Descriptors3D.Asphericity(mol, confId=conf_id), 4
            )
            descriptors["eccentricity"] = round(
                Descriptors3D.Eccentricity(mol, confId=conf_id), 4
            )
            descriptors["spherocity_index"] = round(
                Descriptors3D.SpherocityIndex(mol, confId=conf_id), 4
            )
            
            # Plane of Best Fit deviation
            try:
                descriptors["pbf"] = round(CalcPBF(mol, confId=conf_id), 4)
            except:
                descriptors["pbf"] = None
            
        except Exception as e:
            logger.warning(f"3D descriptor calculation error: {e}")
        
        return descriptors
    
    def _cluster_conformers_by_rmsd(
        self, 
        mol, 
        conf_ids: List[int], 
        threshold: float = 2.0
    ) -> List[List[int]]:
        """Cluster conformers by RMSD using a greedy algorithm.
        
        Args:
            mol: RDKit molecule with conformers
            conf_ids: List of conformer IDs (should be energy-sorted)
            threshold: RMSD threshold in Angstroms
            
        Returns:
            List of clusters, each cluster is a list of conformer IDs
        """
        if len(conf_ids) == 0:
            return []
        
        clusters = []
        assigned = set()
        
        for i, conf_id_i in enumerate(conf_ids):
            if conf_id_i in assigned:
                continue
            
            cluster = [conf_id_i]
            assigned.add(conf_id_i)
            
            for j, conf_id_j in enumerate(conf_ids[i+1:], start=i+1):
                if conf_id_j in assigned:
                    continue
                
                try:
                    rmsd = rdMolAlign.GetBestRMS(
                        mol, mol, 
                        prbId=conf_id_i, 
                        refId=conf_id_j
                    )
                    
                    if rmsd < threshold:
                        cluster.append(conf_id_j)
                        assigned.add(conf_id_j)
                except Exception as e:
                    logger.warning(f"RMSD calculation failed: {e}")
                    continue
            
            clusters.append(cluster)
        
        return clusters
    
    def generate_single_conformer(self, smiles: str) -> Optional[Dict]:
        """Generate a single optimized 3D conformer.
        
        Useful for quick 3D descriptor calculation without full ensemble.
        
        Args:
            smiles: SMILES string
            
        Returns:
            Dictionary with 3D descriptors or None if failed
        """
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        
        mol = Chem.AddHs(mol)
        
        try:
            # Embed single conformer
            AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
            
            # Optimize
            AllChem.MMFFOptimizeMolecule(mol)
            
            # Calculate descriptors
            return self._calculate_3d_descriptors(mol, 0)
            
        except Exception as e:
            logger.error(f"Single conformer generation failed: {e}")
            return None
    
    def get_conformer_sdf(self, smiles: str, num_conformers: int = 10) -> Optional[str]:
        """Generate conformers and return as SDF string.
        
        Useful for visualization or downstream docking.
        
        Args:
            smiles: SMILES string
            num_conformers: Number of conformers to include in SDF
            
        Returns:
            SDF format string with conformers or None if failed
        """
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        
        mol = Chem.AddHs(mol)
        
        try:
            # Generate conformers
            params = AllChem.ETKDGv3()
            params.randomSeed = self.random_seed
            AllChem.EmbedMultipleConfs(mol, numConfs=num_conformers, params=params)
            
            # Optimize
            for conf_id in range(mol.GetNumConformers()):
                AllChem.MMFFOptimizeMolecule(mol, confId=conf_id)
            
            # Write to SDF string
            from rdkit.Chem import SDWriter
            from io import StringIO
            
            sdf_buffer = StringIO()
            writer = SDWriter(sdf_buffer)
            
            for conf_id in range(mol.GetNumConformers()):
                writer.write(mol, confId=conf_id)
            
            writer.close()
            return sdf_buffer.getvalue()
            
        except Exception as e:
            logger.error(f"SDF generation failed: {e}")
            return None


class BatchConformerGenerator:
    """Batch processing for conformer generation.
    
    Handles multiple compounds efficiently with progress tracking.
    """
    
    def __init__(self, generator: Optional[ConformerGenerator] = None):
        """Initialize batch generator.
        
        Args:
            generator: ConformerGenerator instance (creates default if None)
        """
        self.generator = generator or ConformerGenerator()
        self.results: List[Dict] = []
    
    def process_batch(
        self, 
        compounds: List[Dict[str, str]],
        progress_callback: Optional[callable] = None
    ) -> Dict[str, Any]:
        """Process batch of compounds.
        
        Args:
            compounds: List of dicts with 'name' and 'smiles' keys
            progress_callback: Optional callback(current, total, compound_name)
            
        Returns:
            Summary dictionary with results
        """
        self.results = []
        total = len(compounds)
        success_count = 0
        failed_count = 0
        
        for i, compound in enumerate(compounds):
            name = compound.get("name", f"Compound_{i}")
            smiles = compound.get("smiles", "")
            
            if progress_callback:
                progress_callback(i + 1, total, name)
            
            logger.info(f"Processing {i+1}/{total}: {name}")
            
            result = self.generator.generate_conformers(smiles)
            
            result_dict = {
                "name": name,
                "smiles": smiles,
                "success": result.success,
                "num_conformers": result.num_conformers_generated,
                "descriptors_3d": result.descriptors_3d,
                "conformer_metadata": result.conformer_metadata,
                "error": result.error,
                "generation_time": result.generation_time_seconds
            }
            
            self.results.append(result_dict)
            
            if result.success:
                success_count += 1
            else:
                failed_count += 1
        
        return {
            "total_compounds": total,
            "success_count": success_count,
            "failed_count": failed_count,
            "success_rate": round(success_count / total * 100, 1) if total > 0 else 0,
            "results": self.results
        }
    
    def get_failed_compounds(self) -> List[Dict]:
        """Get list of compounds that failed conformer generation."""
        return [r for r in self.results if not r["success"]]
    
    def get_successful_compounds(self) -> List[Dict]:
        """Get list of compounds with successful conformer generation."""
        return [r for r in self.results if r["success"]]


# Convenience functions
def generate_conformers_for_smiles(smiles: str, num_conformers: int = 50) -> ConformerResult:
    """Quick conformer generation for a single SMILES.
    
    Args:
        smiles: SMILES string
        num_conformers: Number of conformers to generate
        
    Returns:
        ConformerResult object
    """
    generator = ConformerGenerator(num_conformers=num_conformers)
    return generator.generate_conformers(smiles)


def calculate_3d_descriptors(smiles: str) -> Optional[Dict]:
    """Quick 3D descriptor calculation for a SMILES.
    
    Args:
        smiles: SMILES string
        
    Returns:
        Dictionary of 3D descriptors or None if failed
    """
    generator = ConformerGenerator()
    return generator.generate_single_conformer(smiles)
