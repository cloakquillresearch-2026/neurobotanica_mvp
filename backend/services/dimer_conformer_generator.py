"""
Dimeric Conformer Generator Service
Extended 3D conformer generation for dimeric cannabinoid structures

Supports NeuroBotanica patent claims:
- 3D conformer generation for novel dimeric structures
- Intermonomer distance and orientation analysis
- Structural stability assessment
- Dimer-specific geometry optimization

Reference: NeuroBotanica MVP Development Plan - Week 4 Task 4.2
"""
import logging
from typing import Optional, List, Dict, Any, Tuple
from dataclasses import dataclass, field
from datetime import datetime

import numpy as np

logger = logging.getLogger(__name__)

# RDKit imports with fallback
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, Descriptors, Descriptors3D
    from rdkit.Chem import rdMolAlign, rdMolTransforms
    from rdkit.Chem.rdMolDescriptors import CalcPBF
    RDKIT_AVAILABLE = True
except ImportError:
    logger.warning("RDKit not available - dimeric conformer generation will be limited")
    RDKIT_AVAILABLE = False

from .conformer_generator import ConformerGenerator, ConformerResult


@dataclass
class DimerGeometryAnalysis:
    """Analysis of dimer-specific geometry features."""
    # Intermonomer metrics
    intermonomer_distance_A: float  # Distance between monomer centroids
    intermonomer_distance_range: Tuple[float, float]  # Min/max across conformers
    
    # Orientation metrics
    relative_orientation_deg: float  # Angle between monomer planes
    orientation_variance: float  # Variance across conformers
    
    # Structural stability
    stability_score: float  # 0.0-1.0 based on geometry
    strain_energy_estimate: float  # Estimated strain from suboptimal geometry
    
    # Linker analysis
    linker_strain: float  # Strain in the linker region
    linker_torsion_angles: List[float]  # Key torsion angles
    
    # Overall assessment
    geometry_quality: str  # "excellent", "good", "acceptable", "poor"
    geometry_warnings: List[str]  # Any detected issues


@dataclass
class DimericConformerResult:
    """Extended result for dimeric conformer generation."""
    # Base conformer result
    base_result: ConformerResult
    
    # Dimer-specific analysis
    dimer_analysis: Optional[DimerGeometryAnalysis] = None
    
    # Monomer identification
    monomer_1_atoms: List[int] = field(default_factory=list)
    monomer_2_atoms: List[int] = field(default_factory=list)
    linker_atoms: List[int] = field(default_factory=list)
    
    # Per-conformer dimer metrics
    conformer_dimer_metrics: List[Dict[str, float]] = field(default_factory=list)
    
    # Best conformers for different applications
    best_for_binding: int = 0  # Conf ID optimal for receptor binding
    best_for_stability: int = 0  # Conf ID most structurally stable
    best_overall: int = 0  # Combined score
    
    @property
    def success(self) -> bool:
        return self.base_result.success
    
    @property
    def num_conformers(self) -> int:
        return self.base_result.num_conformers_generated


class DimericConformerGenerator(ConformerGenerator):
    """Generate 3D conformers for dimeric cannabinoid structures.
    
    Extends ConformerGenerator with:
    - Dimer-specific geometry analysis
    - Intermonomer distance optimization
    - Relative orientation assessment
    - Structural stability scoring
    - Linker strain analysis
    
    Usage:
        generator = DimericConformerGenerator(num_conformers=30)
        result = generator.generate_dimer_conformers(
            dimer_smiles="...",
            monomer1_smiles="...",
            monomer2_smiles="..."
        )
    """
    
    # Optimal intermonomer distances (Angstroms)
    OPTIMAL_DISTANCES = {
        "methylene": (3.5, 5.5),  # Methylene bridge
        "ether": (3.0, 4.5),      # Ether linkage
        "direct": (2.5, 4.0),     # Direct C-C bond
        "default": (3.0, 6.0)
    }
    
    # Geometry quality thresholds
    STABILITY_THRESHOLDS = {
        "excellent": 0.85,
        "good": 0.70,
        "acceptable": 0.50,
        "poor": 0.0
    }
    
    def __init__(
        self,
        num_conformers: int = 30,
        max_attempts: int = 100,
        random_seed: int = 42,
        optimize_for_dimer: bool = True,
        **kwargs
    ):
        """Initialize dimeric conformer generator.
        
        Args:
            num_conformers: Number of conformers to generate
            max_attempts: Maximum embedding attempts
            random_seed: Random seed for reproducibility
            optimize_for_dimer: Whether to apply dimer-specific optimization
            **kwargs: Additional args passed to ConformerGenerator
        """
        super().__init__(
            num_conformers=num_conformers,
            max_attempts=max_attempts,
            random_seed=random_seed,
            **kwargs
        )
        self.optimize_for_dimer = optimize_for_dimer
        
        logger.info(f"DimericConformerGenerator initialized with {num_conformers} conformers")
    
    def generate_dimer_conformers(
        self,
        dimer_smiles: str,
        monomer1_smiles: Optional[str] = None,
        monomer2_smiles: Optional[str] = None,
        linkage_type: str = "default"
    ) -> DimericConformerResult:
        """Generate conformers for a dimeric structure with analysis.
        
        Args:
            dimer_smiles: SMILES of the dimer
            monomer1_smiles: SMILES of first monomer (for analysis)
            monomer2_smiles: SMILES of second monomer (for analysis)
            linkage_type: Type of linkage ("methylene", "ether", "direct", etc.)
            
        Returns:
            DimericConformerResult with conformers and dimer analysis
        """
        # Generate base conformers
        base_result = self.generate_conformers(dimer_smiles)
        
        if not base_result.success:
            return DimericConformerResult(
                base_result=base_result
            )
        
        # Get molecule for analysis
        mol = Chem.MolFromSmiles(dimer_smiles)
        mol = Chem.AddHs(mol)
        
        # Re-generate for analysis (we need the molecule with conformers)
        try:
            params = AllChem.ETKDGv3()
        except AttributeError:
            params = AllChem.ETKDG()
        params.randomSeed = self.random_seed
        AllChem.EmbedMultipleConfs(mol, numConfs=self.num_conformers, params=params)
        
        # Identify monomer regions
        monomer_1_atoms, monomer_2_atoms, linker_atoms = self._identify_monomer_regions(
            mol, monomer1_smiles, monomer2_smiles
        )
        
        # Analyze dimer geometry for each conformer
        conformer_metrics = []
        for conf_id in range(mol.GetNumConformers()):
            metrics = self._analyze_single_conformer_geometry(
                mol, conf_id,
                monomer_1_atoms, monomer_2_atoms, linker_atoms
            )
            conformer_metrics.append(metrics)
        
        # Calculate aggregate dimer analysis
        dimer_analysis = self._calculate_aggregate_analysis(
            conformer_metrics, linkage_type
        )
        
        # Find best conformers for different purposes
        best_binding, best_stability, best_overall = self._select_best_conformers(
            conformer_metrics
        )
        
        return DimericConformerResult(
            base_result=base_result,
            dimer_analysis=dimer_analysis,
            monomer_1_atoms=monomer_1_atoms,
            monomer_2_atoms=monomer_2_atoms,
            linker_atoms=linker_atoms,
            conformer_dimer_metrics=conformer_metrics,
            best_for_binding=best_binding,
            best_for_stability=best_stability,
            best_overall=best_overall
        )
    
    def _identify_monomer_regions(
        self,
        mol,
        monomer1_smiles: Optional[str],
        monomer2_smiles: Optional[str]
    ) -> Tuple[List[int], List[int], List[int]]:
        """Identify atoms belonging to each monomer and the linker."""
        if not RDKIT_AVAILABLE:
            return [], [], []
        
        total_atoms = mol.GetNumAtoms()
        
        if monomer1_smiles and monomer2_smiles:
            # Try to match monomers by substructure
            mol1 = Chem.MolFromSmiles(monomer1_smiles)
            mol2 = Chem.MolFromSmiles(monomer2_smiles)
            
            if mol1 and mol2:
                matches1 = mol.GetSubstructMatches(mol1)
                matches2 = mol.GetSubstructMatches(mol2)
                
                if matches1 and matches2:
                    monomer_1_atoms = list(matches1[0])
                    
                    # Find non-overlapping match for second monomer
                    for match in matches2:
                        if not set(match) & set(monomer_1_atoms):
                            monomer_2_atoms = list(match)
                            break
                    else:
                        monomer_2_atoms = list(matches2[0])
                    
                    # Linker atoms are everything else
                    all_monomer = set(monomer_1_atoms) | set(monomer_2_atoms)
                    linker_atoms = [i for i in range(total_atoms) if i not in all_monomer]
                    
                    return monomer_1_atoms, monomer_2_atoms, linker_atoms
        
        # Fallback: Split molecule roughly in half
        mid = total_atoms // 2
        monomer_1_atoms = list(range(mid))
        monomer_2_atoms = list(range(mid, total_atoms))
        linker_atoms = []
        
        return monomer_1_atoms, monomer_2_atoms, linker_atoms
    
    def _analyze_single_conformer_geometry(
        self,
        mol,
        conf_id: int,
        monomer_1_atoms: List[int],
        monomer_2_atoms: List[int],
        linker_atoms: List[int]
    ) -> Dict[str, float]:
        """Analyze geometry of a single conformer."""
        metrics = {}
        
        if not RDKIT_AVAILABLE or mol.GetNumConformers() == 0:
            return {"intermonomer_distance": 0.0, "stability_score": 0.5}
        
        conf = mol.GetConformer(conf_id)
        
        # Calculate centroids of each monomer
        centroid1 = self._calculate_centroid(conf, monomer_1_atoms)
        centroid2 = self._calculate_centroid(conf, monomer_2_atoms)
        
        # Intermonomer distance
        distance = np.linalg.norm(np.array(centroid1) - np.array(centroid2))
        metrics["intermonomer_distance"] = float(distance)
        
        # Calculate planes of best fit for each monomer
        if len(monomer_1_atoms) >= 3 and len(monomer_2_atoms) >= 3:
            normal1 = self._calculate_plane_normal(conf, monomer_1_atoms)
            normal2 = self._calculate_plane_normal(conf, monomer_2_atoms)
            
            # Relative orientation (angle between planes)
            if normal1 is not None and normal2 is not None:
                cos_angle = np.clip(np.dot(normal1, normal2), -1, 1)
                angle = np.degrees(np.arccos(abs(cos_angle)))
                metrics["relative_orientation"] = float(angle)
            else:
                metrics["relative_orientation"] = 90.0  # Default to perpendicular
        else:
            metrics["relative_orientation"] = 90.0
        
        # Linker torsion analysis
        if linker_atoms:
            torsions = self._analyze_linker_torsions(mol, conf, linker_atoms)
            metrics["linker_strain"] = torsions.get("strain", 0.0)
            metrics["linker_torsion"] = torsions.get("primary_torsion", 0.0)
        else:
            metrics["linker_strain"] = 0.0
            metrics["linker_torsion"] = 0.0
        
        # Calculate stability score
        stability = self._calculate_stability_score(metrics)
        metrics["stability_score"] = stability
        
        return metrics
    
    def _calculate_centroid(
        self,
        conf,
        atom_indices: List[int]
    ) -> Tuple[float, float, float]:
        """Calculate centroid of a set of atoms."""
        if not atom_indices:
            return (0.0, 0.0, 0.0)
        
        positions = []
        for idx in atom_indices:
            try:
                pos = conf.GetAtomPosition(idx)
                positions.append([pos.x, pos.y, pos.z])
            except (IndexError, RuntimeError):
                continue
        
        if not positions:
            return (0.0, 0.0, 0.0)
        
        centroid = np.mean(positions, axis=0)
        return tuple(centroid)
    
    def _calculate_plane_normal(
        self,
        conf,
        atom_indices: List[int]
    ) -> Optional[np.ndarray]:
        """Calculate normal vector to the plane of best fit."""
        if len(atom_indices) < 3:
            return None
        
        positions = []
        for idx in atom_indices[:20]:  # Limit to 20 atoms for efficiency
            try:
                pos = conf.GetAtomPosition(idx)
                positions.append([pos.x, pos.y, pos.z])
            except (IndexError, RuntimeError):
                continue
        
        if len(positions) < 3:
            return None
        
        # Center the points
        points = np.array(positions)
        centroid = np.mean(points, axis=0)
        centered = points - centroid
        
        # SVD to find plane normal
        try:
            _, _, vh = np.linalg.svd(centered)
            normal = vh[-1]  # Last row is normal to plane
            return normal / np.linalg.norm(normal)
        except np.linalg.LinAlgError:
            return None
    
    def _analyze_linker_torsions(
        self,
        mol,
        conf,
        linker_atoms: List[int]
    ) -> Dict[str, float]:
        """Analyze torsion angles in the linker region."""
        result = {"strain": 0.0, "primary_torsion": 0.0}
        
        if len(linker_atoms) < 2:
            return result
        
        # Find connected atoms
        torsion_angles = []
        for idx in linker_atoms:
            atom = mol.GetAtomWithIdx(idx)
            neighbors = [n.GetIdx() for n in atom.GetNeighbors()]
            
            # Look for rotatable bonds
            if len(neighbors) >= 2:
                for n_idx in neighbors:
                    if n_idx not in linker_atoms:
                        # This is a bond to a monomer
                        try:
                            # Get torsion angle if possible
                            n2_neighbors = [
                                n.GetIdx() for n in mol.GetAtomWithIdx(n_idx).GetNeighbors()
                                if n.GetIdx() != idx
                            ]
                            if n2_neighbors:
                                angle = rdMolTransforms.GetDihedralDeg(
                                    conf, n2_neighbors[0], n_idx, idx, neighbors[0]
                                )
                                torsion_angles.append(abs(angle))
                        except (IndexError, RuntimeError):
                            continue
        
        if torsion_angles:
            result["primary_torsion"] = float(torsion_angles[0])
            # Strain from deviation from ideal torsions
            ideal_torsions = [60, 180, 300]  # Gauche and anti
            min_deviations = []
            for angle in torsion_angles:
                min_dev = min(abs(angle - ideal) for ideal in ideal_torsions)
                min_deviations.append(min_dev)
            result["strain"] = float(np.mean(min_deviations) / 60)  # Normalize
        
        return result
    
    def _calculate_stability_score(self, metrics: Dict[str, float]) -> float:
        """Calculate overall stability score from geometry metrics."""
        score = 1.0
        
        # Distance penalty (too close or too far)
        distance = metrics.get("intermonomer_distance", 4.0)
        optimal_range = self.OPTIMAL_DISTANCES.get("default", (3.0, 6.0))
        
        if distance < optimal_range[0]:
            # Too close - steric clash
            score *= 0.5 * (distance / optimal_range[0])
        elif distance > optimal_range[1]:
            # Too far - weak interaction
            score *= max(0.5, 1 - (distance - optimal_range[1]) / 5)
        
        # Orientation bonus for favorable orientations
        orientation = metrics.get("relative_orientation", 90)
        # Parallel (0) or perpendicular (90) are both acceptable
        if orientation < 30 or orientation > 60:
            score *= 1.0  # Good orientation
        else:
            score *= 0.9  # Intermediate orientation
        
        # Linker strain penalty
        strain = metrics.get("linker_strain", 0.0)
        score *= max(0.5, 1 - strain * 0.3)
        
        return min(1.0, max(0.0, score))
    
    def _calculate_aggregate_analysis(
        self,
        conformer_metrics: List[Dict[str, float]],
        linkage_type: str
    ) -> DimerGeometryAnalysis:
        """Calculate aggregate analysis across all conformers."""
        if not conformer_metrics:
            return self._empty_analysis()
        
        # Extract metrics arrays
        distances = [m.get("intermonomer_distance", 0) for m in conformer_metrics]
        orientations = [m.get("relative_orientation", 0) for m in conformer_metrics]
        stabilities = [m.get("stability_score", 0) for m in conformer_metrics]
        linker_strains = [m.get("linker_strain", 0) for m in conformer_metrics]
        torsions = [m.get("linker_torsion", 0) for m in conformer_metrics]
        
        # Calculate aggregates
        avg_distance = float(np.mean(distances))
        distance_range = (float(min(distances)), float(max(distances)))
        avg_orientation = float(np.mean(orientations))
        orientation_var = float(np.var(orientations))
        avg_stability = float(np.mean(stabilities))
        avg_linker_strain = float(np.mean(linker_strains))
        
        # Estimate strain energy (simplified model)
        strain_energy = avg_linker_strain * 2.0  # kcal/mol estimate
        
        # Determine quality
        if avg_stability >= self.STABILITY_THRESHOLDS["excellent"]:
            quality = "excellent"
        elif avg_stability >= self.STABILITY_THRESHOLDS["good"]:
            quality = "good"
        elif avg_stability >= self.STABILITY_THRESHOLDS["acceptable"]:
            quality = "acceptable"
        else:
            quality = "poor"
        
        # Generate warnings
        warnings = []
        optimal_range = self.OPTIMAL_DISTANCES.get(linkage_type, self.OPTIMAL_DISTANCES["default"])
        
        if avg_distance < optimal_range[0]:
            warnings.append(f"Monomers may be too close ({avg_distance:.1f} Å) - possible steric clash")
        if avg_distance > optimal_range[1]:
            warnings.append(f"Monomers may be too far ({avg_distance:.1f} Å) - weak interaction")
        if orientation_var > 400:  # High variance in orientation
            warnings.append("High conformational flexibility - may affect binding specificity")
        if avg_linker_strain > 0.5:
            warnings.append("Significant linker strain detected - may affect stability")
        
        return DimerGeometryAnalysis(
            intermonomer_distance_A=avg_distance,
            intermonomer_distance_range=distance_range,
            relative_orientation_deg=avg_orientation,
            orientation_variance=orientation_var,
            stability_score=avg_stability,
            strain_energy_estimate=strain_energy,
            linker_strain=avg_linker_strain,
            linker_torsion_angles=torsions[:5],  # First 5 conformers
            geometry_quality=quality,
            geometry_warnings=warnings
        )
    
    def _select_best_conformers(
        self,
        conformer_metrics: List[Dict[str, float]]
    ) -> Tuple[int, int, int]:
        """Select best conformers for different applications."""
        if not conformer_metrics:
            return 0, 0, 0
        
        # Best for binding: moderate distance, good stability
        binding_scores = []
        for i, m in enumerate(conformer_metrics):
            dist = m.get("intermonomer_distance", 4.0)
            stab = m.get("stability_score", 0.5)
            # Optimal binding distance around 4-5 Å
            dist_score = 1 - abs(dist - 4.5) / 4
            binding_scores.append(dist_score * 0.4 + stab * 0.6)
        best_binding = int(np.argmax(binding_scores))
        
        # Best for stability: highest stability score
        stabilities = [m.get("stability_score", 0) for m in conformer_metrics]
        best_stability = int(np.argmax(stabilities))
        
        # Best overall: combined metric
        overall_scores = []
        for i, m in enumerate(conformer_metrics):
            score = (
                m.get("stability_score", 0) * 0.4 +
                binding_scores[i] * 0.3 +
                (1 - m.get("linker_strain", 0)) * 0.3
            )
            overall_scores.append(score)
        best_overall = int(np.argmax(overall_scores))
        
        return best_binding, best_stability, best_overall
    
    def _empty_analysis(self) -> DimerGeometryAnalysis:
        """Return empty analysis for failed generation."""
        return DimerGeometryAnalysis(
            intermonomer_distance_A=0.0,
            intermonomer_distance_range=(0.0, 0.0),
            relative_orientation_deg=0.0,
            orientation_variance=0.0,
            stability_score=0.0,
            strain_energy_estimate=0.0,
            linker_strain=0.0,
            linker_torsion_angles=[],
            geometry_quality="unknown",
            geometry_warnings=["Analysis failed"]
        )
