"""
GenomePath Data Pipeline - Correlation Generator
================================================

Generates synthetic TK↔Genomic correlations using trade secret algorithms
from correlation.py. Creates training pairs with confidence scores and
quality ratings.

Usage:
    python scripts/generate_correlations.py

Outputs:
    data/processed/training_correlations.json - TK↔Genomic correlation pairs
    data/processed/correlation_statistics.json - Quality metrics and distributions
"""

import json
import os
import sys
from pathlib import Path
from typing import Dict, List, Tuple
from datetime import datetime
from collections import defaultdict

HIGH_EVIDENCE_LEVELS = {
    "published_knowledge",
    "peer_reviewed",
    "clinical_research",
    "clinical_trials",
    "modern_research",
}

# Add backend to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from backend.services.genomepath.bridge import (
    GenomePathBridge,
    TKEncoder,
    GenomicSequenceEncoder
)
from backend.services.genomepath.correlation import (
    TKGenomicCorrelator,
    CorrelationQuality
)


class CorrelationGenerator:
    """Generates synthetic TK↔Genomic correlations for training."""
    
    def __init__(
        self,
        tk_practices_path: str = "data/processed/tk_practices.json",
        genomic_targets_path: str = "data/processed/genomic_targets.json"
    ):
        self.tk_practices_path = Path(tk_practices_path)
        self.genomic_targets_path = Path(genomic_targets_path)
        
        self.tk_practices = []
        self.practice_by_id: Dict[str, Dict] = {}
        self.genomic_targets = {}
        self.correlations = []
        
        self.bridge = GenomePathBridge()
        self.correlator = TKGenomicCorrelator()
        
        self.stats = {
            "total_correlations": 0,
            "tk_to_genomic": 0,
            "genomic_to_tk": 0,
            "quality_distribution": defaultdict(int),
            "average_confidence": 0.0,
            "bidirectional_consistent": 0
        }
    
    def load_tk_practices(self):
        """Load TK practice dataset."""
        
        if not self.tk_practices_path.exists():
            raise FileNotFoundError(
                f"TK practices not found: {self.tk_practices_path}\n"
                "Run scripts/build_tk_dataset.py first"
            )
        
        with open(self.tk_practices_path, 'r', encoding='utf-8') as f:
            data = json.load(f)
            self.tk_practices = data.get("practices", [])
            self.practice_by_id = {
                practice.get("practice_id"): practice
                for practice in self.tk_practices
                if practice.get("practice_id")
            }
        
        print(f"📚 Loaded {len(self.tk_practices)} TK practices")
    
    def load_genomic_targets(self):
        """Load genomic target database."""
        
        if not self.genomic_targets_path.exists():
            raise FileNotFoundError(
                f"Genomic targets not found: {self.genomic_targets_path}\n"
                "Run scripts/extract_genomic_targets.py first"
            )
        
        with open(self.genomic_targets_path, 'r', encoding='utf-8') as f:
            data = json.load(f)
            self.genomic_targets = data.get("genes", {})
        
        print(f"🧬 Loaded {len(self.genomic_targets)} genomic targets")
    
    def generate_tk_to_genomic_correlations(self):
        """Generate TK → Genomic hypothesis correlations."""
        
        print()
        print("🔄 Generating TK → Genomic correlations...")
        
        for practice in self.tk_practices:
            # Skip sacred knowledge
            if practice.get("ceremonial_significance"):
                print(f"   ⚠️  Skipping sacred practice: {practice['practice_name']}")
                continue

            # Skip practices flagged for additional metadata or consent
            if practice.get("suppress_from_training"):
                print(f"   ⏭️  Skipping flagged practice: {practice['practice_name']}")
                continue

            evidence_level = (practice.get("evidence_level") or "").lower()
            high_dose_threshold = 0.48 if evidence_level in HIGH_EVIDENCE_LEVELS else 0.60
            
            try:
                # Encode TK practice
                tk_vector, bridge_result = self.bridge.correlate_tk_to_genomic(
                    practice_name=practice["practice_name"],
                    source_community_id=practice["source_community_id"],
                    knowledge_domain=practice["knowledge_domain"],
                    preparation_method=practice.get("preparation_method", ""),
                    indications=practice.get("indications", []),
                    cultural_context=practice.get("cultural_context", ""),
                    ceremonial_significance=practice.get("ceremonial_significance", False)
                )
                
                # Generate genomic hypotheses
                correlation_result = self.correlator.correlate_tk_to_genomic(
                    tk_vector=tk_vector,
                    bridge_result=bridge_result,
                    traditional_indications=practice.get("indications", [])
                )
                
                # Store correlations with dosage variations
                # Dosage-specific target modulation (low/medium/high)
                dosage_profiles = [
                    {
                        'dosage': 'low',
                        'range': '2.5-5mg THC',
                        'confidence_modifier': 0.05,  # Better studied
                        'target_limit': 5  # Low dose: fewer targets
                    },
                    {
                        'dosage': 'medium',
                        'range': '10-25mg THC',
                        'confidence_modifier': 0.0,  # Baseline
                        'target_limit': 8  # Medium: broader engagement
                    },
                    {
                        'dosage': 'high',
                        'range': '50mg+ THC',
                        'confidence_modifier': -0.03,  # Less studied
                        'target_limit': 10  # High: all targets
                    }
                ]
                
                # Create dosage-specific correlations
                sorted_hypotheses = sorted(
                    correlation_result.genomic_hypotheses,
                    key=lambda h: h.overall_confidence,
                    reverse=True
                )
                
                for profile in dosage_profiles:
                    for hypothesis in sorted_hypotheses[:profile['target_limit']]:
                        # Gate medium/high dose variants when base confidence is weak
                        if (
                            hypothesis.overall_confidence < high_dose_threshold and
                            profile['dosage'] in {"medium", "high"}
                        ):
                            continue

                        # Adjust confidence based on dosage profile
                        adjusted_confidence = min(
                            hypothesis.overall_confidence + profile['confidence_modifier'],
                            0.95
                        )
                        
                        correlation = {
                            "correlation_id": f"TK2G_{len(self.correlations):05d}",
                            "direction": "tk_to_genomic",
                            "tk_practice_id": practice.get("practice_id"),
                            "tk_practice_name": practice["practice_name"],
                            "source_community": practice["source_community_id"],
                            "genomic_target": hypothesis.target_gene_id,
                            "target_type": hypothesis.target_type.value,
                            "confidence": adjusted_confidence,
                            "quality": hypothesis.correlation_quality.value,
                            "mechanism": hypothesis.predicted_mechanism.value if hypothesis.predicted_mechanism else None,
                            "tissue_confidence": hypothesis.tissue_expression_confidence,
                            "pathway_confidence": hypothesis.pathway_confidence,
                            "disease_confidence": hypothesis.disease_association_confidence,
                            "literature_confidence": hypothesis.literature_support_confidence,
                            "dosage": profile['dosage'],
                            "dosage_range": profile['range'],
                            "evidence_sources": practice.get("literature_sources", []),
                            "generated_date": datetime.now().isoformat()
                        }
                        
                        self.correlations.append(correlation)
                        self.stats["tk_to_genomic"] += 1
                        quality_key = hypothesis.correlation_quality.value.upper()
                        self.stats["quality_distribution"][quality_key] += 1
                
                print(f"   ✅ {practice['practice_name']}: {len(correlation_result.genomic_hypotheses)} genomic hypotheses")
                
            except ValueError as e:
                print(f"   ⚠️  Error processing {practice['practice_name']}: {e}")
            except Exception as e:
                print(f"   ❌ Unexpected error: {e}")
    
    def generate_genomic_to_tk_correlations(self):
        """Generate Genomic → TK practice correlations."""
        
        print()
        print("🔄 Generating Genomic → TK correlations...")
        
        # Select representative genes (those with most condition associations)
        top_genes = sorted(
            self.genomic_targets.items(),
            key=lambda x: x[1].get("total_studies", 0),
            reverse=True
        )[:46]  # All genes for maximum correlation coverage
        
        for gene_id, gene_data in top_genes:
            try:
                # Build tissue expression and pathway data
                tissues = gene_data.get("tissues", [])
                pathways = gene_data.get("pathways", [])
                
                # Get known TK correlations from associated conditions
                known_tk = []
                for condition_data in gene_data.get("associated_conditions", []):
                    condition = condition_data.get("condition")
                    # Find TK practices for this condition
                    for practice in self.tk_practices:
                        if condition in practice.get("indications", []):
                            known_tk.append(practice.get("practice_id", ""))
                
                # Encode genomic sequence
                genomic_vector, bridge_result = self.bridge.correlate_genomic_to_tk(
                    gene_id=gene_id,
                    tissue_expression=tissues[:8],  # Increased from 5 to 8 tissues
                    pathway_involvement=pathways[:8],  # Increased from 5 to 8 pathways
                    known_tk_correlations=known_tk[:6]  # Increased from 3 to 6 TK practices
                )
                
                # Generate TK practice correlations
                correlation_result = self.correlator.correlate_genomic_to_tk(
                    genomic_vector=genomic_vector,
                    bridge_result=bridge_result
                )
                
                # Store correlations
                for tk_corr in correlation_result.traditional_correlations:
                    practice_id = tk_corr.correlated_practice_id
                    practice_meta = self.practice_by_id.get(practice_id)
                    practice_name = (
                        practice_meta.get("practice_name")
                        if practice_meta else tk_corr.correlated_practice_name
                    )
                    practice_community = (
                        practice_meta.get("source_community_id")
                        if practice_meta else tk_corr.correlated_community_id
                    )
                    indications = (
                        practice_meta.get("indications", [])
                        if practice_meta else tk_corr.traditional_indications
                    )
                    preparation_notes = (
                        practice_meta.get("preparation_method")
                        if practice_meta else None
                    )
                    correlation = {
                        "correlation_id": f"G2TK_{len(self.correlations):05d}",
                        "direction": "genomic_to_tk",
                        "genomic_target": gene_id,
                        "tk_practice_id": practice_id,
                        "tk_practice_predicted": practice_name,
                        "predicted_community": practice_community,
                        "confidence": tk_corr.overall_confidence,
                        "quality": tk_corr.correlation_quality.value,
                        "requires_validation": tk_corr.community_validation_required,
                        "traditional_indications": indications,
                        "recommended_preparation": preparation_notes,
                        "evidence_sources": gene_data.get("evidence_sources", []),
                        "generated_date": datetime.now().isoformat()
                    }
                    
                    self.correlations.append(correlation)
                    self.stats["genomic_to_tk"] += 1
                    quality_key = tk_corr.correlation_quality.value.upper()
                    self.stats["quality_distribution"][quality_key] += 1
                
                print(f"   ✅ {gene_id}: {len(correlation_result.traditional_correlations)} TK practice predictions")
                
            except Exception as e:
                print(f"   ⚠️  Error processing {gene_id}: {e}")
    
    def calculate_statistics(self):
        """Calculate correlation statistics and quality metrics."""
        
        self.stats["total_correlations"] = len(self.correlations)
        
        if self.correlations:
            confidences = [c["confidence"] for c in self.correlations]
            self.stats["average_confidence"] = sum(confidences) / len(confidences)
            self.stats["min_confidence"] = min(confidences)
            self.stats["max_confidence"] = max(confidences)
        
        # Count bidirectional consistency
        tk_to_genomic_pairs = {
            (c["tk_practice_name"], c["genomic_target"])
            for c in self.correlations if c["direction"] == "tk_to_genomic"
        }
        
        genomic_to_tk_pairs = {
            (c["tk_practice_predicted"], c["genomic_target"])
            for c in self.correlations if c["direction"] == "genomic_to_tk"
        }
        
        self.stats["bidirectional_consistent"] = len(
            tk_to_genomic_pairs.intersection(genomic_to_tk_pairs)
        )
    
    def export_correlations(self, output_path: str = "data/processed/training_correlations.json"):
        """Export correlation dataset."""
        
        output = {
            "metadata": {
                "version": "1.0",
                "created_date": datetime.now().isoformat(),
                "total_correlations": self.stats["total_correlations"],
                "tk_to_genomic": self.stats["tk_to_genomic"],
                "genomic_to_tk": self.stats["genomic_to_tk"],
                "average_confidence": round(self.stats["average_confidence"], 4),
                "bidirectional_consistent": self.stats["bidirectional_consistent"],
                "quality_distribution": dict(self.stats["quality_distribution"]),
                "generation_method": "GenomePath TS-GP-001 trade secret algorithms"
            },
            "correlations": self.correlations
        }
        
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        
        with open(output_path, 'w', encoding='utf-8') as f:
            json.dump(output, f, indent=2, ensure_ascii=False)
        
        print(f"✅ Exported {self.stats['total_correlations']} correlations to {output_path}")
    
    def export_statistics(self, output_path: str = "data/processed/correlation_statistics.json"):
        """Export detailed statistics."""
        
        statistics = {
            "generation_date": datetime.now().isoformat(),
            "overall_statistics": self.stats,
            "quality_breakdown": {
                quality: {
                    "count": self.stats["quality_distribution"][quality],
                    "percentage": round(
                        100 * self.stats["quality_distribution"][quality] / self.stats["total_correlations"],
                        2
                    ) if self.stats["total_correlations"] > 0 else 0
                }
                for quality in ["EXCELLENT", "GOOD", "MODERATE", "POOR"]
            },
            "direction_breakdown": {
                "tk_to_genomic": {
                    "count": self.stats["tk_to_genomic"],
                    "percentage": round(
                        100 * self.stats["tk_to_genomic"] / self.stats["total_correlations"],
                        2
                    ) if self.stats["total_correlations"] > 0 else 0
                },
                "genomic_to_tk": {
                    "count": self.stats["genomic_to_tk"],
                    "percentage": round(
                        100 * self.stats["genomic_to_tk"] / self.stats["total_correlations"],
                        2
                    ) if self.stats["total_correlations"] > 0 else 0
                }
            }
        }
        
        with open(output_path, 'w', encoding='utf-8') as f:
            json.dump(statistics, f, indent=2, ensure_ascii=False)
        
        print(f"✅ Exported statistics to {output_path}")
    
    def print_summary(self):
        """Print generation summary."""
        
        print()
        print("=" * 70)
        print("Correlation Generation Summary")
        print("=" * 70)
        print(f"Total correlations: {self.stats['total_correlations']}")
        print(f"  TK → Genomic: {self.stats['tk_to_genomic']}")
        print(f"  Genomic → TK: {self.stats['genomic_to_tk']}")
        print(f"Bidirectional consistent: {self.stats['bidirectional_consistent']}")
        print(f"Average confidence: {self.stats['average_confidence']:.4f}")
        print()
        print("Quality Distribution:")
        for quality, count in sorted(self.stats["quality_distribution"].items()):
            percentage = 100 * count / self.stats["total_correlations"] if self.stats["total_correlations"] > 0 else 0
            print(f"  {quality}: {count} ({percentage:.1f}%)")
        print()


def main():
    """Main correlation generation pipeline."""
    
    print("=" * 70)
    print("GenomePath Correlation Generator")
    print("=" * 70)
    print()
    
    generator = CorrelationGenerator()
    
    # Load data
    print("📂 Loading datasets...")
    generator.load_tk_practices()
    generator.load_genomic_targets()
    
    # Generate correlations
    generator.generate_tk_to_genomic_correlations()
    generator.generate_genomic_to_tk_correlations()
    
    # Calculate statistics
    print()
    print("📊 Calculating statistics...")
    generator.calculate_statistics()
    
    # Export results
    print()
    print("💾 Exporting results...")
    generator.export_correlations()
    generator.export_statistics()
    
    # Print summary
    generator.print_summary()
    
    print("=" * 70)
    print("✅ Correlation generation complete!")
    print("=" * 70)
    print()
    print("Next steps:")
    print("  1. Review data/processed/training_correlations.json")
    print("  2. Run scripts/validate_dataset.py to check quality")
    print("  3. Manually validate top 50 correlations against literature")
    print()


if __name__ == "__main__":
    main()
