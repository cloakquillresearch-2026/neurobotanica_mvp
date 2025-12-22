"""
GenomePath Data Pipeline - Dataset Validator
============================================

Validates data completeness, consistency, and quality across all datasets.
Checks for sacred knowledge violations, attribution integrity, and
correlation quality metrics.

Usage:
    python scripts/validate_dataset.py

Outputs:
    data/processed/validation_report.json - Comprehensive validation report
"""

import json
import os
from pathlib import Path
from typing import Dict, List, Set
from datetime import datetime
from collections import defaultdict


class DatasetValidator:
    """Validates GenomePath datasets for quality and compliance."""
    
    def __init__(self):
        self.tk_practices = []
        self.genomic_targets = {}
        self.correlations = []
        
        self.validation_results = {
            "validation_date": datetime.now().isoformat(),
            "overall_status": "PASS",
            "errors": [],
            "warnings": [],
            "sacred_knowledge_checks": {},
            "attribution_checks": {},
            "quality_metrics": {},
            "completeness_checks": {}
        }
    
    def load_datasets(self):
        """Load all datasets for validation."""
        
        # Load TK practices
        tk_path = Path("data/processed/tk_practices.json")
        if tk_path.exists():
            with open(tk_path, 'r', encoding='utf-8') as f:
                data = json.load(f)
                self.tk_practices = data.get("practices", [])
            print(f"✅ Loaded {len(self.tk_practices)} TK practices")
        else:
            self.validation_results["warnings"].append("TK practices file not found")
            print("⚠️  TK practices not found")
        
        # Load genomic targets
        genomic_path = Path("data/processed/genomic_targets.json")
        if genomic_path.exists():
            with open(genomic_path, 'r', encoding='utf-8') as f:
                data = json.load(f)
                self.genomic_targets = data.get("genes", {})
            print(f"✅ Loaded {len(self.genomic_targets)} genomic targets")
        else:
            self.validation_results["warnings"].append("Genomic targets file not found")
            print("⚠️  Genomic targets not found")
        
        # Load correlations
        corr_path = Path("data/processed/training_correlations.json")
        if corr_path.exists():
            with open(corr_path, 'r', encoding='utf-8') as f:
                data = json.load(f)
                self.correlations = data.get("correlations", [])
            print(f"✅ Loaded {len(self.correlations)} correlations")
        else:
            self.validation_results["warnings"].append("Correlations file not found")
            print("⚠️  Correlations not found")
    
    def validate_sacred_knowledge_protection(self):
        """Validate sacred knowledge is properly flagged and protected."""
        
        print()
        print("🔒 Validating sacred knowledge protection...")
        
        sacred_practices = [
            p for p in self.tk_practices 
            if p.get("ceremonial_significance") is True
        ]
        
        violations = []
        for practice in sacred_practices:
            # Check consent status
            if practice.get("community_consent_status") != "explicit_consent_granted":
                violations.append({
                    "practice_id": practice.get("practice_id"),
                    "practice_name": practice.get("practice_name"),
                    "issue": "Sacred practice without explicit consent - MUST NOT BE USED"
                })
            
            # Check if practice appears in correlations
            practice_id = practice.get("practice_id")
            for corr in self.correlations:
                if corr.get("tk_practice_id") == practice_id:
                    violations.append({
                        "practice_id": practice_id,
                        "correlation_id": corr.get("correlation_id"),
                        "issue": "CRITICAL: Sacred practice used in correlation generation"
                    })
        
        self.validation_results["sacred_knowledge_checks"] = {
            "total_sacred_practices": len(sacred_practices),
            "violations": violations,
            "status": "FAIL" if violations else "PASS"
        }
        
        if violations:
            self.validation_results["overall_status"] = "FAIL"
            self.validation_results["errors"].extend([
                f"Sacred knowledge violation: {v['issue']}" for v in violations
            ])
            print(f"   ❌ {len(violations)} CRITICAL violations found")
        else:
            print(f"   ✅ All {len(sacred_practices)} sacred practices properly protected")
    
    def validate_community_attribution(self):
        """Validate all practices have proper community attribution."""
        
        print()
        print("📝 Validating community attribution...")
        
        missing_attribution = []
        for practice in self.tk_practices:
            if not practice.get("source_community_id"):
                missing_attribution.append({
                    "practice_id": practice.get("practice_id"),
                    "practice_name": practice.get("practice_name"),
                    "issue": "Missing source_community_id"
                })
            
            if not practice.get("attribution_notes"):
                missing_attribution.append({
                    "practice_id": practice.get("practice_id"),
                    "practice_name": practice.get("practice_name"),
                    "issue": "Missing attribution_notes"
                })
        
        self.validation_results["attribution_checks"] = {
            "total_practices": len(self.tk_practices),
            "missing_attribution": missing_attribution,
            "status": "FAIL" if missing_attribution else "PASS"
        }
        
        if missing_attribution:
            self.validation_results["warnings"].extend([
                f"Missing attribution: {a['practice_name']} - {a['issue']}"
                for a in missing_attribution
            ])
            print(f"   ⚠️  {len(missing_attribution)} practices with incomplete attribution")
        else:
            print(f"   ✅ All {len(self.tk_practices)} practices properly attributed")
    
    def validate_correlation_quality(self):
        """Validate correlation quality metrics and distributions."""
        
        print()
        print("📊 Validating correlation quality...")
        
        quality_dist = defaultdict(int)
        confidence_scores = []
        low_quality_correlations = []
        
        for corr in self.correlations:
            quality = corr.get("quality", "UNKNOWN")
            confidence = corr.get("confidence", 0.0)
            
            quality_dist[quality] += 1
            confidence_scores.append(confidence)
            
            # Flag low quality correlations
            if quality == "POOR" or confidence < 0.60:
                low_quality_correlations.append({
                    "correlation_id": corr.get("correlation_id"),
                    "quality": quality,
                    "confidence": confidence,
                    "direction": corr.get("direction")
                })
        
        avg_confidence = sum(confidence_scores) / len(confidence_scores) if confidence_scores else 0
        
        self.validation_results["quality_metrics"] = {
            "total_correlations": len(self.correlations),
            "average_confidence": round(avg_confidence, 4),
            "quality_distribution": dict(quality_dist),
            "low_quality_count": len(low_quality_correlations),
            "low_quality_correlations": low_quality_correlations[:20],  # Top 20
            "status": "PASS" if avg_confidence >= 0.65 else "WARNING"
        }
        
        print(f"   Average confidence: {avg_confidence:.4f}")
        print(f"   Quality distribution: {dict(quality_dist)}")
        
        if avg_confidence < 0.65:
            self.validation_results["warnings"].append(
                f"Average confidence ({avg_confidence:.4f}) below recommended threshold (0.65)"
            )
            print(f"   ⚠️  Average confidence below recommended threshold")
        else:
            print(f"   ✅ Correlation quality meets standards")
    
    def validate_completeness(self):
        """Validate dataset completeness and coverage."""
        
        print()
        print("📋 Validating dataset completeness...")
        
        # Check TK practice coverage
        tk_conditions = set()
        for practice in self.tk_practices:
            tk_conditions.update(practice.get("indications", []))
        
        # Check genomic target coverage
        genomic_conditions = set()
        for gene_data in self.genomic_targets.values():
            for condition_data in gene_data.get("associated_conditions", []):
                genomic_conditions.add(condition_data.get("condition"))
        
        # Check correlation coverage
        corr_tk_practices = {c.get("tk_practice_id") for c in self.correlations}
        corr_genomic_targets = {c.get("genomic_target") for c in self.correlations}
        
        missing_tk_in_corr = []
        for practice in self.tk_practices:
            if not practice.get("ceremonial_significance"):  # Exclude sacred
                if practice.get("practice_id") not in corr_tk_practices:
                    missing_tk_in_corr.append(practice.get("practice_name"))
        
        self.validation_results["completeness_checks"] = {
            "tk_practices": {
                "total": len(self.tk_practices),
                "conditions_covered": len(tk_conditions),
                "conditions": sorted(list(tk_conditions))
            },
            "genomic_targets": {
                "total": len(self.genomic_targets),
                "conditions_covered": len(genomic_conditions),
                "conditions": sorted(list(genomic_conditions))
            },
            "correlations": {
                "total": len(self.correlations),
                "tk_practices_covered": len(corr_tk_practices),
                "genomic_targets_covered": len(corr_genomic_targets),
                "missing_tk_practices": missing_tk_in_corr[:10]  # Top 10
            },
            "status": "PASS"
        }
        
        print(f"   TK practices: {len(self.tk_practices)} ({len(tk_conditions)} conditions)")
        print(f"   Genomic targets: {len(self.genomic_targets)} ({len(genomic_conditions)} conditions)")
        print(f"   Correlations: {len(self.correlations)}")
        print(f"   Coverage: {len(corr_tk_practices)}/{len(self.tk_practices)} TK practices")
        
        if missing_tk_in_corr:
            print(f"   ⚠️  {len(missing_tk_in_corr)} TK practices not in correlations")
        else:
            print(f"   ✅ Complete coverage achieved")
    
    def export_validation_report(self, output_path: str = "data/processed/validation_report.json"):
        """Export validation report."""
        
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        
        with open(output_path, 'w', encoding='utf-8') as f:
            json.dump(self.validation_results, f, indent=2, ensure_ascii=False)
        
        print()
        print(f"✅ Validation report exported to {output_path}")
    
    def print_summary(self):
        """Print validation summary."""
        
        print()
        print("=" * 70)
        print("Validation Summary")
        print("=" * 70)
        print(f"Overall Status: {self.validation_results['overall_status']}")
        print()
        
        if self.validation_results["errors"]:
            print(f"❌ ERRORS ({len(self.validation_results['errors'])}):")
            for error in self.validation_results["errors"]:
                print(f"   - {error}")
            print()
        
        if self.validation_results["warnings"]:
            print(f"⚠️  WARNINGS ({len(self.validation_results['warnings'])}):")
            for warning in self.validation_results["warnings"]:
                print(f"   - {warning}")
            print()
        
        print("Check Results:")
        print(f"  Sacred Knowledge Protection: {self.validation_results['sacred_knowledge_checks'].get('status', 'NOT RUN')}")
        print(f"  Community Attribution: {self.validation_results['attribution_checks'].get('status', 'NOT RUN')}")
        print(f"  Correlation Quality: {self.validation_results['quality_metrics'].get('status', 'NOT RUN')}")
        print(f"  Dataset Completeness: {self.validation_results['completeness_checks'].get('status', 'NOT RUN')}")
        print()


def main():
    """Main validation pipeline."""
    
    print("=" * 70)
    print("GenomePath Dataset Validator")
    print("=" * 70)
    print()
    
    validator = DatasetValidator()
    
    # Load datasets
    print("📂 Loading datasets...")
    validator.load_datasets()
    
    # Run validation checks
    validator.validate_sacred_knowledge_protection()
    validator.validate_community_attribution()
    validator.validate_correlation_quality()
    validator.validate_completeness()
    
    # Export report
    validator.export_validation_report()
    
    # Print summary
    validator.print_summary()
    
    print("=" * 70)
    
    if validator.validation_results["overall_status"] == "FAIL":
        print("❌ VALIDATION FAILED")
        print("=" * 70)
        print()
        print("CRITICAL: Fix errors before using this data for training/deployment")
        return 1
    else:
        print("✅ VALIDATION PASSED")
        print("=" * 70)
        print()
        print("Dataset ready for training and deployment")
        print()
        print("Recommended next steps:")
        print("  1. Manually validate top 50 correlations against literature")
        print("  2. Review any warnings and improve data quality")
        print("  3. Add more TK practices to increase coverage")
        print("  4. Consider expanding genomic target database")
        return 0


if __name__ == "__main__":
    exit(main())
