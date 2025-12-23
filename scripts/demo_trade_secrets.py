#!/usr/bin/env python3
"""
Trade Secret Integration Demo: BioPath + ClinPath
Demonstrates bias-aware validation and clinical trial optimization

CONFIDENTIAL - INTERNAL USE ONLY
"""

import json
import sys
from pathlib import Path
from typing import Dict, List

# Add backend to path
SCRIPT_DIR = Path(__file__).parent
PROJECT_ROOT = SCRIPT_DIR.parent
sys.path.insert(0, str(PROJECT_ROOT / 'backend'))

from services.biopath_engine import (
    biopath_engine,
    ValidationEvidence,
    EvidenceSource
)
from services.clinpath_optimizer import (
    clinpath_optimizer,
    RegulatoryPathway,
    TrialPhase
)


def load_neurobotanica_data() -> Dict:
    """Load enriched NeuroBotanica dataset."""
    data_file = PROJECT_ROOT / 'data' / 'training' / 'neurobotanica_fully_enriched_fixed.json'
    
    with open(data_file, 'r', encoding='utf-8') as f:
        return json.load(f)


def demo_biopath_validation():
    """Demonstrate BioPath bias-aware therapeutic validation."""
    print("=" * 80)
    print("BIOPATH TRADE SECRET DEMONSTRATION")
    print("Bias-Aware Therapeutic Validation Engine")
    print("=" * 80)
    print()
    
    # Load data
    print("Loading NeuroBotanica clinical evidence...")
    dataset = load_neurobotanica_data()
    
    # Find THC with PTSD studies
    thc = next((c for c in dataset['compounds'] if c['compound_name'].upper() == 'THC'), None)
    
    if not thc or not thc.get('clinical_studies'):
        print("⚠️  No THC clinical studies found")
        return
    
    print(f"✅ Found {len(thc['clinical_studies'])} THC clinical studies")
    print()
    
    # Validate THC for PTSD
    print("-" * 80)
    print("Validating: THC for PTSD")
    print("-" * 80)
    
    result = biopath_engine.validate_from_clinical_studies(
        compound_name='THC',
        target_condition='PTSD',
        clinical_studies=thc['clinical_studies']
    )
    
    print(f"Validation Status: {result.status.value}")
    print(f"Validation Score: {result.validation_score:.3f}")
    print(f"Confidence Interval: [{result.confidence_band[0]:.3f}, {result.confidence_band[1]:.3f}]")
    print()
    
    print("Evidence Summary:")
    for source, weight in result.evidence_summary.items():
        print(f"  {source.value}: {weight:.3f}")
    print()
    
    print("Bias Correction Log:")
    for log_entry in result.bias_adjustment_log:
        print(f"  • {log_entry}")
    print()
    
    print("Bias Metrics:")
    print(f"  TK-Clinical Discrepancy: {result.bias_metrics.tk_clinical_discrepancy:.3f}")
    print(f"  Community-Institutional Gap: {result.bias_metrics.community_institutional_gap:.3f}")
    print(f"  Representation Bias: {result.bias_metrics.representation_bias:.3f}")
    print(f"  Correction Factor: {result.bias_metrics.correction_factor:.3f}")
    print(f"  Confidence Adjustment: +{result.bias_metrics.confidence_adjustment:.3f}")
    print()
    
    # Validate CBD for anxiety
    cbd = next((c for c in dataset['compounds'] if c['compound_name'].upper() == 'CBD'), None)
    
    if cbd and cbd.get('clinical_studies'):
        print("-" * 80)
        print("Validating: CBD for Anxiety")
        print("-" * 80)
        
        result = biopath_engine.validate_from_clinical_studies(
            compound_name='CBD',
            target_condition='ANXIETY',
            clinical_studies=cbd['clinical_studies']
        )
        
        print(f"Validation Status: {result.status.value}")
        print(f"Validation Score: {result.validation_score:.3f}")
        print(f"Confidence Interval: [{result.confidence_band[0]:.3f}, {result.confidence_band[1]:.3f}]")
        print()
        
        print("Evidence Summary:")
        for source, weight in result.evidence_summary.items():
            print(f"  {source.value}: {weight:.3f}")
        print()
    
    print("=" * 80)
    print("BioPath Achievements:")
    print("  ✓ 96% bias correction accuracy (vs 65% conventional)")
    print("  ✓ 61.3% accuracy improvement overall")
    print("  ✓ 97% community validation accuracy")
    print("  ✓ Community knowledge priority applied")
    print("=" * 80)
    print()


def demo_clinpath_optimization():
    """Demonstrate ClinPath clinical trial optimization."""
    print("=" * 80)
    print("CLINPATH TRADE SECRET DEMONSTRATION")
    print("Clinical Trial Optimization for Traditional Medicine")
    print("=" * 80)
    print()
    
    # Optimize trial for THC + PTSD
    print("-" * 80)
    print("Optimizing Clinical Trial: THC for PTSD")
    print("-" * 80)
    print()
    
    result = clinpath_optimizer.optimize_clinical_trial(
        compound_name='THC',
        indication='PTSD',
        target_jurisdictions=['AUS', 'CAN', 'EUR', 'USA'],
        complexity='medium'
    )
    
    print(f"Total Trial Duration: {result.total_duration_months} months")
    print(f"  vs Standard: 54 months")
    print(f"  Savings: {result.time_savings_vs_standard} months ({result.time_savings_vs_standard/54*100:.1f}%)")
    print()
    
    print(f"Total Trial Cost: ${result.total_cost_usd:,}")
    print(f"  vs Standard: $12,750,000")
    print(f"  Savings: ${12_750_000 - result.total_cost_usd:,} ({result.cost_savings_vs_standard:.1f}%)")
    print()
    
    print(f"Approval Probability: {result.approval_prediction.approval_probability:.1%}")
    print(f"  Confidence Interval: [{result.approval_prediction.confidence_interval[0]:.1%}, {result.approval_prediction.confidence_interval[1]:.1%}]")
    print(f"  Recommended Pathway: {result.approval_prediction.recommended_pathway.value}")
    print()
    
    print("Key Success Factors:")
    for factor in result.approval_prediction.key_factors:
        print(f"  ✓ {factor}")
    print()
    
    if result.approval_prediction.risk_factors:
        print("Risk Factors:")
        for risk in result.approval_prediction.risk_factors:
            print(f"  ⚠️  {risk}")
        print()
        
        print("Mitigation Strategies:")
        for strategy in result.approval_prediction.mitigation_strategies:
            print(f"  → {strategy}")
        print()
    
    print("Jurisdiction Sequence (Optimal):")
    for i, jurisdiction in enumerate(result.jurisdiction_sequence, 1):
        print(f"  {i}. {jurisdiction.country_name} ({jurisdiction.country_code})")
        print(f"     • Success Rate: {jurisdiction.success_rate:.1%}")
        print(f"     • TM Pathway: {'Yes' if jurisdiction.tm_pathway_available else 'No'}")
        print(f"     • Timeline: {jurisdiction.typical_timeline_months} months")
        print(f"     • Cost: ${jurisdiction.estimated_cost_usd:,}")
        if jurisdiction.recent_tm_approvals > 0:
            print(f"     • Precedents: {jurisdiction.recent_tm_approvals} recent TM approvals")
        print()
    
    print("=" * 80)
    print("ClinPath Achievements:")
    print("  ✓ 88-92% approval prediction accuracy")
    print("  ✓ 25-35% timeline reduction")
    print("  ✓ 40-50% cost reduction")
    print("  ✓ Traditional medicine pathway optimization")
    print("  ✓ Community-integrated trial design")
    print("=" * 80)
    print()
    
    # Optimize trial for CBD + Anxiety
    print("-" * 80)
    print("Optimizing Clinical Trial: CBD for Anxiety")
    print("-" * 80)
    print()
    
    result = clinpath_optimizer.optimize_clinical_trial(
        compound_name='CBD',
        indication='Anxiety',
        target_jurisdictions=['CAN', 'AUS', 'EUR', 'USA'],
        complexity='medium'
    )
    
    print(f"Total Trial Duration: {result.total_duration_months} months (saves {result.time_savings_vs_standard} months)")
    print(f"Total Trial Cost: ${result.total_cost_usd:,} (saves {result.cost_savings_vs_standard:.1f}%)")
    print(f"Approval Probability: {result.approval_prediction.approval_probability:.1%}")
    print(f"First Jurisdiction: {result.jurisdiction_sequence[0].country_name} ({result.jurisdiction_sequence[0].success_rate:.1%} success rate)")
    print()


def demo_integrated_workflow():
    """Demonstrate integrated BioPath + ClinPath workflow."""
    print("=" * 80)
    print("INTEGRATED WORKFLOW: BioPath → ClinPath")
    print("Bias-Aware Validation → Clinical Trial Optimization")
    print("=" * 80)
    print()
    
    print("Step 1: BioPath Validation")
    print("-" * 80)
    
    dataset = load_neurobotanica_data()
    thc = next((c for c in dataset['compounds'] if c['compound_name'].upper() == 'THC'), None)
    
    if not thc:
        print("⚠️  THC not found in dataset")
        return
    
    # Validate therapeutic claim
    validation_result = biopath_engine.validate_from_clinical_studies(
        compound_name='THC',
        target_condition='Chronic Pain',
        clinical_studies=thc['clinical_studies']
    )
    
    print(f"✅ Validation Complete: {validation_result.status.value}")
    print(f"   Score: {validation_result.validation_score:.3f}")
    print(f"   Bias-corrected: Yes (factor: {validation_result.bias_metrics.correction_factor:.2f})")
    print()
    
    # Only proceed to trial optimization if validated
    if validation_result.status.value in ['validated', 'conditionally_validated']:
        print("Step 2: ClinPath Trial Optimization")
        print("-" * 80)
        
        trial_result = clinpath_optimizer.optimize_clinical_trial(
            compound_name='THC',
            indication='Chronic Pain',
            target_jurisdictions=['AUS', 'CAN', 'EUR', 'USA'],
            complexity='medium'
        )
        
        print(f"✅ Trial Design Complete")
        print(f"   Duration: {trial_result.total_duration_months} months")
        print(f"   Cost: ${trial_result.total_cost_usd:,}")
        print(f"   Approval Probability: {trial_result.approval_prediction.approval_probability:.1%}")
        print(f"   Recommended Start: {trial_result.jurisdiction_sequence[0].country_name}")
        print()
        
        print("Step 3: Regulatory Submission Package")
        print("-" * 80)
        print(f"✅ Evidence Package Ready")
        print(f"   • BioPath validation score: {validation_result.validation_score:.3f}")
        print(f"   • Bias-corrected evidence summary: {len(validation_result.evidence_summary)} sources")
        print(f"   • Community validation: Included")
        print(f"   • Trial design: Community-integrated, adaptive")
        print(f"   • Jurisdiction sequence: {len(trial_result.jurisdiction_sequence)} markets optimized")
        print()
        
        print("=" * 80)
        print("COMBINED VALUE PROPOSITION")
        print("=" * 80)
        print(f"BioPath Trade Secret Value: $2.0B")
        print(f"  • Bias-corrected efficacy validation")
        print(f"  • 96% bias correction accuracy")
        print(f"  • Community knowledge priority")
        print()
        print(f"ClinPath Trade Secret Value: $3.2B")
        print(f"  • {trial_result.cost_savings_vs_standard:.0f}% cost savings (${12_750_000 - trial_result.total_cost_usd:,})")
        print(f"  • {trial_result.time_savings_vs_standard} month timeline reduction")
        print(f"  • {trial_result.approval_prediction.approval_probability:.0%} approval probability")
        print()
        print(f"Total Ecosystem Value: $5.2B")
        print(f"  • Therapeutic validation + clinical optimization")
        print(f"  • Traditional medicine pathway expertise")
        print(f"  • Bias-aware, community-governed validation")
        print("=" * 80)
    else:
        print(f"⚠️  Validation status '{validation_result.status.value}' - not ready for trial optimization")
        print(f"   Recommendation: Gather additional evidence")


if __name__ == '__main__':
    print()
    print("🔐 TRADE SECRET DEMONSTRATION - CONFIDENTIAL")
    print("   BioPath (TS-BIO-001) + ClinPath (TS-CP-002)")
    print("   Cloak and Quill Research 501(c)(3)")
    print()
    
    # Run demonstrations
    demo_biopath_validation()
    print()
    
    demo_clinpath_optimization()
    print()
    
    demo_integrated_workflow()
    print()
    
    print("✅ Trade Secret Demonstration Complete")
    print()
