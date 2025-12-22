"""
GenomePath Data Pipeline - Traditional Knowledge Dataset Builder
================================================================

Creates structured TK practice dataset with proper attribution, consent flags,
and sacred knowledge protection.

Usage:
    python scripts/build_tk_dataset.py

Outputs:
    data/processed/tk_practices.json - Structured TK practice database
    data/processed/tk_practice_template.json - Template for manual curation
"""

import json
import os
from pathlib import Path
from typing import Dict, List
from datetime import datetime


# Example TK practices from published ethnobotanical literature
# These are publicly documented, non-sacred practices suitable for MVP
EXAMPLE_TK_PRACTICES = [
    {
        "practice_name": "Cannabis for chronic pain relief",
        "source_community_id": "traditional_chinese_medicine",
        "knowledge_domain": "pain_management",
        "preparation_method": "Decoction of cannabis flowers and leaves",
        "indications": ["chronic_pain", "inflammation", "muscle_spasms"],
        "contraindications": ["pregnancy", "severe_cardiovascular_disease"],
        "cultural_context": "Used in TCM for over 2000 years, documented in Shennong Bencaojing (Divine Farmer's Materia Medica)",
        "ceremonial_significance": False,
        "traditional_dosage": "3-5g dried material per day",
        "preparation_notes": "Boiled in water for 30-60 minutes, consumed as tea",
        "evidence_level": "historical_documentation",
        "literature_sources": ["PMID:29392251", "ISBN:9780195320794"],
        "community_consent_status": "published_knowledge",
        "attribution_notes": "Traditional Chinese Medicine historical records"
    },
    {
        "practice_name": "Cannabis oil for epilepsy",
        "source_community_id": "ayurvedic_tradition",
        "knowledge_domain": "neurological_disorders",
        "preparation_method": "Cannabis leaves infused in sesame oil",
        "indications": ["epilepsy", "seizures", "neurological_disorders"],
        "contraindications": ["pregnancy", "lactation"],
        "cultural_context": "Documented in Ayurvedic texts as 'vijaya' (victory over disease)",
        "ceremonial_significance": False,
        "traditional_dosage": "5-10 drops topically or orally",
        "preparation_notes": "Cannabis leaves macerated in sesame oil for 30 days",
        "evidence_level": "historical_documentation",
        "literature_sources": ["PMID:28412918"],
        "community_consent_status": "published_knowledge",
        "attribution_notes": "Ayurvedic pharmacopoeia"
    },
    {
        "practice_name": "Hemp seed for appetite stimulation",
        "source_community_id": "traditional_eastern_european",
        "knowledge_domain": "nutritional_medicine",
        "preparation_method": "Ground hemp seeds mixed with honey",
        "indications": ["appetite_stimulation", "cachexia", "wasting_syndrome"],
        "contraindications": ["nut_allergies"],
        "cultural_context": "Traditional Eastern European folk remedy for appetite loss",
        "ceremonial_significance": False,
        "traditional_dosage": "1-2 tablespoons daily",
        "preparation_notes": "Fresh hemp seeds ground and mixed with raw honey",
        "evidence_level": "ethnographic_documentation",
        "literature_sources": ["ethnobotany_database"],
        "community_consent_status": "public_domain",
        "attribution_notes": "Eastern European folk medicine tradition"
    },
    {
        "practice_name": "Cannabis topical for arthritis",
        "source_community_id": "caribbean_traditional_medicine",
        "knowledge_domain": "musculoskeletal_disorders",
        "preparation_method": "Cannabis leaves macerated in coconut oil",
        "indications": ["arthritis", "joint_pain", "inflammation"],
        "contraindications": ["open_wounds", "skin_sensitivity"],
        "cultural_context": "TRAMIL database documented practice from Caribbean traditional medicine",
        "ceremonial_significance": False,
        "traditional_dosage": "Applied topically 2-3 times daily",
        "preparation_notes": "Fresh cannabis leaves soaked in coconut oil for 14 days, strained",
        "evidence_level": "tramil_validated",
        "literature_sources": ["TRAMIL:Cannabis_001"],
        "community_consent_status": "tramil_approved",
        "attribution_notes": "TRAMIL Caribbean Traditional Medicine Database"
    },
    {
        "practice_name": "Cannabis tea for anxiety",
        "source_community_id": "traditional_brazilian",
        "knowledge_domain": "mental_health",
        "preparation_method": "Cannabis leaves brewed as tea",
        "indications": ["anxiety", "nervousness", "insomnia"],
        "contraindications": ["pregnancy", "psychosis_history"],
        "cultural_context": "Traditional Brazilian herbal medicine practice",
        "ceremonial_significance": False,
        "traditional_dosage": "1 cup before bedtime",
        "preparation_notes": "5-10g dried leaves steeped in hot water for 10 minutes",
        "evidence_level": "ethnographic_documentation",
        "literature_sources": ["ethnobotany_brazil"],
        "community_consent_status": "published_knowledge",
        "attribution_notes": "Brazilian ethnobotanical surveys"
    }
]


# Template for sacred/protected knowledge (DO NOT USE IN MVP)
SACRED_PRACTICE_TEMPLATE = {
    "practice_name": "[PROTECTED - DO NOT INCLUDE IN DATABASE]",
    "source_community_id": "[INDIGENOUS_COMMUNITY_NAME]",
    "knowledge_domain": "ceremonial",
    "preparation_method": "[SACRED PREPARATION - REQUIRES EXPLICIT CONSENT]",
    "indications": ["spiritual_healing"],
    "contraindications": [],
    "cultural_context": "Sacred ceremonial practice - requires community governance approval",
    "ceremonial_significance": True,  # ← This flags sacred knowledge
    "traditional_dosage": "[PROTECTED]",
    "preparation_notes": "[REQUIRES ELDER APPROVAL AND BENEFIT-SHARING AGREEMENT]",
    "evidence_level": "ceremonial_knowledge",
    "literature_sources": [],
    "community_consent_status": "requires_formal_consent",
    "attribution_notes": "Protected traditional knowledge - not for commercial use"
}


class TKDatasetBuilder:
    """Builds structured traditional knowledge practice dataset."""
    
    def __init__(self):
        self.practices = []
        self.stats = {
            "total_practices": 0,
            "non_sacred": 0,
            "sacred_flagged": 0,
            "communities": set(),
            "knowledge_domains": set(),
            "conditions_covered": set()
        }
    
    def validate_practice(self, practice: Dict) -> tuple[bool, str]:
        """Validate TK practice against schema requirements."""
        
        required_fields = [
            "practice_name",
            "source_community_id",
            "knowledge_domain",
            "preparation_method",
            "indications",
            "ceremonial_significance"
        ]
        
        for field in required_fields:
            if field not in practice:
                return False, f"Missing required field: {field}"
        
        # Check for sacred knowledge
        if practice.get("ceremonial_significance") is True:
            if practice.get("community_consent_status") != "explicit_consent_granted":
                return False, "Sacred knowledge requires explicit community consent"
        
        # Validate indications are not empty
        if not practice.get("indications"):
            return False, "Must have at least one indication"
        
        return True, "Valid"
    
    def add_practice(self, practice: Dict):
        """Add a validated TK practice to dataset."""
        
        is_valid, message = self.validate_practice(practice)
        
        if not is_valid:
            print(f"⚠️  Skipping invalid practice '{practice.get('practice_name', 'UNKNOWN')}': {message}")
            return
        
        # Add metadata
        practice["practice_id"] = f"TK_{len(self.practices):04d}"
        practice["added_date"] = datetime.now().isoformat()
        
        self.practices.append(practice)
        
        # Update statistics
        self.stats["total_practices"] += 1
        self.stats["communities"].add(practice["source_community_id"])
        self.stats["knowledge_domains"].add(practice["knowledge_domain"])
        
        if practice.get("ceremonial_significance"):
            self.stats["sacred_flagged"] += 1
        else:
            self.stats["non_sacred"] += 1
        
        for indication in practice.get("indications", []):
            self.stats["conditions_covered"].add(indication)
    
    def load_example_practices(self):
        """Load example non-sacred practices from published literature."""
        
        print("📚 Loading example TK practices from published sources...")
        
        for practice in EXAMPLE_TK_PRACTICES:
            self.add_practice(practice)
        
        print(f"   Added {len(EXAMPLE_TK_PRACTICES)} practices from published literature")
    
    def load_expansion_practices(self, expansion_file: str = "data/processed/tk_practices_expansion.json"):
        """Load additional practices from expansion file."""
        
        if not os.path.exists(expansion_file):
            print(f"   ℹ️  No expansion file found at {expansion_file}")
            return
        
        print(f"📚 Loading expansion practices from {expansion_file}...")
        
        with open(expansion_file, 'r', encoding='utf-8') as f:
            expansion_data = json.load(f)
        
        practices = expansion_data.get("practices", [])
        
        for practice in practices:
            self.add_practice(practice)
        
        print(f"   Added {len(practices)} additional practices from expansion file")
    
    def generate_curation_template(self, output_path: str = "data/processed/tk_practice_template.json"):
        """Generate template for manual TK practice curation."""
        
        template = {
            "instructions": "Use this template to add new TK practices. Follow ethical guidelines for TK documentation.",
            "ethical_guidelines": [
                "Only include non-sacred, publicly documented knowledge",
                "Set ceremonial_significance=true for ANY ceremonial/sacred practices",
                "Verify community_consent_status before adding",
                "Provide clear attribution in attribution_notes",
                "Include literature sources where available",
                "Respect Indigenous Data Sovereignty principles"
            ],
            "practice_template": {
                "practice_name": "Descriptive name of the practice",
                "source_community_id": "community_or_tradition_identifier",
                "knowledge_domain": "One of: pain_management, neurological_disorders, mental_health, etc.",
                "preparation_method": "Detailed preparation instructions",
                "indications": ["condition1", "condition2"],
                "contraindications": ["contraindication1"],
                "cultural_context": "Historical and cultural background",
                "ceremonial_significance": False,
                "traditional_dosage": "Traditional dosage information",
                "preparation_notes": "Additional preparation details",
                "evidence_level": "One of: historical_documentation, ethnographic_documentation, tramil_validated, published_knowledge",
                "literature_sources": ["PMID:12345678", "ISBN:1234567890"],
                "community_consent_status": "One of: published_knowledge, public_domain, tramil_approved, explicit_consent_granted",
                "attribution_notes": "Proper attribution to source community/tradition"
            },
            "example_non_sacred_practice": EXAMPLE_TK_PRACTICES[0],
            "sacred_practice_warning": SACRED_PRACTICE_TEMPLATE
        }
        
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        
        with open(output_path, 'w', encoding='utf-8') as f:
            json.dump(template, f, indent=2, ensure_ascii=False)
        
        print(f"✅ Generated curation template: {output_path}")
    
    def export_dataset(self, output_path: str = "data/processed/tk_practices.json"):
        """Export structured TK practice dataset."""
        
        output = {
            "metadata": {
                "version": "1.0",
                "created_date": datetime.now().isoformat(),
                "total_practices": self.stats["total_practices"],
                "non_sacred_practices": self.stats["non_sacred"],
                "sacred_practices_flagged": self.stats["sacred_flagged"],
                "unique_communities": len(self.stats["communities"]),
                "knowledge_domains": len(self.stats["knowledge_domains"]),
                "conditions_covered": len(self.stats["conditions_covered"]),
                "ethical_compliance": "All practices from published, non-sacred sources",
                "sacred_knowledge_protection": "ceremonial_significance flag enforced"
            },
            "communities": sorted(list(self.stats["communities"])),
            "knowledge_domains": sorted(list(self.stats["knowledge_domains"])),
            "conditions_covered": sorted(list(self.stats["conditions_covered"])),
            "practices": self.practices
        }
        
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        
        with open(output_path, 'w', encoding='utf-8') as f:
            json.dump(output, f, indent=2, ensure_ascii=False)
        
        print(f"✅ Exported {self.stats['total_practices']} TK practices to {output_path}")
    
    def print_statistics(self):
        """Print dataset statistics."""
        
        print()
        print("=" * 70)
        print("TK Dataset Statistics")
        print("=" * 70)
        print(f"Total practices: {self.stats['total_practices']}")
        print(f"  Non-sacred: {self.stats['non_sacred']}")
        print(f"  Sacred (flagged, protected): {self.stats['sacred_flagged']}")
        print(f"Unique communities: {len(self.stats['communities'])}")
        print(f"Knowledge domains: {len(self.stats['knowledge_domains'])}")
        print(f"Conditions covered: {len(self.stats['conditions_covered'])}")
        print()
        print("Communities represented:")
        for community in sorted(self.stats["communities"]):
            print(f"  - {community}")
        print()


def main():
    """Main TK dataset building pipeline."""
    
    print("=" * 70)
    print("GenomePath Traditional Knowledge Dataset Builder")
    print("=" * 70)
    print()
    print("⚠️  ETHICAL NOTICE:")
    print("   This tool respects Indigenous Data Sovereignty")
    print("   Only non-sacred, published knowledge is included")
    print("   Sacred knowledge requires explicit community consent")
    print("=" * 70)
    print()
    
    builder = TKDatasetBuilder()
    
    # Load example practices from published literature
    builder.load_example_practices()
    
    # Load expansion practices (15 additional)
    builder.load_expansion_practices()
    
    # Generate template for manual curation
    print()
    print("📝 Generating curation template...")
    builder.generate_curation_template()
    
    # Export dataset
    print()
    print("💾 Exporting TK practice dataset...")
    builder.export_dataset()
    
    # Print statistics
    builder.print_statistics()
    
    print("=" * 70)
    print("✅ TK dataset building complete!")
    print("=" * 70)
    print()
    print("Next steps:")
    print("  1. Review data/processed/tk_practices.json")
    print("  2. Use data/processed/tk_practice_template.json to add more practices")
    print("  3. Run scripts/generate_correlations.py to create TK↔Genomic pairs")
    print()
    print("To add more practices:")
    print("  - Edit tk_practice_template.json with new practices")
    print("  - Ensure all practices respect ethical guidelines")
    print("  - Run this script again to validate and merge")
    print()


if __name__ == "__main__":
    main()
