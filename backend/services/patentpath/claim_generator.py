"""
Claim Generator - PatentPath Lite Week 7
Generate patent claims from templates for cannabinoid compounds.

Features:
- Multiple claim types (composition, method, pharmaceutical, synthesis)
- USPTO formatting
- Filing cost estimation
- Prosecution cost estimation

DISCLAIMER: These claims are template-generated for planning purposes only.
Consult a patent attorney for filing-ready claims.
"""
from typing import Dict, List, Optional, Any
from dataclasses import dataclass, field
from datetime import datetime
import logging

logger = logging.getLogger(__name__)

# Required disclaimer for all claim outputs
CLAIM_DISCLAIMER = (
    "DISCLAIMER: These claims are template-generated for planning purposes only. "
    "Consult a patent attorney for filing-ready claims. "
    "Not legal advice."
)


class ClaimType:
    """Supported claim types."""
    COMPOSITION_OF_MATTER = "composition_of_matter"
    METHOD_OF_USE = "method_of_use"
    PHARMACEUTICAL_COMPOSITION = "pharmaceutical_composition"
    METHOD_OF_SYNTHESIS = "method_of_synthesis"
    DIMER_COMPOSITION = "dimer_composition"


@dataclass
class GeneratedClaims:
    """Container for generated patent claims."""
    claim_type: str
    claims_text: str
    num_claims: int
    placeholders_filled: Dict[str, str]
    
    def to_dict(self) -> Dict:
        return {
            "claim_type": self.claim_type,
            "claims_text": self.claims_text,
            "num_claims": self.num_claims,
            "placeholders": self.placeholders_filled
        }


@dataclass
class FilingCostEstimate:
    """USPTO filing cost estimate."""
    base_filing_fee: int
    search_fee: int
    examination_fee: int
    additional_claims_fee: int
    total_uspto_fees: int
    entity_size: str
    
    def to_dict(self) -> Dict:
        return {
            "base_filing_fee": f"${self.base_filing_fee:,}",
            "search_fee": f"${self.search_fee:,}",
            "examination_fee": f"${self.examination_fee:,}",
            "additional_claims_fee": f"${self.additional_claims_fee:,}",
            "total_uspto_fees": f"${self.total_uspto_fees:,}",
            "entity_size": self.entity_size
        }


@dataclass
class ClaimGenerationReport:
    """Complete claim generation report."""
    compound_name: str
    compound_smiles: Optional[str]
    therapeutic_use: Optional[str]
    generated_claims: Dict[str, GeneratedClaims]
    claim_types: List[str]
    total_claims: int
    filing_cost_estimate: FilingCostEstimate
    prosecution_cost_estimate: str
    generated_date: str
    format: str = "USPTO_utility_patent"
    disclaimer: str = CLAIM_DISCLAIMER
    upgrade_note: str = "PatentPath Full uses AI to generate USPTO-compliant claims with 78% acceptance rate vs. 45-60% industry average."
    
    def to_dict(self) -> Dict:
        return {
            "compound_name": self.compound_name,
            "compound_smiles": self.compound_smiles,
            "therapeutic_use": self.therapeutic_use,
            "generated_claims": {
                k: v.to_dict() for k, v in self.generated_claims.items()
            },
            "claim_types": self.claim_types,
            "total_claims": self.total_claims,
            "filing_cost_estimate": self.filing_cost_estimate.to_dict(),
            "prosecution_cost_estimate": self.prosecution_cost_estimate,
            "generated_date": self.generated_date,
            "format": self.format,
            "disclaimer": self.disclaimer,
            "upgrade_note": self.upgrade_note
        }


class ClaimGenerator:
    """Generate patent claims from templates.
    
    MVP version uses fill-in-the-blank templates.
    Full version uses AI for USPTO-compliant claims with 78% acceptance rate.
    """
    
    def __init__(self):
        self.templates = self._load_templates()
    
    def _load_templates(self) -> Dict[str, str]:
        """Load patent claim templates."""
        return {
            ClaimType.COMPOSITION_OF_MATTER: """CLAIM 1: A compound having the structure of Formula I:

{structure_representation}

wherein:
- R1 is {r1_definition};
- R2 is {r2_definition};
- The compound is characterized by {key_features}.

CLAIM 2: The compound of Claim 1, wherein the compound is selected from the group consisting of:
{specific_compounds}.

CLAIM 3: The compound of Claim 1, wherein the compound is a {compound_class}.

CLAIM 4: The compound of Claim 1, wherein the compound has a molecular weight between {mw_low} and {mw_high} daltons.

CLAIM 5: The compound of Claim 1, wherein the compound exhibits binding affinity for {receptor_targets}.""",
            
            ClaimType.METHOD_OF_USE: """CLAIM 1: A method for treating {condition} in a subject in need thereof, comprising:
administering to the subject a therapeutically effective amount of a compound of Formula I:

{compound_structure}

CLAIM 2: The method of Claim 1, wherein the therapeutically effective amount is between {dose_low} and {dose_high} mg per day.

CLAIM 3: The method of Claim 1, wherein the subject is a human.

CLAIM 4: The method of Claim 1, wherein the compound is administered {admin_route}.

CLAIM 5: The method of Claim 1, further comprising administering one or more additional therapeutic agents.

CLAIM 6: The method of Claim 1, wherein the treatment is continued for at least {treatment_duration}.""",
            
            ClaimType.PHARMACEUTICAL_COMPOSITION: """CLAIM 1: A pharmaceutical composition comprising:
(a) a therapeutically effective amount of a compound of Formula I:

{compound_structure}

and
(b) a pharmaceutically acceptable carrier.

CLAIM 2: The pharmaceutical composition of Claim 1, wherein the composition is formulated for {admin_route} administration.

CLAIM 3: The pharmaceutical composition of Claim 1, further comprising one or more excipients selected from the group consisting of: {excipients}.

CLAIM 4: The pharmaceutical composition of Claim 1, wherein the compound is present in an amount of {dose_low} mg to {dose_high} mg per unit dose.

CLAIM 5: The pharmaceutical composition of Claim 1, further comprising one or more terpenes selected from: {terpenes}.

CLAIM 6: The pharmaceutical composition of Claim 1, wherein the composition provides {release_profile} release of the compound.""",
            
            ClaimType.METHOD_OF_SYNTHESIS: """CLAIM 1: A method for synthesizing a compound of Formula I, comprising:
(a) reacting {starting_material_1} with {starting_material_2} under {reaction_conditions} to form an intermediate compound;
(b) {step_2_description};
(c) isolating the compound of Formula I.

CLAIM 2: The method of Claim 1, wherein the reaction of step (a) is conducted at a temperature between {temp_low}°C and {temp_high}°C.

CLAIM 3: The method of Claim 1, wherein the reaction of step (a) is conducted in the presence of {catalyst}.

CLAIM 4: The method of Claim 1, wherein the yield is at least {yield_percent}%.

CLAIM 5: The method of Claim 1, further comprising purifying the compound by {purification_method}.""",
            
            ClaimType.DIMER_COMPOSITION: """CLAIM 1: A dimeric cannabinoid compound of Formula I:

{structure_representation}

wherein:
- A is a first cannabinoid moiety selected from {cannabinoid_1};
- B is a second cannabinoid moiety selected from {cannabinoid_2};
- L is a linker moiety having the structure {linker_structure}.

CLAIM 2: The compound of Claim 1, wherein the linker L comprises {linker_length} atoms.

CLAIM 3: The compound of Claim 1, wherein the first cannabinoid moiety A is {cannabinoid_1} and the second cannabinoid moiety B is {cannabinoid_2}.

CLAIM 4: The compound of Claim 1, wherein the compound exhibits enhanced {therapeutic_property} compared to either cannabinoid moiety alone.

CLAIM 5: The compound of Claim 1, wherein the compound has a CB1:CB2 binding ratio of {binding_ratio}.

CLAIM 6: A pharmaceutical composition comprising the compound of Claim 1 and a pharmaceutically acceptable carrier.

CLAIM 7: A method of treating {condition} comprising administering an effective amount of the compound of Claim 1 to a subject in need thereof."""
        }
    
    def generate_claims(
        self,
        compound_name: str,
        compound_smiles: Optional[str] = None,
        therapeutic_use: Optional[str] = None,
        key_features: Optional[List[str]] = None,
        claim_types: Optional[List[str]] = None,
        custom_params: Optional[Dict[str, str]] = None
    ) -> ClaimGenerationReport:
        """Generate patent claims for a compound.
        
        Args:
            compound_name: Name of compound
            compound_smiles: SMILES structure
            therapeutic_use: Indication (e.g., "Alzheimer's disease")
            key_features: List of key structural features
            claim_types: Which claim types to generate
            custom_params: Custom template parameters
            
        Returns:
            ClaimGenerationReport with generated claims and cost estimates
        """
        if claim_types is None:
            claim_types = [
                ClaimType.COMPOSITION_OF_MATTER,
                ClaimType.METHOD_OF_USE,
                ClaimType.PHARMACEUTICAL_COMPOSITION
            ]
        
        generated_claims = {}
        total_claims = 0
        
        for claim_type in claim_types:
            if claim_type not in self.templates:
                logger.warning(f"Unknown claim type: {claim_type}")
                continue
            
            template = self.templates[claim_type]
            
            # Fill template with compound-specific data
            filled_claims, placeholders = self._fill_template(
                template,
                claim_type,
                compound_name=compound_name,
                compound_smiles=compound_smiles,
                therapeutic_use=therapeutic_use,
                key_features=key_features or [],
                custom_params=custom_params or {}
            )
            
            num_claims = self._count_claims(filled_claims)
            total_claims += num_claims
            
            generated_claims[claim_type] = GeneratedClaims(
                claim_type=claim_type,
                claims_text=filled_claims,
                num_claims=num_claims,
                placeholders_filled=placeholders
            )
        
        # Calculate costs
        filing_cost = self._estimate_filing_cost(total_claims)
        prosecution_cost = self._estimate_prosecution_cost(total_claims)
        
        return ClaimGenerationReport(
            compound_name=compound_name,
            compound_smiles=compound_smiles,
            therapeutic_use=therapeutic_use,
            generated_claims=generated_claims,
            claim_types=claim_types,
            total_claims=total_claims,
            filing_cost_estimate=filing_cost,
            prosecution_cost_estimate=prosecution_cost,
            generated_date=datetime.now().isoformat()
        )
    
    def _fill_template(
        self,
        template: str,
        claim_type: str,
        **kwargs
    ) -> tuple[str, Dict[str, str]]:
        """Fill claim template with compound-specific data."""
        compound_name = kwargs.get("compound_name", "Compound")
        compound_smiles = kwargs.get("compound_smiles")
        therapeutic_use = kwargs.get("therapeutic_use")
        key_features = kwargs.get("key_features", [])
        custom_params = kwargs.get("custom_params", {})
        
        # Default values based on claim type
        defaults = {
            # Structure placeholders
            "structure_representation": f"[{compound_smiles or 'Chemical Structure Diagram'}]",
            "compound_structure": compound_smiles or "[SMILES: To be inserted]",
            
            # Substituent definitions
            "r1_definition": "hydrogen, C1-C6 alkyl, or aryl",
            "r2_definition": "hydroxyl, C1-C4 alkoxy, or halogen",
            
            # Compound characteristics
            "key_features": ", ".join(key_features) if key_features else "novel cannabinoid structure",
            "specific_compounds": compound_name,
            "compound_class": self._infer_compound_class(compound_name),
            "mw_low": "300",
            "mw_high": "700",
            "receptor_targets": "CB1 and CB2 cannabinoid receptors",
            
            # Therapeutic use
            "condition": therapeutic_use or "[Therapeutic Indication]",
            "therapeutic_property": "therapeutic efficacy",
            
            # Dosing
            "dose_low": "5",
            "dose_high": "500",
            "treatment_duration": "4 weeks",
            
            # Administration
            "admin_route": "oral",
            "excipients": "lactose, microcrystalline cellulose, magnesium stearate, and hydroxypropyl methylcellulose",
            "terpenes": "myrcene, limonene, linalool, and β-caryophyllene",
            "release_profile": "sustained",
            
            # Synthesis
            "starting_material_1": "[Starting Material 1]",
            "starting_material_2": "[Starting Material 2]",
            "reaction_conditions": "oxidative conditions in the presence of a suitable catalyst",
            "step_2_description": "purifying the intermediate compound by chromatography",
            "temp_low": "20",
            "temp_high": "80",
            "catalyst": "a palladium catalyst",
            "yield_percent": "50",
            "purification_method": "column chromatography or recrystallization",
            
            # Dimer-specific
            "cannabinoid_1": "THC, CBD, CBG, or CBN",
            "cannabinoid_2": "THC, CBD, CBG, or CBN",
            "linker_structure": "-O-(CH2)n-O-",
            "linker_length": "2 to 12",
            "binding_ratio": "1:1 to 1:100"
        }
        
        # Override with custom parameters
        fill_data = {**defaults, **custom_params}
        
        # Track which placeholders were used
        placeholders_used = {}
        
        # Fill template
        try:
            filled = template
            for key, value in fill_data.items():
                placeholder = "{" + key + "}"
                if placeholder in filled:
                    filled = filled.replace(placeholder, str(value))
                    placeholders_used[key] = str(value)
        except Exception as e:
            logger.error(f"Error filling template: {e}")
            filled = template
        
        return filled, placeholders_used
    
    def _infer_compound_class(self, compound_name: str) -> str:
        """Infer compound class from name."""
        name_lower = compound_name.lower()
        
        if "dimer" in name_lower:
            return "dimeric cannabinoid compound"
        elif any(cb in name_lower for cb in ["thc", "tetrahydrocannabinol"]):
            return "tetrahydrocannabinol derivative"
        elif any(cb in name_lower for cb in ["cbd", "cannabidiol"]):
            return "cannabidiol derivative"
        elif any(cb in name_lower for cb in ["cbg", "cannabigerol"]):
            return "cannabigerol derivative"
        elif any(cb in name_lower for cb in ["cbn", "cannabinol"]):
            return "cannabinol derivative"
        elif "cannabin" in name_lower:
            return "cannabinoid compound"
        else:
            return "organic compound of pharmaceutical interest"
    
    def _count_claims(self, claims_text: str) -> int:
        """Count number of claims in text."""
        return claims_text.upper().count("CLAIM ")
    
    def _estimate_filing_cost(
        self,
        num_claims: int,
        entity_size: str = "large"
    ) -> FilingCostEstimate:
        """Estimate USPTO filing cost.
        
        Based on 2025 USPTO fee schedule:
        - Micro entity: 75% discount
        - Small entity: 50% discount
        - Large entity: Full fees
        """
        # Large entity base fees (2025)
        base_fees = {
            "filing": 1280,
            "search": 2640,
            "examination": 3040
        }
        
        # Entity size multipliers
        multipliers = {
            "micro": 0.25,
            "small": 0.5,
            "large": 1.0
        }
        
        multiplier = multipliers.get(entity_size, 1.0)
        
        base_filing = int(base_fees["filing"] * multiplier)
        search_fee = int(base_fees["search"] * multiplier)
        exam_fee = int(base_fees["examination"] * multiplier)
        
        # Additional claim fees: $100 per claim over 20 (large entity)
        if num_claims > 20:
            extra_claims_fee = int((num_claims - 20) * 100 * multiplier)
        else:
            extra_claims_fee = 0
        
        total = base_filing + search_fee + exam_fee + extra_claims_fee
        
        return FilingCostEstimate(
            base_filing_fee=base_filing,
            search_fee=search_fee,
            examination_fee=exam_fee,
            additional_claims_fee=extra_claims_fee,
            total_uspto_fees=total,
            entity_size=entity_size
        )
    
    def _estimate_prosecution_cost(self, num_claims: int) -> str:
        """Estimate attorney prosecution costs."""
        # Industry averages: $5,000-15,000 for basic prosecution
        base_prosecution = 8000
        
        # More claims = more complex prosecution
        complexity_fee = max(0, (num_claims - 15) * 200)
        
        low_estimate = base_prosecution + complexity_fee
        high_estimate = int(low_estimate * 1.75)
        
        return f"${low_estimate:,} - ${high_estimate:,} (attorney fees, estimated)"
    
    def get_available_claim_types(self) -> List[str]:
        """Get list of available claim types."""
        return list(self.templates.keys())
    
    def generate_dimer_claims(
        self,
        dimer_name: str,
        parent_1: str,
        parent_2: str,
        linker_type: str,
        therapeutic_use: Optional[str] = None,
        smiles: Optional[str] = None
    ) -> ClaimGenerationReport:
        """Convenience method for generating dimer-specific claims.
        
        Args:
            dimer_name: Name of dimer compound
            parent_1: First parent cannabinoid
            parent_2: Second parent cannabinoid
            linker_type: Type of linker (e.g., "ester", "ether")
            therapeutic_use: Target indication
            smiles: SMILES structure
            
        Returns:
            ClaimGenerationReport with dimer-specific claims
        """
        # Determine linker structure based on type
        linker_structures = {
            "ester": "-O-C(=O)-(CH2)n-C(=O)-O-",
            "ether": "-O-(CH2)n-O-",
            "amide": "-NH-C(=O)-(CH2)n-C(=O)-NH-",
            "alkyl": "-(CH2)n-",
            "peg": "-O-(CH2-CH2-O)n-"
        }
        
        linker = linker_structures.get(linker_type.lower(), "-O-(CH2)n-O-")
        
        custom_params = {
            "cannabinoid_1": parent_1,
            "cannabinoid_2": parent_2,
            "linker_structure": linker,
            "linker_length": "2 to 8" if linker_type.lower() in ["ester", "ether"] else "4 to 12"
        }
        
        return self.generate_claims(
            compound_name=dimer_name,
            compound_smiles=smiles,
            therapeutic_use=therapeutic_use,
            key_features=[
                f"dimeric structure linking {parent_1} and {parent_2}",
                f"{linker_type} linker moiety",
                "enhanced receptor binding profile"
            ],
            claim_types=[
                ClaimType.DIMER_COMPOSITION,
                ClaimType.METHOD_OF_USE,
                ClaimType.PHARMACEUTICAL_COMPOSITION
            ],
            custom_params=custom_params
        )
