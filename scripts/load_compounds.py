"""
Load 63 Cannabinoid Compounds from JSON to Database
Populates the compound library for FDA CMC documentation

Source: data/training/neurobotanica_complete_dataset_63compounds.json
"""
import json
import sys
from pathlib import Path

# Add project root to path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

from sqlalchemy.orm import Session
from backend.models.database import engine, SessionLocal, Base
from backend.models.compound import Cannabinoid, DimericPrediction


def determine_compound_class(compound_data: dict) -> str:
    """Determine compound classification."""
    name = compound_data.get("name", "").lower()
    
    # Known synthetic cannabinoids
    synthetics = ["nabilone", "dronabinol", "win-55212", "cp55940", "jwh"]
    if any(s in name for s in synthetics):
        return "Synthetic"
    
    # Endocannabinoids
    endos = ["anandamide", "2-ag", "arachidonoyl"]
    if any(e in name for e in endos):
        return "Endocannabinoid"
    
    # Default to phytocannabinoid
    return "Phytocannabinoid"


def check_fda_approved(compound_data: dict) -> tuple:
    """Check if compound is FDA-approved and return drug name."""
    name = compound_data.get("name", "").lower()
    abbrev = compound_data.get("abbreviation", "").lower()
    
    fda_drugs = {
        "cbd": ("Epidiolex", ["Dravet syndrome", "Lennox-Gastaut syndrome", "Tuberous sclerosis"], "V"),
        "cannabidiol": ("Epidiolex", ["Dravet syndrome", "Lennox-Gastaut syndrome", "Tuberous sclerosis"], "V"),
        "dronabinol": ("Marinol", ["CINV", "AIDS wasting"], "III"),
        "thc": (None, None, None),  # THC itself not approved, only dronabinol
        "nabilone": ("Cesamet", ["CINV"], "II"),
    }
    
    for key, (drug, indications, schedule) in fda_drugs.items():
        if key in name or key in abbrev:
            if drug:
                return True, drug, indications, schedule
    
    return False, None, None, None


def load_compounds(json_path: str, db: Session) -> int:
    """Load compounds from JSON file into database."""
    print(f"📂 Loading compounds from: {json_path}")
    
    with open(json_path, 'r', encoding='utf-8') as f:
        data = json.load(f)
    
    # Handle different JSON structures
    if isinstance(data, list):
        compounds = data
    elif "compounds" in data:
        compounds = data["compounds"]
    else:
        compounds = [data]
    
    print(f"🧬 Found {len(compounds)} compounds in JSON")
    
    loaded_count = 0
    skipped_count = 0
    seen_names = set()  # Track names within this batch
    
    for compound_data in compounds:
        # Support both 'name' and 'compound_name' field names
        name = compound_data.get("compound_name") or compound_data.get("name")
        
        if not name:
            continue
        
        # Skip if already seen in this batch
        if name in seen_names:
            skipped_count += 1
            continue
        seen_names.add(name)
        
        # Check if already exists
        existing = db.query(Cannabinoid).filter(
            Cannabinoid.name == name
        ).first()
        
        if existing:
            skipped_count += 1
            continue
        
        # Check FDA status
        is_approved, drug_name, indications, schedule = check_fda_approved(compound_data)
        
        # Extract receptor affinities
        affinities = compound_data.get("receptor_affinities", {})
        
        # Extract RDKit descriptors from nested structure
        rdkit = compound_data.get("rdkit_descriptors", {})
        
        # Create compound record
        compound = Cannabinoid(
            name=name,
            abbreviation=compound_data.get("abbreviation") or name[:5].upper(),
            smiles=compound_data.get("smiles", ""),
            inchi=compound_data.get("inchi"),
            inchi_key=compound_data.get("inchi_key"),
            
            # Classification
            compound_class=determine_compound_class(compound_data),
            is_natural=compound_data.get("is_natural", True),
            is_metabolite=compound_data.get("is_metabolite", False),
            is_dimeric=compound_data.get("is_dimeric", False),
            
            # 2D Descriptors - check both top-level and rdkit_descriptors
            molecular_weight=compound_data.get("molecular_weight") or rdkit.get("molecular_weight"),
            exact_mass=compound_data.get("exact_mass"),
            logp=compound_data.get("logp") or compound_data.get("logP") or rdkit.get("logP"),
            tpsa=compound_data.get("tpsa") or compound_data.get("TPSA") or rdkit.get("tpsa"),
            h_bond_donors=compound_data.get("h_bond_donors") or compound_data.get("HBD") or rdkit.get("h_bond_donors"),
            h_bond_acceptors=compound_data.get("h_bond_acceptors") or compound_data.get("HBA") or rdkit.get("h_bond_acceptors"),
            rotatable_bonds=compound_data.get("rotatable_bonds") or rdkit.get("rotatable_bonds"),
            aromatic_rings=compound_data.get("aromatic_rings") or rdkit.get("num_aromatic_rings"),
            heavy_atoms=compound_data.get("heavy_atoms") or rdkit.get("num_heavy_atoms"),
            fraction_csp3=compound_data.get("fraction_csp3") or rdkit.get("sp3_fraction"),
            rdkit_descriptors=rdkit if rdkit else compound_data.get("descriptors"),
            
            # Receptor affinities
            cb1_affinity_ki=affinities.get("CB1", {}).get("ki") if isinstance(affinities.get("CB1"), dict) else affinities.get("CB1"),
            cb1_source=affinities.get("CB1", {}).get("source") if isinstance(affinities.get("CB1"), dict) else None,
            cb2_affinity_ki=affinities.get("CB2", {}).get("ki") if isinstance(affinities.get("CB2"), dict) else affinities.get("CB2"),
            cb2_source=affinities.get("CB2", {}).get("source") if isinstance(affinities.get("CB2"), dict) else None,
            trpv1_affinity=affinities.get("TRPV1"),
            gpr55_affinity=affinities.get("GPR55"),
            receptor_affinities=affinities,
            
            # Pharmacology
            mechanism_of_action=compound_data.get("mechanism_of_action"),
            therapeutic_categories=compound_data.get("therapeutic_categories"),
            conditions_treated=compound_data.get("conditions_treated"),
            
            # FDA Status
            fda_approved=is_approved,
            fda_drug_name=drug_name,
            fda_indications=indications,
            schedule_classification=schedule,
            
            # External IDs
            chembl_id=compound_data.get("chembl_id"),
            pubchem_cid=compound_data.get("pubchem_cid"),
            drugbank_id=compound_data.get("drugbank_id"),
            
            # Dimeric predictions
            dimeric_potential_score=compound_data.get("dimeric_potential"),
            predicted_synergy_score=compound_data.get("synergy_score"),
            dimeric_predictions=compound_data.get("dimeric_predictions")
        )
        
        db.add(compound)
        loaded_count += 1
    
    db.commit()
    print(f"✅ Loaded {loaded_count} new compounds (skipped {skipped_count} existing)")
    
    return loaded_count


def load_dimeric_predictions(json_path: str, db: Session) -> int:
    """Load dimeric predictions from JSON file."""
    print(f"📂 Loading dimeric predictions from: {json_path}")
    
    with open(json_path, 'r', encoding='utf-8') as f:
        data = json.load(f)
    
    predictions = data if isinstance(data, list) else data.get("predictions", [])
    print(f"🔮 Found {len(predictions)} dimeric predictions")
    
    loaded_count = 0
    
    for pred_data in predictions:
        dimer_name = pred_data.get("dimer_name") or pred_data.get("name")
        
        if not dimer_name:
            continue
        
        # Check if exists
        existing = db.query(DimericPrediction).filter(
            DimericPrediction.dimer_name == dimer_name
        ).first()
        
        if existing:
            continue
        
        # Get therapeutic potential as JSON string if needed
        therapeutic = pred_data.get("therapeutic_potential")
        if isinstance(therapeutic, dict):
            therapeutic = json.dumps(therapeutic)
        
        # Get supporting studies as JSON string
        supporting = pred_data.get("supporting_studies")
        if isinstance(supporting, (list, dict)):
            supporting = json.dumps(supporting)
        
        prediction = DimericPrediction(
            dimer_name=dimer_name,
            # Support both 'component_1/2' and 'parent_1/2' field names
            parent_1_name=pred_data.get("component_1") or pred_data.get("parent_1") or pred_data.get("parent1") or "Unknown",
            parent_2_name=pred_data.get("component_2") or pred_data.get("parent_2") or pred_data.get("parent2") or "Unknown",
            predicted_smiles=pred_data.get("smiles"),
            
            triangulation_score=pred_data.get("triangulation_score") or pred_data.get("formation_probability"),
            synergy_prediction=pred_data.get("synergy_prediction") or pred_data.get("synergy_score"),
            novelty_score=pred_data.get("novelty_score"),
            therapeutic_potential=therapeutic,
            
            predicted_cb1_affinity=pred_data.get("predicted_cb1"),
            predicted_cb2_affinity=pred_data.get("predicted_cb2"),
            predicted_selectivity=pred_data.get("selectivity"),
            
            supporting_studies=supporting,
            confidence_score=pred_data.get("confidence"),
            
            prior_art_score=pred_data.get("prior_art_score"),
            novelty_assessment=pred_data.get("novelty_assessment"),
            fto_status=pred_data.get("fto_status")
        )
        
        db.add(prediction)
        loaded_count += 1
    
    db.commit()
    print(f"✅ Loaded {loaded_count} dimeric predictions")
    
    return loaded_count


def main():
    """Main entry point."""
    print("🌿 NeuroBotanica Compound Library Loader")
    print("=" * 50)
    
    # Initialize database
    print("📦 Creating database tables...")
    Base.metadata.create_all(bind=engine)
    
    project_root = Path(__file__).parent.parent
    
    # Find compound JSON
    compound_paths = [
        project_root / "data" / "training" / "neurobotanica_complete_dataset_63compounds.json",
        project_root / "data" / "processed" / "neurobotanica_descriptors_validated.json"
    ]
    
    # Find dimeric predictions JSON
    dimeric_paths = [
        project_root / "data" / "training" / "neurobotanica_dimeric_predictions.json"
    ]
    
    db = SessionLocal()
    
    try:
        # Load compounds
        for path in compound_paths:
            if path.exists():
                load_compounds(str(path), db)
                break
        else:
            print("⚠️  No compound JSON found, skipping compounds")
        
        # Load dimeric predictions
        for path in dimeric_paths:
            if path.exists():
                load_dimeric_predictions(str(path), db)
                break
        else:
            print("⚠️  No dimeric predictions JSON found, skipping")
        
        # Print summary
        compound_count = db.query(Cannabinoid).count()
        dimer_count = db.query(DimericPrediction).count()
        fda_count = db.query(Cannabinoid).filter(Cannabinoid.fda_approved == True).count()
        
        print("\n" + "=" * 50)
        print("📊 DATABASE SUMMARY")
        print("=" * 50)
        print(f"Total compounds: {compound_count}")
        print(f"FDA-approved: {fda_count}")
        print(f"Dimeric predictions: {dimer_count}")
        
    finally:
        db.close()
    
    print("\n✅ Compound loading complete!")


if __name__ == "__main__":
    main()
