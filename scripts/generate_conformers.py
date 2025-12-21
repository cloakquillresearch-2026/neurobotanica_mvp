"""
Generate 3D Conformers for All Cannabinoids
NeuroBotanica MVP Development Plan - Week 1 Task 4.3

This script generates conformer ensembles for all cannabinoids in the database
using the ETKDG method with MMFF94 optimization.

Usage:
    python scripts/generate_conformers.py
    python scripts/generate_conformers.py --limit 5  # Test with 5 compounds
    python scripts/generate_conformers.py --compound "THC"  # Single compound
"""
import sys
import os
import argparse
from datetime import datetime
from pathlib import Path

# Add backend to path
sys.path.insert(0, str(Path(__file__).parent.parent / "backend"))

from sqlalchemy.orm import Session
from models.database import SessionLocal, engine, Base
from models.compound import Cannabinoid

# Import conformer generator
try:
    from services.conformer_generator import ConformerGenerator, ConformerResult
    GENERATOR_AVAILABLE = True
except ImportError as e:
    print(f"Warning: Could not import ConformerGenerator: {e}")
    GENERATOR_AVAILABLE = False


def generate_all_conformers(
    num_conformers: int = 50,
    limit: Optional[int] = None,
    compound_name: Optional[str] = None,
    force_regenerate: bool = False
) -> dict:
    """Generate 3D conformers for all cannabinoids in database.
    
    Args:
        num_conformers: Number of conformers to generate per compound
        limit: Optional limit on number of compounds to process
        compound_name: Optional specific compound name to process
        force_regenerate: If True, regenerate conformers even if they exist
        
    Returns:
        Summary dictionary with results
    """
    if not GENERATOR_AVAILABLE:
        print("❌ ConformerGenerator not available. Please install RDKit:")
        print("   pip install rdkit")
        return {"error": "RDKit not available"}
    
    db = SessionLocal()
    generator = ConformerGenerator(num_conformers=num_conformers)
    
    try:
        # Build query
        query = db.query(Cannabinoid)
        
        if compound_name:
            query = query.filter(
                (Cannabinoid.name.ilike(f"%{compound_name}%")) |
                (Cannabinoid.abbreviation.ilike(f"%{compound_name}%"))
            )
        elif not force_regenerate:
            query = query.filter(Cannabinoid.has_conformers == False)
        
        if limit:
            query = query.limit(limit)
        
        cannabinoids = query.all()
        
        if len(cannabinoids) == 0:
            print("ℹ️  No cannabinoids found matching criteria.")
            if not force_regenerate:
                print("   Use --force to regenerate existing conformers.")
            return {"message": "No compounds to process"}
        
        print(f"\n{'='*70}")
        print(f"🧬 NEUROBOTANICA 3D CONFORMER GENERATION")
        print(f"{'='*70}")
        print(f"Compounds to process: {len(cannabinoids)}")
        print(f"Conformers per compound: {num_conformers}")
        print(f"Method: ETKDG + MMFF94 optimization")
        print(f"{'='*70}\n")
        
        success_count = 0
        failed_count = 0
        failed_compounds = []
        results = []
        
        start_time = datetime.now()
        
        for idx, cb in enumerate(cannabinoids, 1):
            print(f"\n[{idx}/{len(cannabinoids)}] Processing: {cb.name}")
            print(f"    Abbreviation: {cb.abbreviation or 'N/A'}")
            print(f"    SMILES: {cb.smiles[:60]}{'...' if len(cb.smiles) > 60 else ''}")
            
            compound_start = datetime.now()
            result = generator.generate_conformers(cb.smiles)
            compound_time = (datetime.now() - compound_start).total_seconds()
            
            if result.success:
                # Update database
                cb.has_conformers = True
                cb.conformer_generation_method = "ETKDG"
                cb.num_conformers_generated = result.num_conformers_generated
                cb.conformer_metadata = result.conformer_metadata
                cb.rdkit_descriptors_3d = result.descriptors_3d
                cb.conformers_generated_at = datetime.now()
                
                db.commit()
                
                print(f"    ✅ SUCCESS: {result.num_conformers_generated} conformers")
                print(f"    📊 Energy range: {result.energy_range:.2f} kcal/mol")
                print(f"    🔮 RMSD clusters: {result.num_clusters_rmsd_2A}")
                print(f"    ⏱️  Time: {compound_time:.1f}s")
                
                # Show some 3D descriptors
                if result.descriptors_3d:
                    print(f"    📐 Asphericity: {result.descriptors_3d.get('asphericity', 'N/A')}")
                    print(f"    📐 Radius of gyration: {result.descriptors_3d.get('radius_of_gyration', 'N/A')}")
                
                success_count += 1
                results.append({
                    "name": cb.name,
                    "success": True,
                    "conformers": result.num_conformers_generated,
                    "time_seconds": compound_time
                })
            else:
                print(f"    ❌ FAILED: {result.error}")
                failed_count += 1
                failed_compounds.append(cb.name)
                results.append({
                    "name": cb.name,
                    "success": False,
                    "error": result.error,
                    "time_seconds": compound_time
                })
        
        total_time = (datetime.now() - start_time).total_seconds()
        
        # Print summary
        print(f"\n{'='*70}")
        print(f"📋 CONFORMER GENERATION SUMMARY")
        print(f"{'='*70}")
        print(f"✅ Successful: {success_count}/{len(cannabinoids)}")
        print(f"❌ Failed: {failed_count}/{len(cannabinoids)}")
        print(f"⏱️  Total time: {total_time:.1f}s ({total_time/60:.1f} minutes)")
        print(f"📊 Average time per compound: {total_time/len(cannabinoids):.1f}s")
        
        if failed_compounds:
            print(f"\n⚠️  Failed compounds:")
            for name in failed_compounds:
                print(f"   - {name}")
        
        # Show database summary
        total_with_conformers = db.query(Cannabinoid).filter(
            Cannabinoid.has_conformers == True
        ).count()
        total_compounds = db.query(Cannabinoid).count()
        
        print(f"\n📈 Database Status:")
        print(f"   Compounds with conformers: {total_with_conformers}/{total_compounds}")
        print(f"   Coverage: {total_with_conformers/total_compounds*100:.1f}%")
        print(f"{'='*70}\n")
        
        return {
            "total_processed": len(cannabinoids),
            "success_count": success_count,
            "failed_count": failed_count,
            "failed_compounds": failed_compounds,
            "total_time_seconds": total_time,
            "results": results
        }
        
    finally:
        db.close()


def show_conformer_status():
    """Show current conformer generation status for all compounds."""
    db = SessionLocal()
    
    try:
        print(f"\n{'='*70}")
        print(f"📊 CONFORMER STATUS REPORT")
        print(f"{'='*70}")
        
        compounds = db.query(Cannabinoid).order_by(Cannabinoid.name).all()
        
        with_conformers = 0
        without_conformers = 0
        
        print(f"\n{'Name':<30} {'Abbrev':<10} {'Conformers':<12} {'Method':<10}")
        print(f"{'-'*30} {'-'*10} {'-'*12} {'-'*10}")
        
        for cb in compounds:
            status = "✅" if cb.has_conformers else "⏳"
            conformer_count = cb.num_conformers_generated or 0
            method = cb.conformer_generation_method or "-"
            
            print(f"{cb.name[:30]:<30} {(cb.abbreviation or '-')[:10]:<10} "
                  f"{status} {conformer_count:<10} {method:<10}")
            
            if cb.has_conformers:
                with_conformers += 1
            else:
                without_conformers += 1
        
        print(f"\n📈 Summary:")
        print(f"   ✅ With conformers: {with_conformers}")
        print(f"   ⏳ Without conformers: {without_conformers}")
        print(f"   📊 Total compounds: {len(compounds)}")
        print(f"{'='*70}\n")
        
    finally:
        db.close()


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Generate 3D conformers for NeuroBotanica cannabinoid database"
    )
    parser.add_argument(
        "--num-conformers", "-n",
        type=int,
        default=50,
        help="Number of conformers to generate per compound (default: 50)"
    )
    parser.add_argument(
        "--limit", "-l",
        type=int,
        default=None,
        help="Limit number of compounds to process"
    )
    parser.add_argument(
        "--compound", "-c",
        type=str,
        default=None,
        help="Process specific compound by name"
    )
    parser.add_argument(
        "--force", "-f",
        action="store_true",
        help="Force regeneration of existing conformers"
    )
    parser.add_argument(
        "--status", "-s",
        action="store_true",
        help="Show conformer generation status only"
    )
    
    args = parser.parse_args()
    
    if args.status:
        show_conformer_status()
        return
    
    # Run conformer generation
    result = generate_all_conformers(
        num_conformers=args.num_conformers,
        limit=args.limit,
        compound_name=args.compound,
        force_regenerate=args.force
    )
    
    if "error" in result:
        sys.exit(1)


# Make Optional available
from typing import Optional

if __name__ == "__main__":
    main()
