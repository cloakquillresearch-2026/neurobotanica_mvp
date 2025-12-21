#!/usr/bin/env python
"""
External Database Sync Script - NeuroBotanica Week 5
Syncs cannabinoid data with ChEMBL and PubChem databases.

Usage:
    python scripts/sync_external_databases.py [--compounds N] [--dry-run]
"""
import sys
import os
import json
import argparse
from datetime import datetime
from typing import Dict, List, Optional
import logging

# Add project root to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from backend.services.chembl_client import ChEMBLClient
from backend.services.pubchem_client import PubChemClient

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class ExternalDatabaseSync:
    """Sync cannabinoid data with external databases."""
    
    def __init__(self, dry_run: bool = False):
        """Initialize sync service.
        
        Args:
            dry_run: If True, don't write changes
        """
        self.chembl = ChEMBLClient()
        self.pubchem = PubChemClient()
        self.dry_run = dry_run
        
        # Load training dataset
        self.dataset_path = os.path.join(
            os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
            "data", "training", "neurobotanica_complete_dataset_63compounds.json"
        )
        
        # Output path for enriched data
        self.output_path = os.path.join(
            os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
            "data", "processed", "cannabinoids_enriched.json"
        )
        
        self.stats = {
            "total_compounds": 0,
            "chembl_matches": 0,
            "pubchem_matches": 0,
            "activities_added": 0,
            "properties_added": 0,
            "errors": []
        }
    
    def load_dataset(self) -> List[Dict]:
        """Load cannabinoid dataset."""
        if not os.path.exists(self.dataset_path):
            logger.error(f"Dataset not found: {self.dataset_path}")
            return []
        
        with open(self.dataset_path, 'r') as f:
            data = json.load(f)
        
        compounds = data.get("compounds", data.get("cannabinoids", []))
        self.stats["total_compounds"] = len(compounds)
        logger.info(f"Loaded {len(compounds)} compounds from dataset")
        return compounds
    
    def sync_with_chembl(self, compound: Dict) -> Dict:
        """Sync single compound with ChEMBL.
        
        Args:
            compound: Compound dictionary with name/smiles
            
        Returns:
            Updated compound dict with ChEMBL data
        """
        name = compound.get("name", "Unknown")
        smiles = compound.get("smiles", "")
        
        logger.info(f"Syncing {name} with ChEMBL...")
        
        # Search by SMILES first, then name
        chembl_id = None
        
        if smiles:
            chembl_id = self.chembl.search_by_smiles(smiles)
        
        if not chembl_id:
            # Try known cannabinoid IDs
            known_ids = self.chembl.get_known_cannabinoid_ids()
            name_upper = name.upper()
            for known_name, known_id in known_ids.items():
                if known_name.upper() in name_upper or name_upper in known_name.upper():
                    chembl_id = known_id
                    break
        
        if not chembl_id:
            # Try name search
            chembl_id = self.chembl.search_by_name(name)
        
        if chembl_id:
            self.stats["chembl_matches"] += 1
            compound["chembl_id"] = chembl_id
            
            # Get molecule details
            molecule = self.chembl.get_molecule(chembl_id)
            if molecule:
                compound["chembl_data"] = molecule.to_dict()
            
            # Get cannabinoid receptor activities
            activities = self.chembl.get_cannabinoid_target_activities(chembl_id)
            if activities:
                compound["chembl_activities"] = {}
                for target, acts in activities.items():
                    compound["chembl_activities"][target] = [
                        a.to_provenance_dict() for a in acts
                    ]
                    self.stats["activities_added"] += len(acts)
                    logger.info(f"  Found {len(acts)} {target} activities")
        else:
            logger.info(f"  No ChEMBL match found")
        
        return compound
    
    def sync_with_pubchem(self, compound: Dict) -> Dict:
        """Sync single compound with PubChem.
        
        Args:
            compound: Compound dictionary with name/smiles
            
        Returns:
            Updated compound dict with PubChem data
        """
        name = compound.get("name", "Unknown")
        smiles = compound.get("smiles", "")
        
        logger.info(f"Syncing {name} with PubChem...")
        
        # Search by SMILES first
        cid = None
        
        if smiles:
            cid = self.pubchem.search_by_smiles(smiles)
        
        if not cid:
            # Try known cannabinoid CIDs
            known_cids = self.pubchem.get_known_cannabinoid_cids()
            name_upper = name.upper()
            for known_name, known_cid in known_cids.items():
                if known_name.upper() in name_upper or name_upper in known_name.upper():
                    cid = known_cid
                    break
        
        if not cid:
            # Try name search
            cid = self.pubchem.search_by_name(name)
        
        if cid:
            self.stats["pubchem_matches"] += 1
            compound["pubchem_cid"] = cid
            
            # Get properties
            properties = self.pubchem.get_compound_properties(cid)
            if properties:
                compound["pubchem_properties"] = properties.to_dict()
                self.stats["properties_added"] += 1
            
            # Get synonyms
            synonyms = self.pubchem.get_synonyms(cid)
            if synonyms:
                compound["pubchem_synonyms"] = synonyms.get_common_names(20)
            
            # Get cannabinoid-relevant bioassays
            bioassays = self.pubchem.get_cannabinoid_relevant_assays(cid)
            if bioassays:
                compound["pubchem_bioassays"] = [
                    a.to_provenance_dict() for a in bioassays[:50]
                ]
                logger.info(f"  Found {len(bioassays)} relevant bioassays")
        else:
            logger.info(f"  No PubChem match found")
        
        return compound
    
    def merge_receptor_affinities(self, compound: Dict) -> Dict:
        """Merge receptor affinity data from multiple sources.
        
        Creates a unified receptor_affinities_enriched field with
        provenance from ChEMBL activities.
        """
        enriched_affinities = {}
        
        # Get existing affinities
        existing = compound.get("receptor_affinities", {})
        
        # Get ChEMBL activities
        chembl_activities = compound.get("chembl_activities", {})
        
        # Merge CB1 data
        cb1_data = []
        if "CB1" in existing:
            cb1_data.append({
                "affinity_value": existing["CB1"],
                "affinity_unit": "Ki (nM)",
                "source": "original_dataset"
            })
        if "CB1" in chembl_activities:
            cb1_data.extend(chembl_activities["CB1"])
        if cb1_data:
            enriched_affinities["CB1"] = cb1_data
        
        # Merge CB2 data
        cb2_data = []
        if "CB2" in existing:
            cb2_data.append({
                "affinity_value": existing["CB2"],
                "affinity_unit": "Ki (nM)",
                "source": "original_dataset"
            })
        if "CB2" in chembl_activities:
            cb2_data.extend(chembl_activities["CB2"])
        if cb2_data:
            enriched_affinities["CB2"] = cb2_data
        
        # Add other targets from ChEMBL
        for target in ["GPR55", "TRPV1", "FAAH", "MAGL", "PPARα", "PPARγ", "5-HT1A"]:
            if target in chembl_activities:
                enriched_affinities[target] = chembl_activities[target]
        
        if enriched_affinities:
            compound["receptor_affinities_enriched"] = enriched_affinities
        
        return compound
    
    def sync_all(self, max_compounds: Optional[int] = None) -> Dict:
        """Sync all compounds with external databases.
        
        Args:
            max_compounds: Maximum number of compounds to sync (for testing)
            
        Returns:
            Statistics dictionary
        """
        logger.info("=" * 60)
        logger.info("NeuroBotanica External Database Sync")
        logger.info("=" * 60)
        
        compounds = self.load_dataset()
        if not compounds:
            return self.stats
        
        if max_compounds:
            compounds = compounds[:max_compounds]
            logger.info(f"Limited to {max_compounds} compounds")
        
        enriched_compounds = []
        
        for i, compound in enumerate(compounds, 1):
            name = compound.get("name", "Unknown")
            logger.info(f"\n[{i}/{len(compounds)}] Processing: {name}")
            
            try:
                # Sync with ChEMBL
                compound = self.sync_with_chembl(compound)
                
                # Sync with PubChem
                compound = self.sync_with_pubchem(compound)
                
                # Merge receptor affinities
                compound = self.merge_receptor_affinities(compound)
                
                # Add sync timestamp
                compound["last_external_sync"] = datetime.now().isoformat()
                
                enriched_compounds.append(compound)
                
            except Exception as e:
                logger.error(f"Error processing {name}: {e}")
                self.stats["errors"].append({"compound": name, "error": str(e)})
                enriched_compounds.append(compound)
        
        # Save enriched data
        if not self.dry_run:
            self._save_enriched_data(enriched_compounds)
        
        # Print summary
        self._print_summary()
        
        return self.stats
    
    def _save_enriched_data(self, compounds: List[Dict]):
        """Save enriched compound data to file."""
        output_data = {
            "metadata": {
                "generated_at": datetime.now().isoformat(),
                "source_dataset": self.dataset_path,
                "num_compounds": len(compounds),
                "chembl_matches": self.stats["chembl_matches"],
                "pubchem_matches": self.stats["pubchem_matches"]
            },
            "compounds": compounds
        }
        
        os.makedirs(os.path.dirname(self.output_path), exist_ok=True)
        
        with open(self.output_path, 'w') as f:
            json.dump(output_data, f, indent=2, default=str)
        
        logger.info(f"\nEnriched data saved to: {self.output_path}")
    
    def _print_summary(self):
        """Print sync summary."""
        print("\n" + "=" * 60)
        print("SYNC SUMMARY")
        print("=" * 60)
        print(f"Total compounds:    {self.stats['total_compounds']}")
        print(f"ChEMBL matches:     {self.stats['chembl_matches']}")
        print(f"PubChem matches:    {self.stats['pubchem_matches']}")
        print(f"Activities added:   {self.stats['activities_added']}")
        print(f"Properties added:   {self.stats['properties_added']}")
        print(f"Errors:             {len(self.stats['errors'])}")
        
        if self.stats['errors']:
            print("\nErrors:")
            for err in self.stats['errors'][:5]:
                print(f"  - {err['compound']}: {err['error']}")
        
        print("=" * 60)


def main():
    parser = argparse.ArgumentParser(
        description="Sync cannabinoid data with ChEMBL and PubChem"
    )
    parser.add_argument(
        "--compounds", "-n",
        type=int,
        default=None,
        help="Maximum number of compounds to sync (for testing)"
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Don't save output file"
    )
    
    args = parser.parse_args()
    
    sync = ExternalDatabaseSync(dry_run=args.dry_run)
    stats = sync.sync_all(max_compounds=args.compounds)
    
    # Exit with error code if too many failures
    error_rate = len(stats['errors']) / max(stats['total_compounds'], 1)
    if error_rate > 0.5:
        logger.error("Too many errors during sync")
        sys.exit(1)


if __name__ == "__main__":
    main()
