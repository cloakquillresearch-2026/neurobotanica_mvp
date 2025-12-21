#!/usr/bin/env python3
"""
Generate Dimeric Predictions Script
Creates 120+ dimeric cannabinoid predictions with triangulation scoring

Supports NeuroBotanica patent claims:
- Novel dimeric cannabinoid structures
- Synergy prediction for compound combinations
- Therapeutic potential predictions
- Multi-source validation triangulation

Reference: NeuroBotanica MVP Development Plan - Week 4
"""
import json
import sys
import os
import logging
from pathlib import Path
from datetime import datetime
from typing import Dict, Any, List, Optional

# Add backend to path
sys.path.insert(0, str(Path(__file__).parent.parent / "backend"))

from services.dimer_predictor import (
    DimericPredictor,
    DimerPrediction,
    LinkageType,
    DimerType
)
from services.triangulation_scorer import (
    TriangulationScorer,
    TriangulationResult
)
from services.dimer_conformer_generator import (
    DimericConformerGenerator,
    DimericConformerResult
)

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def load_compound_data(data_path: Path) -> Dict[str, str]:
    """Load compound SMILES from the training dataset.
    
    Args:
        data_path: Path to the compound dataset JSON
        
    Returns:
        Dictionary of compound_name -> SMILES
    """
    with open(data_path, 'r') as f:
        data = json.load(f)
    
    compounds = {}
    for compound in data.get('compounds', []):
        name = compound.get('compound_name')
        smiles = compound.get('smiles')
        compound_type = compound.get('compound_type', '')
        
        # Only include cannabinoids for dimer predictions
        if name and smiles and compound_type == 'cannabinoid':
            compounds[name] = smiles
    
    logger.info(f"Loaded {len(compounds)} cannabinoids from dataset")
    return compounds


def generate_predictions(
    compounds: Dict[str, str],
    include_homodimers: bool = True,
    include_heterodimers: bool = True
) -> List[DimerPrediction]:
    """Generate all dimer predictions.
    
    Args:
        compounds: Dictionary of compound_name -> SMILES
        include_homodimers: Whether to include A-A dimers
        include_heterodimers: Whether to include A-B dimers
        
    Returns:
        List of DimerPrediction objects
    """
    predictor = DimericPredictor()
    
    logger.info(f"Generating predictions for {len(compounds)} cannabinoids...")
    logger.info(f"Homodimers: {include_homodimers}, Heterodimers: {include_heterodimers}")
    
    # Convert dict to list of tuples for the API
    monomers_list = [(name, smiles) for name, smiles in compounds.items()]
    
    # Generate predictions (API only supports include_homodimers flag)
    # Heterodimers are always generated along with homodimers
    predictions = predictor.generate_all_combinations(
        monomers_list,
        include_homodimers=include_homodimers
    )
    
    logger.info(f"Generated {len(predictions)} total predictions")
    return predictions


def calculate_triangulation_scores(
    predictions: List[DimerPrediction]
) -> List[Dict[str, Any]]:
    """Calculate triangulation scores for all predictions.
    
    Args:
        predictions: List of dimer predictions
        
    Returns:
        List of predictions with triangulation results
    """
    scorer = TriangulationScorer()
    results = []
    
    for pred in predictions:
        # Calculate therapeutic score (average of top 3 indications)
        therapeutic_scores = list(pred.therapeutic_potential.values())
        top_therapeutic = sorted(therapeutic_scores, reverse=True)[:3]
        ml_therapeutic = sum(top_therapeutic) / len(top_therapeutic) if top_therapeutic else 0.5
        
        # Calculate triangulation
        tri_result = scorer.calculate_triangulation_score(
            formation_probability=pred.formation_probability,
            structural_similarity=1 - pred.novelty_score,  # Lower novelty = higher similarity
            ml_therapeutic_score=ml_therapeutic,
            experimental_validation="predicted",  # All are predictions
            literature_score=pred.synergy_prediction * 0.8  # Synergy from literature
        )
        
        results.append({
            "prediction": pred,
            "triangulation": tri_result
        })
    
    return results


def generate_conformers_for_top_dimers(
    predictions_with_tri: List[Dict[str, Any]],
    top_n: int = 20
) -> List[Dict[str, Any]]:
    """Generate 3D conformers for top-ranked dimers.
    
    Args:
        predictions_with_tri: Predictions with triangulation scores
        top_n: Number of top dimers to generate conformers for
        
    Returns:
        List with conformer data added
    """
    generator = DimericConformerGenerator()
    
    # Sort by triangulation score
    sorted_preds = sorted(
        predictions_with_tri,
        key=lambda x: x['triangulation'].triangulation_score,
        reverse=True
    )[:top_n]
    
    logger.info(f"Generating conformers for top {len(sorted_preds)} dimers...")
    
    for item in sorted_preds:
        pred = item['prediction']
        
        try:
            # Get monomer SMILES - use parent_1_smiles and parent_2_smiles
            monomer_a = pred.parent_1_smiles if hasattr(pred, 'parent_1_smiles') else None
            monomer_b = pred.parent_2_smiles if hasattr(pred, 'parent_2_smiles') else None
            
            # Skip if no valid SMILES
            if not pred.predicted_smiles:
                item['conformers'] = None
                logger.warning(f"No SMILES for {pred.dimer_name}")
                continue
            
            # Use predicted_smiles for the dimer structure
            result = generator.generate_dimer_conformers(
                dimer_smiles=pred.predicted_smiles,
                monomer1_smiles=monomer_a,
                monomer2_smiles=monomer_b
            )
            
            if result and result.success:
                # Access attributes through base_result
                base = result.base_result
                item['conformers'] = {
                    'num_conformers': base.num_conformers_generated,
                    'successful': base.num_conformers_generated,
                    'lowest_energy': base.lowest_energy,
                    'energies': [],  # Not stored individually
                    'best_for_binding': result.best_for_binding,
                    'best_for_stability': result.best_for_stability
                }
                logger.info(f"Generated {base.num_conformers_generated} conformers for {pred.dimer_name}")
            else:
                item['conformers'] = None
                logger.warning(f"Could not generate conformers for {pred.dimer_name}")
                
        except Exception as e:
            item['conformers'] = None
            logger.error(f"Error generating conformers for {pred.dimer_name}: {e}")
    
    return sorted_preds


def format_output(
    predictions_with_tri: List[Dict[str, Any]],
    top_with_conformers: List[Dict[str, Any]]
) -> Dict[str, Any]:
    """Format predictions for JSON output.
    
    Args:
        predictions_with_tri: All predictions with triangulation
        top_with_conformers: Top predictions with conformer data
        
    Returns:
        Formatted dictionary for JSON serialization
    """
    # Create ID mapping for top predictions
    top_ids = {pred['prediction'].dimer_name for pred in top_with_conformers}
    
    # Format all predictions
    formatted_predictions = []
    
    for item in predictions_with_tri:
        pred = item['prediction']
        tri = item['triangulation']
        
        formatted = {
            "id": len(formatted_predictions) + 1,
            "dimer_name": pred.dimer_name,
            "dimer_smiles": pred.predicted_smiles,
            "dimer_type": pred.dimer_type.value,
            "linkage_type": pred.linkage_type.value,
            "monomer_a": pred.parent_1_name,
            "monomer_b": pred.parent_2_name,
            "scores": {
                "formation_probability": round(pred.formation_probability, 4),
                "synergy_score": round(pred.synergy_prediction, 4),
                "novelty_score": round(pred.novelty_score, 4),
                "structural_validity": round(pred.structural_validity, 4)
            },
            "triangulation": {
                "score": tri.triangulation_score,
                "uncertainty": tri.uncertainty,
                "confidence_interval": {
                    "lower": tri.confidence_interval_lower,
                    "upper": tri.confidence_interval_upper
                },
                "grade": tri.validation_grade,
                "num_sources": tri.num_sources,
                "consensus": tri.consensus_strength
            },
            "predicted_properties": {
                "molecular_weight": pred.predicted_mw,
                "logp": pred.predicted_logp,
                "cb1_affinity": pred.predicted_cb1_affinity,
                "cb2_affinity": pred.predicted_cb2_affinity,
                "selectivity": pred.predicted_selectivity
            },
            "therapeutic_potential": {
                k: round(v, 3) for k, v in pred.therapeutic_potential.items()
            }
        }
        
        # Add conformer data if available
        if pred.dimer_name in top_ids:
            for top_item in top_with_conformers:
                if top_item['prediction'].dimer_name == pred.dimer_name:
                    formatted['conformers'] = top_item.get('conformers')
                    formatted['has_3d_structure'] = top_item.get('conformers') is not None
                    break
        else:
            formatted['has_3d_structure'] = False
        
        formatted_predictions.append(formatted)
    
    # Sort by triangulation score
    formatted_predictions.sort(
        key=lambda x: x['triangulation']['score'],
        reverse=True
    )
    
    # Add rank
    for i, pred in enumerate(formatted_predictions, 1):
        pred['rank'] = i
    
    # Calculate summary statistics
    scores = [p['triangulation']['score'] for p in formatted_predictions]
    grades = [p['triangulation']['grade'] for p in formatted_predictions]
    
    return {
        "metadata": {
            "version": "1.0",
            "generated_at": datetime.now().isoformat(),
            "total_predictions": len(formatted_predictions),
            "predictions_with_conformers": len(top_with_conformers),
            "generation_method": "DimericPredictor v1.0"
        },
        "summary": {
            "avg_triangulation_score": round(sum(scores) / len(scores), 4) if scores else 0,
            "min_score": round(min(scores), 4) if scores else 0,
            "max_score": round(max(scores), 4) if scores else 0,
            "grade_distribution": {
                "A": grades.count("A"),
                "B": grades.count("B"),
                "C": grades.count("C"),
                "D": grades.count("D")
            },
            "dimer_types": {
                "homodimer": sum(1 for p in formatted_predictions if p['dimer_type'] == 'homodimer'),
                "heterodimer": sum(1 for p in formatted_predictions if p['dimer_type'] == 'heterodimer')
            }
        },
        "predictions": formatted_predictions
    }


def main():
    """Main entry point for dimeric prediction generation."""
    logger.info("=" * 60)
    logger.info("NeuroBotanica Dimeric Prediction Generator")
    logger.info("=" * 60)
    
    # Define paths
    project_root = Path(__file__).parent.parent
    data_path = project_root / "data" / "training" / "neurobotanica_complete_dataset_63compounds.json"
    output_path = project_root / "data" / "training" / "neurobotanica_dimeric_predictions.json"
    
    # Step 1: Load compound data
    logger.info("\nStep 1: Loading compound data...")
    compounds = load_compound_data(data_path)
    
    if not compounds:
        logger.error("No compounds loaded! Check the data file.")
        return
    
    # Step 2: Generate predictions
    logger.info("\nStep 2: Generating dimer predictions...")
    predictions = generate_predictions(
        compounds,
        include_homodimers=True,
        include_heterodimers=True
    )
    
    if not predictions:
        logger.error("No predictions generated!")
        return
    
    # Step 3: Calculate triangulation scores
    logger.info("\nStep 3: Calculating triangulation scores...")
    predictions_with_tri = calculate_triangulation_scores(predictions)
    
    # Step 4: Generate conformers for top dimers
    logger.info("\nStep 4: Generating 3D conformers for top dimers...")
    top_with_conformers = generate_conformers_for_top_dimers(
        predictions_with_tri,
        top_n=20
    )
    
    # Step 5: Format and save output
    logger.info("\nStep 5: Formatting and saving results...")
    output = format_output(predictions_with_tri, top_with_conformers)
    
    with open(output_path, 'w') as f:
        json.dump(output, f, indent=2)
    
    logger.info(f"\nResults saved to: {output_path}")
    
    # Print summary
    logger.info("\n" + "=" * 60)
    logger.info("GENERATION COMPLETE")
    logger.info("=" * 60)
    logger.info(f"Total predictions: {output['metadata']['total_predictions']}")
    logger.info(f"Predictions with conformers: {output['metadata']['predictions_with_conformers']}")
    logger.info(f"Average triangulation score: {output['summary']['avg_triangulation_score']}")
    logger.info(f"Grade distribution: {output['summary']['grade_distribution']}")
    logger.info(f"Dimer types: {output['summary']['dimer_types']}")
    
    # Print top 10 predictions
    logger.info("\n" + "-" * 60)
    logger.info("TOP 10 PREDICTIONS")
    logger.info("-" * 60)
    
    for pred in output['predictions'][:10]:
        logger.info(
            f"#{pred['rank']}: {pred['dimer_name']} "
            f"(Grade {pred['triangulation']['grade']}, "
            f"Score: {pred['triangulation']['score']:.3f})"
        )


if __name__ == "__main__":
    main()
