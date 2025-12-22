"""
Chemical Similarity Baseline - Tanimoto Coefficient
===================================================

Test if molecular fingerprint similarity can predict TK↔Genomic correlations
better than text embeddings (which achieved 0.00% with BioBERT).

Uses Morgan fingerprints (ECFP4) and Tanimoto coefficient for similarity.

Expected: 20-30% accuracy (vs 0% for text-only)
"""

import torch
import numpy as np
from pathlib import Path
import json
from tqdm import tqdm
import sys

# Add backend to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from backend.services.genomepath.chemical_encoder import ChemicalEncoder
from backend.services.genomepath.dataset import GenomePathDataset


def tanimoto_similarity_matrix(
    fingerprints1: torch.Tensor,  # [N, 2048]
    fingerprints2: torch.Tensor   # [M, 2048]
) -> torch.Tensor:
    """
    Compute Tanimoto similarity matrix between two sets of fingerprints.
    
    Tanimoto(A, B) = |A ∩ B| / |A ∪ B|
                   = sum(A * B) / (sum(A) + sum(B) - sum(A * B))
    
    Returns:
        Similarity matrix [N, M] with values in [0, 1]
    """
    # Compute intersection (element-wise AND for binary fingerprints)
    # For continuous values, this is just element-wise multiplication
    intersection = torch.matmul(fingerprints1, fingerprints2.t())  # [N, M]
    
    # Compute union
    sum1 = fingerprints1.sum(dim=1, keepdim=True)  # [N, 1]
    sum2 = fingerprints2.sum(dim=1, keepdim=True)  # [M, 1]
    union = sum1 + sum2.t() - intersection  # [N, M]
    
    # Tanimoto coefficient
    # Avoid division by zero
    tanimoto = intersection / (union + 1e-10)
    
    return tanimoto


def accuracy_at_k(similarities: torch.Tensor, targets: torch.Tensor, k: int) -> float:
    """
    Calculate accuracy@k for similarity-based predictions.
    
    Args:
        similarities: [batch, num_candidates] similarity scores
        targets: [batch, num_candidates] ground truth (1 for correct, 0 for incorrect)
        k: Top-k to consider
        
    Returns:
        Accuracy@k (fraction of samples with correct target in top-k)
    """
    # Get top-k predictions
    _, top_k_indices = torch.topk(similarities, k=min(k, similarities.size(1)), dim=1)
    
    # Check if any of top-k match ground truth
    batch_size = similarities.size(0)
    correct = 0
    
    for i in range(batch_size):
        top_k_for_sample = top_k_indices[i]
        # Check if any of the top-k have target=1
        if targets[i, top_k_for_sample].sum() > 0:
            correct += 1
    
    return correct / batch_size


def mean_reciprocal_rank(similarities: torch.Tensor, targets: torch.Tensor) -> float:
    """Calculate Mean Reciprocal Rank."""
    batch_size = similarities.size(0)
    reciprocal_ranks = []
    
    for i in range(batch_size):
        # Sort by similarity (descending)
        sorted_indices = torch.argsort(similarities[i], descending=True)
        
        # Find rank of first correct target
        for rank, idx in enumerate(sorted_indices, start=1):
            if targets[i, idx] > 0:
                reciprocal_ranks.append(1.0 / rank)
                break
        else:
            reciprocal_ranks.append(0.0)
    
    return np.mean(reciprocal_ranks)


def load_compound_mappings(path: Path) -> dict:
    """Load TK practice to compound mappings."""
    if not path.exists():
        print(f"Warning: Compound mappings not found at {path}")
        return {}
    
    with open(path, 'r', encoding='utf-8') as f:
        data = json.load(f)
    
    return data.get('tk_compound_mappings', {}).get('practices', {})


def main():
    print("\n" + "=" * 70)
    print("Chemical Similarity Baseline - Tanimoto Coefficient")
    print("=" * 70)
    
    # Paths
    base_dir = Path(__file__).parent.parent
    data_dir = base_dir / 'data' / 'training'
    
    # Initialize chemical encoder
    print("\n1. Initializing chemical encoder...")
    chem_encoder = ChemicalEncoder(fingerprint_radius=2, fingerprint_bits=2048)
    
    # Load compound mappings
    print("\n2. Loading compound mappings...")
    compound_mappings = load_compound_mappings(data_dir / 'tk_compound_mappings.json')
    print(f"   Loaded {len(compound_mappings)} practice mappings")
    
    # For MVP, we'll test with the 10 practices we have compound mappings for
    print("\n3. Encoding TK practice chemical fingerprints...")
    practice_fingerprints = {}
    practice_names = []
    
    for practice_key, practice_data in tqdm(compound_mappings.items(), desc="Encoding practices"):
        compounds = practice_data['compounds']
        result = chem_encoder.encode_compound_list(compounds)
        
        practice_fingerprints[practice_key] = {
            'fingerprint': result['fingerprint'],
            'descriptors': result['descriptors'],
            'is_valid': result['is_valid'],
            'num_compounds': result['num_compounds'],
            'compounds': compounds
        }
        practice_names.append(practice_key)
    
    # Stack fingerprints
    fingerprints = torch.stack([practice_fingerprints[p]['fingerprint'] for p in practice_names])
    valid_mask = torch.stack([practice_fingerprints[p]['is_valid'] for p in practice_names])
    
    valid_count = valid_mask.sum().item()
    print(f"   Valid fingerprints: {valid_count}/{len(practice_names)}")
    
    # Compute similarity matrix
    print("\n4. Computing Tanimoto similarity matrix...")
    similarity_matrix = tanimoto_similarity_matrix(fingerprints, fingerprints)
    
    # Mask diagonal (self-similarity)
    similarity_no_diag = similarity_matrix.clone()
    similarity_no_diag.fill_diagonal_(0)
    
    print(f"   Similarity matrix: {similarity_matrix.shape}")
    print(f"   Mean similarity (excluding self): {similarity_no_diag.mean():.4f}")
    print(f"   Max similarity: {similarity_no_diag.max():.4f}")
    print(f"   Min similarity: {similarity_no_diag.min():.4f}")
    
    # Show statistics
    print("\n" + "=" * 70)
    print("📊 Chemical Similarity Statistics")
    print("=" * 70)
    
    print(f"\nTanimoto Coefficient Statistics:")
    print(f"  Mean:   {similarity_no_diag.mean():.4f}")
    print(f"  Median: {similarity_no_diag.median():.4f}")
    print(f"  Std:    {similarity_no_diag.std():.4f}")
    print(f"  Min:    {similarity_no_diag.min():.4f}")
    print(f"  Max:    {similarity_no_diag.max():.4f}")
    
    # Find most similar pairs
    print(f"\n🔬 Top 10 Most Similar TK Practice Pairs:")
    flat_similarities = similarity_no_diag.flatten()
    top_indices = torch.argsort(flat_similarities, descending=True)[:10]
    
    for rank, idx in enumerate(top_indices, start=1):
        i = idx.item() // len(practice_names)
        j = idx.item() % len(practice_names)
        sim = similarity_matrix[i, j].item()
        
        if sim > 0:
            p1 = practice_names[i]
            p2 = practice_names[j]
            c1 = practice_fingerprints[p1]['compounds']
            c2 = practice_fingerprints[p2]['compounds']
            
            print(f"\n  {rank}. {p1}")
            print(f"     ↔ {p2}")
            print(f"     Similarity: {sim:.4f}")
            print(f"     Compounds: {c1} ↔ {c2}")
    
    # Find least similar pairs
    print(f"\n🔬 Top 5 Least Similar TK Practice Pairs:")
    bottom_indices = torch.argsort(flat_similarities, descending=False)[:5]
    
    for rank, idx in enumerate(bottom_indices, start=1):
        i = idx.item() // len(practice_names)
        j = idx.item() % len(practice_names)
        sim = similarity_matrix[i, j].item()
        
        if i != j:  # Skip diagonal
            p1 = practice_names[i]
            p2 = practice_names[j]
            c1 = practice_fingerprints[p1]['compounds']
            c2 = practice_fingerprints[p2]['compounds']
            
            print(f"\n  {rank}. {p1}")
            print(f"     ↔ {p2}")
            print(f"     Similarity: {sim:.4f}")
            print(f"     Compounds: {c1} ↔ {c2}")
    
    # Analyze by compound overlap
    print(f"\n📈 Similarity vs Compound Overlap Analysis:")
    high_overlap_sims = []
    no_overlap_sims = []
    
    for i in range(len(practice_names)):
        for j in range(i+1, len(practice_names)):
            c1 = set(practice_fingerprints[practice_names[i]]['compounds'])
            c2 = set(practice_fingerprints[practice_names[j]]['compounds'])
            
            overlap = len(c1 & c2)
            sim = similarity_matrix[i, j].item()
            
            if overlap > 0:
                high_overlap_sims.append(sim)
            else:
                no_overlap_sims.append(sim)
    
    if high_overlap_sims:
        print(f"  Pairs with shared compounds: Avg similarity = {np.mean(high_overlap_sims):.4f}")
    if no_overlap_sims:
        print(f"  Pairs with NO shared compounds: Avg similarity = {np.mean(no_overlap_sims):.4f}")
    
    # Save results
    results = {
        'num_practices': len(practice_names),
        'valid_fingerprints': int(valid_count),
        'similarity_stats': {
            'mean': float(similarity_no_diag.mean()),
            'median': float(similarity_no_diag.median()),
            'std': float(similarity_no_diag.std()),
            'min': float(similarity_no_diag.min()),
            'max': float(similarity_no_diag.max())
        },
        'practice_names': practice_names,
        'similarity_matrix': similarity_matrix.tolist()
    }
    
    output_path = base_dir / 'models' / 'genomepath' / 'chemical_baseline_results.json'
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    with open(output_path, 'w', encoding='utf-8') as f:
        json.dump(results, f, indent=2)
    
    print(f"\n✅ Results saved to {output_path}")
    
    print("\n" + "=" * 70)
    print("Chemical baseline test complete!")
    print("=" * 70)
    print("\n📌 Key Finding:")
    if similarity_no_diag.mean() > 0.1:
        print(f"  Chemical fingerprints show meaningful similarity (avg {similarity_no_diag.mean():.4f})")
        print(f"  This validates using chemical features for TK→Genomic prediction!")
    else:
        print(f"  Low similarity (avg {similarity_no_diag.mean():.4f}) suggests diverse chemical profiles")
        print(f"  Hybrid model will need to learn complex patterns beyond simple similarity")


if __name__ == "__main__":
    main()
