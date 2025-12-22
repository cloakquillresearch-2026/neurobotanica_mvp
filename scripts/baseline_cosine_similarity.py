"""
Simple Cosine Similarity Baseline for GenomePath
=================================================

Tests if simple cosine similarity between SentenceTransformer embeddings
can match or beat the 30M parameter transformer model.

If this baseline beats 13.47% accuracy, the transformer is overkill.
If this also fails, data quality is the issue.
"""

import sys
from pathlib import Path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

import torch
import torch.nn.functional as F
import json
import numpy as np
from tqdm import tqdm

from backend.services.genomepath.tokenizer import GenomePathTokenizer
from backend.services.genomepath.dataset import GenomePathDataset, create_data_splits


def cosine_similarity_matrix(a, b):
    """
    Compute cosine similarity between all pairs in a and b.
    
    Args:
        a: [N, D] tensor
        b: [M, D] tensor
        
    Returns:
        [N, M] similarity matrix
    """
    a_norm = F.normalize(a, p=2, dim=1)
    b_norm = F.normalize(b, p=2, dim=1)
    return torch.mm(a_norm, b_norm.t())


def accuracy_at_k(similarities, targets, k=5):
    """Calculate accuracy@k."""
    batch_size = similarities.size(0)
    
    # Get top-k predictions
    _, top_k_indices = torch.topk(similarities, k=min(k, similarities.size(1)), dim=1)
    
    correct = 0
    for i in range(batch_size):
        true_indices = torch.where(targets[i] > 0.5)[0]
        if len(true_indices) > 0:
            if torch.any(torch.isin(true_indices, top_k_indices[i])):
                correct += 1
    
    return correct / batch_size if batch_size > 0 else 0.0


def mean_reciprocal_rank(similarities, targets):
    """Calculate Mean Reciprocal Rank."""
    batch_size = similarities.size(0)
    reciprocal_ranks = []
    
    for i in range(batch_size):
        sorted_indices = torch.argsort(similarities[i], descending=True)
        true_indices = torch.where(targets[i] > 0.5)[0]
        
        if len(true_indices) > 0:
            for rank, idx in enumerate(sorted_indices, start=1):
                if idx in true_indices:
                    reciprocal_ranks.append(1.0 / rank)
                    break
        else:
            reciprocal_ranks.append(0.0)
    
    return np.mean(reciprocal_ranks) if reciprocal_ranks else 0.0


def main():
    print("=" * 70)
    print("Cosine Similarity Baseline Test")
    print("=" * 70)
    
    # Load dataset
    print("\nLoading data...")
    tokenizer = GenomePathTokenizer()
    dataset = GenomePathDataset(
        correlations_path='data/processed/training_correlations.json',
        tk_practices_path='data/processed/tk_practices.json',
        genomic_targets_path='data/processed/genomic_targets.json',
        tokenizer=tokenizer,
        augment=False
    )
    
    # Create splits
    train_ds, val_ds, test_ds = create_data_splits(dataset, seed=42)
    
    print(f"\nDataset splits:")
    print(f"  Train: {len(train_ds)} samples")
    print(f"  Val: {len(val_ds)} samples")
    print(f"  Test: {len(test_ds)} samples")
    
    # Pre-compute all TK practice embeddings (mean pool)
    print("\nEncoding TK practices...")
    tk_practice_ids = sorted(dataset.tk_practices.keys())
    tk_embeddings_list = []
    
    for practice_id in tqdm(tk_practice_ids):
        practice = dataset.tk_practices[practice_id]
        result = tokenizer.encode_tk_practice(
            practice_name=practice['practice_name'],
            indications=practice.get('indications', []),
            preparation_method=practice.get('preparation_method', ''),
            cultural_context=practice.get('cultural_context', '')
        )
        tk_emb = result['embeddings']
        tk_mask = result['attention_mask']
        # Mean pool over sequence length
        tk_emb_pooled = (tk_emb * tk_mask.unsqueeze(-1)).sum(dim=0) / tk_mask.sum()
        tk_embeddings_list.append(tk_emb_pooled)
    
    tk_embeddings = torch.stack(tk_embeddings_list)  # [30, embedding_dim]
    
    # Pre-compute all genomic target embeddings
    print("Encoding genomic targets...")
    genomic_embeddings_list = []
    
    for target in tqdm(dataset.genomic_targets):
        # Get all pathways, tissues, diseases
        all_pathways = []
        all_tissues = []
        all_diseases = []
        for cond in target.get('associated_conditions', []):
            all_pathways.extend(cond.get('pathways', []))
            all_tissues.extend(cond.get('tissues', []))
            all_diseases.append(cond.get('condition', ''))
        
        result = tokenizer.encode_genomic_target(
            gene_id=target['gene_id'],
            pathways=all_pathways,
            tissues=all_tissues,
            diseases=all_diseases
        )
        genomic_emb = result['embeddings']
        # Already pooled to [1, embedding_dim * 4]
        genomic_embeddings_list.append(genomic_emb.squeeze(0))
    
    genomic_embeddings = torch.stack(genomic_embeddings_list)  # [46, embedding_dim * 4]
    
    # Get embedding dimension
    tk_dim = tk_embeddings.size(1)
    genomic_dim = genomic_embeddings.size(1)
    
    print(f"\nTK embeddings shape: {tk_embeddings.shape} ({tk_dim}-d)")
    print(f"Genomic embeddings shape: {genomic_embeddings.shape} ({genomic_dim}-d)")
    
    # For fair comparison, we need same dimensionality
    # Option 1: Project genomic down to tk_dim
    # Option 2: Use first tk_dim dimensions of genomic
    # We'll use option 2 for speed
    if genomic_dim > tk_dim:
        genomic_embeddings_projected = genomic_embeddings[:, :tk_dim]
        print(f"Projected genomic to {tk_dim}-d for comparison")
    elif genomic_dim < tk_dim:
        # Pad genomic with zeros
        padding = torch.zeros(genomic_embeddings.size(0), tk_dim - genomic_dim)
        genomic_embeddings_projected = torch.cat([genomic_embeddings, padding], dim=1)
        print(f"Padded genomic to {tk_dim}-d for comparison")
    else:
        genomic_embeddings_projected = genomic_embeddings
        print("Dimensions match, no projection needed")
    
    # Evaluate on test set
    print("\n" + "=" * 70)
    print("Evaluating on Test Set")
    print("=" * 70)
    
    tk_acc_at_1 = []
    tk_acc_at_3 = []
    tk_acc_at_5 = []
    tk_acc_at_10 = []
    g_acc_at_1 = []
    g_acc_at_3 = []
    tk_mrr = []
    g_mrr = []
    
    for idx in tqdm(range(len(test_ds)), desc="Processing test samples"):
        sample = test_ds[idx]
        
        # Get sample embeddings
        sample_tk_emb = sample['tk_embeddings']  # [seq_len, embedding_dim]
        sample_tk_mask = sample['tk_mask']
        sample_genomic_emb = sample['genomic_embeddings']  # [1, embedding_dim * 4]
        
        # Mean pool TK embedding
        sample_tk_pooled = (sample_tk_emb * sample_tk_mask.unsqueeze(-1)).sum(dim=0) / sample_tk_mask.sum()
        sample_tk_pooled = sample_tk_pooled.unsqueeze(0)  # [1, embedding_dim]
        
        # Project genomic embedding to match TK dimension
        sample_genomic_dim = sample_genomic_emb.size(1)
        if sample_genomic_dim > tk_dim:
            sample_genomic_projected = sample_genomic_emb[:, :tk_dim]
        elif sample_genomic_dim < tk_dim:
            padding = torch.zeros(1, tk_dim - sample_genomic_dim)
            sample_genomic_projected = torch.cat([sample_genomic_emb, padding], dim=1)
        else:
            sample_genomic_projected = sample_genomic_emb
        
        # TK → Genomic: compare sample TK with all genomic targets
        tk_to_g_similarities = cosine_similarity_matrix(
            sample_tk_pooled, 
            genomic_embeddings_projected
        )  # [1, 46]
        
        tk_targets = sample['tk_to_genomic_target'].unsqueeze(0)  # [1, 46]
        tk_acc_at_1.append(accuracy_at_k(tk_to_g_similarities, tk_targets, k=1))
        tk_acc_at_3.append(accuracy_at_k(tk_to_g_similarities, tk_targets, k=3))
        tk_acc_at_5.append(accuracy_at_k(tk_to_g_similarities, tk_targets, k=5))
        tk_acc_at_10.append(accuracy_at_k(tk_to_g_similarities, tk_targets, k=10))
        tk_mrr.append(mean_reciprocal_rank(tk_to_g_similarities, tk_targets))
        
        # Genomic → TK: compare sample genomic with all TK practices
        g_to_tk_similarities = cosine_similarity_matrix(
            sample_genomic_projected,
            tk_embeddings
        )  # [1, 30]
        
        g_targets = sample['genomic_to_tk_target'].unsqueeze(0)  # [1, 30]
        g_acc_at_1.append(accuracy_at_k(g_to_tk_similarities, g_targets, k=1))
        g_acc_at_3.append(accuracy_at_k(g_to_tk_similarities, g_targets, k=3))
        g_mrr.append(mean_reciprocal_rank(g_to_tk_similarities, g_targets))
    
    # Results
    results = {
        'tk_acc@1': np.mean(tk_acc_at_1),
        'tk_acc@3': np.mean(tk_acc_at_3),
        'tk_acc@5': np.mean(tk_acc_at_5),
        'tk_acc@10': np.mean(tk_acc_at_10),
        'g_acc@1': np.mean(g_acc_at_1),
        'g_acc@3': np.mean(g_acc_at_3),
        'tk_mrr': np.mean(tk_mrr),
        'g_mrr': np.mean(g_mrr),
    }
    
    print("\n📊 Baseline Results:")
    print(f"\n  TK→Genomic Metrics:")
    print(f"    Acc@1:  {results['tk_acc@1']:.2%}")
    print(f"    Acc@3:  {results['tk_acc@3']:.2%}")
    print(f"    Acc@5:  {results['tk_acc@5']:.2%}")
    print(f"    Acc@10: {results['tk_acc@10']:.2%}")
    print(f"    MRR:    {results['tk_mrr']:.4f}")
    print(f"\n  Genomic→TK Metrics:")
    print(f"    Acc@1:  {results['g_acc@1']:.2%}")
    print(f"    Acc@3:  {results['g_acc@3']:.2%}")
    print(f"    MRR:    {results['g_mrr']:.4f}")
    
    # Compare with transformer
    print("\n" + "=" * 70)
    print("Comparison with 30M Transformer")
    print("=" * 70)
    
    transformer_results = {
        'tk_acc@5': 0.1347,
        'g_acc@3': 0.0632,
    }
    
    print(f"\n  TK→Genomic Acc@5:")
    print(f"    Baseline:    {results['tk_acc@5']:.2%}")
    print(f"    Transformer: {transformer_results['tk_acc@5']:.2%}")
    if results['tk_acc@5'] > transformer_results['tk_acc@5']:
        print(f"    → Baseline WINS by {(results['tk_acc@5'] - transformer_results['tk_acc@5']):.2%}! 🎉")
    else:
        print(f"    → Transformer wins by {(transformer_results['tk_acc@5'] - results['tk_acc@5']):.2%}")
    
    print(f"\n  Genomic→TK Acc@3:")
    print(f"    Baseline:    {results['g_acc@3']:.2%}")
    print(f"    Transformer: {transformer_results['g_acc@3']:.2%}")
    if results['g_acc@3'] > transformer_results['g_acc@3']:
        print(f"    → Baseline WINS by {(results['g_acc@3'] - transformer_results['g_acc@3']):.2%}! 🎉")
    else:
        print(f"    → Transformer wins by {(transformer_results['g_acc@3'] - results['g_acc@3']):.2%}")
    
    # Save results
    output_path = Path('models/genomepath/baseline_results.json')
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, 'w') as f:
        json.dump(results, f, indent=2)
    
    print(f"\n✅ Results saved to {output_path}")
    
    # Interpretation
    print("\n" + "=" * 70)
    print("Interpretation")
    print("=" * 70)
    
    if results['tk_acc@5'] > transformer_results['tk_acc@5']:
        print("\n⚠️  Simple cosine similarity beats the 30M transformer!")
        print("    → The transformer is overfitting and not learning useful patterns")
        print("    → Recommendation: Use simpler model or improve data quality")
    else:
        print("\n✅ Transformer is learning something beyond simple similarity")
        print("    → But accuracy is still too low for production")
        print("    → Recommendation: Continue training or try different architecture")
    
    # Random baseline
    random_tk_acc5 = 5 / 46  # 5 guesses out of 46 targets
    random_g_acc3 = 3 / 30   # 3 guesses out of 30 practices
    
    print(f"\n📊 Random Baseline:")
    print(f"    TK→G Acc@5: {random_tk_acc5:.2%}")
    print(f"    G→TK Acc@3: {random_g_acc3:.2%}")
    
    if results['tk_acc@5'] < random_tk_acc5 * 1.5:
        print(f"\n🔴 WARNING: Baseline is barely better than random!")
        print(f"    → Data quality issue or embeddings not capturing relationships")


if __name__ == "__main__":
    main()
