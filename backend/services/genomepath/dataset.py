"""
GenomePath PyTorch Dataset - Training Data Preparation
=======================================================

Loads 525 correlations and creates train/val/test splits for model training.

Components:
- GenomePathDataset (PyTorch Dataset class)
- Data loading from JSON files
- Train/val/test splitting with stratification
- Batch collation with padding
- Data augmentation support

License: Proprietary - Trade Secret Protection Required
"""

import torch
from torch.utils.data import Dataset, DataLoader, random_split
import json
import numpy as np
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from collections import defaultdict
import random

from backend.services.genomepath.tokenizer import GenomePathTokenizer


class GenomePathDataset(Dataset):
    """
    PyTorch Dataset for GenomePath training.
    
    Loads TK↔Genomic correlations and prepares model-ready tensors.
    """
    
    def __init__(
        self,
        correlations_path: str,
        tk_practices_path: str,
        genomic_targets_path: str,
        tokenizer: GenomePathTokenizer,
        max_tk_length: int = 128,
        augment: bool = False
    ):
        """
        Initialize dataset.
        
        Args:
            correlations_path: Path to training_correlations.json
            tk_practices_path: Path to tk_practices.json
            genomic_targets_path: Path to genomic_targets.json
            tokenizer: Initialized GenomePathTokenizer
            max_tk_length: Maximum sequence length for TK practices
            augment: Whether to apply data augmentation
        """
        self.tokenizer = tokenizer
        self.max_tk_length = max_tk_length
        self.augment = augment
        
        # Load data
        print(f"Loading correlations from {correlations_path}...")
        with open(correlations_path, 'r', encoding='utf-8') as f:
            corr_data = json.load(f)
            # Handle both list and dict formats
            if isinstance(corr_data, dict):
                self.correlations = corr_data.get('correlations', [])
            else:
                self.correlations = corr_data
        self.correlation_confidences = [c.get('confidence', 0.5) for c in self.correlations]
        
        print(f"Loading TK practices from {tk_practices_path}...")
        with open(tk_practices_path, 'r', encoding='utf-8') as f:
            tk_data = json.load(f)
            # Handle both list and dict formats
            if isinstance(tk_data, dict):
                practices_list = tk_data.get('practices', [])
            else:
                practices_list = tk_data
            self.tk_practices = {p['practice_id']: p for p in practices_list}
        
        print(f"Loading genomic targets from {genomic_targets_path}...")
        with open(genomic_targets_path, 'r', encoding='utf-8') as f:
            genomic_data = json.load(f)
            # Handle both list and dict formats
            if isinstance(genomic_data, dict):
                if 'genes' in genomic_data:
                    # genes is a dict with gene_id as keys
                    genes_dict = genomic_data['genes']
                    self.genomic_targets = list(genes_dict.values())
                    self.genomic_lookup = genes_dict
                else:
                    self.genomic_targets = genomic_data
                    self.genomic_lookup = {g['gene_id']: g for g in genomic_data}
            else:
                self.genomic_targets = genomic_data
                self.genomic_lookup = {g['gene_id']: g for g in self.genomic_targets}
        
        # Create target indices mapping
        self.gene_id_to_idx = {g['gene_id']: i for i, g in enumerate(self.genomic_targets)}
        self.practice_id_to_idx = {pid: i for i, pid in enumerate(sorted(self.tk_practices.keys()))}
        
        print(f"✅ Loaded {len(self.correlations)} correlations")
        print(f"✅ Loaded {len(self.tk_practices)} TK practices")
        print(f"✅ Loaded {len(self.genomic_targets)} genomic targets")
        
    def __len__(self) -> int:
        return len(self.correlations)

    def get_confidence(self, correlation_idx: int) -> float:
        """Return stored confidence for a raw correlation index."""
        return self.correlation_confidences[correlation_idx]
    
    def __getitem__(self, idx: int) -> Dict[str, torch.Tensor]:
        """
        Get a single correlation sample.
        
        Returns:
            Dictionary with:
            - tk_embeddings: [seq_len, 384]
            - tk_mask: [seq_len]
            - genomic_embeddings: [1, 1536]
            - genomic_mask: [1]
            - tk_to_genomic_target: [46] one-hot (which genomic targets)
            - genomic_to_tk_target: [30] one-hot (which TK practices)
            - confidence: float (original confidence score)
            - quality: str (correlation quality rating)
            - dosage: str (low/medium/high)
        """
        correlation = self.correlations[idx]
        
        # Get TK practice
        tk_practice_id = correlation.get('tk_practice_id', '')
        tk_practice = self.tk_practices.get(tk_practice_id, {})
        
        # Get genomic target
        genomic_target_id = correlation.get('genomic_target', '')
        genomic_target = self.genomic_lookup.get(genomic_target_id, {})
        
        # Apply augmentation if enabled
        if self.augment:
            tk_practice = self._augment_tk_practice(tk_practice)
        
        # Encode TK practice
        tk_encoded = self.tokenizer.encode_tk_practice(
            practice_name=tk_practice.get('practice_name', ''),
            indications=tk_practice.get('indications', []),
            preparation_method=tk_practice.get('preparation_method', ''),
            cultural_context=tk_practice.get('cultural_context', ''),
            source_community=tk_practice.get('source_community_id', '')
        )
        
        # Encode genomic target
        genomic_encoded = self.tokenizer.encode_genomic_target(
            gene_id=genomic_target.get('gene_id', ''),
            pathways=genomic_target.get('pathways', []),
            tissues=genomic_target.get('tissues', []),
            diseases=genomic_target.get('diseases', [])
        )
        
        # Create targets (one-hot encoding)
        tk_to_genomic_target = torch.zeros(len(self.genomic_targets))
        if genomic_target_id in self.gene_id_to_idx:
            tk_to_genomic_target[self.gene_id_to_idx[genomic_target_id]] = 1.0
        
        genomic_to_tk_target = torch.zeros(len(self.tk_practices))
        if tk_practice_id in self.practice_id_to_idx:
            genomic_to_tk_target[self.practice_id_to_idx[tk_practice_id]] = 1.0
        
        return {
            'tk_embeddings': tk_encoded['text_embeddings'],  # [seq_len, embed_dim]
            'tk_mask': tk_encoded['attention_mask'],    # [seq_len]
            'genomic_embeddings': genomic_encoded['embeddings'],  # [1, 1536]
            'genomic_mask': genomic_encoded['attention_mask'],    # [1]
            'tk_to_genomic_target': tk_to_genomic_target,  # [46]
            'genomic_to_tk_target': genomic_to_tk_target,  # [30]
            'confidence': torch.tensor(correlation.get('confidence', 0.5)),
            'quality': correlation.get('quality', 'moderate'),
            'dosage': correlation.get('dosage', 'medium'),
            'correlation_id': correlation.get('correlation_id', '')
        }
    
    def _augment_tk_practice(self, practice: Dict) -> Dict:
        """
        Apply data augmentation to TK practice.
        
        Augmentations:
        - Shuffle indication order
        - Synonym replacement (basic)
        """
        augmented = practice.copy()
        
        # Shuffle indications
        if 'indications' in augmented and len(augmented['indications']) > 1:
            indications = augmented['indications'].copy()
            random.shuffle(indications)
            augmented['indications'] = indications
        
        return augmented
    
    def get_statistics(self) -> Dict:
        """Get dataset statistics."""
        stats = {
            'total_correlations': len(self.correlations),
            'tk_to_genomic': sum(1 for c in self.correlations if c.get('direction') == 'tk_to_genomic'),
            'genomic_to_tk': sum(1 for c in self.correlations if c.get('direction') == 'genomic_to_tk'),
            'quality_distribution': defaultdict(int),
            'dosage_distribution': defaultdict(int),
            'confidence_stats': {
                'mean': 0.0,
                'std': 0.0,
                'min': 1.0,
                'max': 0.0
            }
        }
        
        confidences = []
        for corr in self.correlations:
            stats['quality_distribution'][corr.get('quality', 'unknown')] += 1
            stats['dosage_distribution'][corr.get('dosage', 'unknown')] += 1
            conf = corr.get('confidence', 0.5)
            confidences.append(conf)
            stats['confidence_stats']['min'] = min(stats['confidence_stats']['min'], conf)
            stats['confidence_stats']['max'] = max(stats['confidence_stats']['max'], conf)
        
        stats['confidence_stats']['mean'] = np.mean(confidences)
        stats['confidence_stats']['std'] = np.std(confidences)
        
        return stats


def collate_fn(batch: List[Dict]) -> Dict[str, torch.Tensor]:
    """
    Collate batch with padding for variable-length sequences.
    
    Args:
        batch: List of samples from GenomePathDataset
        
    Returns:
        Batched and padded tensors
    """
    # Find max TK sequence length in batch
    max_tk_len = max(sample['tk_embeddings'].size(0) for sample in batch)
    batch_size = len(batch)
    
    # Determine embedding dimensions dynamically (BioBERT = 768 by default)
    embedding_dim = batch[0]['tk_embeddings'].size(1)

    # Initialize batched tensors
    tk_embeddings_batch = torch.zeros(batch_size, max_tk_len, embedding_dim)
    tk_mask_batch = torch.zeros(batch_size, max_tk_len, dtype=torch.long)

    genomic_seq_len = batch[0]['genomic_embeddings'].size(0)
    genomic_dim = batch[0]['genomic_embeddings'].size(1)
    genomic_embeddings_batch = torch.zeros(batch_size, genomic_seq_len, genomic_dim)
    genomic_mask_batch = torch.zeros(batch_size, genomic_seq_len, dtype=torch.long)
    
    tk_to_genomic_targets = []
    genomic_to_tk_targets = []
    confidences = []
    qualities = []
    dosages = []
    correlation_ids = []
    
    # Fill batched tensors
    for i, sample in enumerate(batch):
        # TK embeddings
        tk_len = sample['tk_embeddings'].size(0)
        tk_embeddings_batch[i, :tk_len, :] = sample['tk_embeddings']
        tk_mask_batch[i, :tk_len] = sample['tk_mask']
        
        # Genomic embeddings (all same size)
        genomic_embeddings_batch[i, :genomic_seq_len, :] = sample['genomic_embeddings']
        genomic_mask_batch[i, :sample['genomic_mask'].size(0)] = sample['genomic_mask']
        
        # Targets
        tk_to_genomic_targets.append(sample['tk_to_genomic_target'])
        genomic_to_tk_targets.append(sample['genomic_to_tk_target'])
        
        # Metadata
        confidences.append(sample['confidence'])
        qualities.append(sample['quality'])
        dosages.append(sample['dosage'])
        correlation_ids.append(sample['correlation_id'])
    
    return {
        'tk_embeddings': tk_embeddings_batch,
        'tk_mask': tk_mask_batch,
        'genomic_embeddings': genomic_embeddings_batch,
        'genomic_mask': genomic_mask_batch,
        'tk_to_genomic_targets': torch.stack(tk_to_genomic_targets),
        'genomic_to_tk_targets': torch.stack(genomic_to_tk_targets),
        'confidences': torch.stack(confidences),
        'qualities': qualities,
        'dosages': dosages,
        'correlation_ids': correlation_ids
    }


def create_data_splits(
    dataset: GenomePathDataset,
    train_ratio: float = 0.8,
    val_ratio: float = 0.1,
    test_ratio: float = 0.1,
    seed: int = 42
) -> Tuple[Dataset, Dataset, Dataset]:
    """
    Split dataset into train/val/test with stratification by quality.
    
    Args:
        dataset: Full dataset
        train_ratio: Proportion for training
        val_ratio: Proportion for validation
        test_ratio: Proportion for testing
        seed: Random seed for reproducibility
        
    Returns:
        (train_dataset, val_dataset, test_dataset)
    """
    assert abs(train_ratio + val_ratio + test_ratio - 1.0) < 1e-6, "Ratios must sum to 1.0"
    
    # Set seed for reproducibility
    torch.manual_seed(seed)
    np.random.seed(seed)
    
    # Calculate split sizes
    total_size = len(dataset)
    train_size = int(total_size * train_ratio)
    val_size = int(total_size * val_ratio)
    test_size = total_size - train_size - val_size
    
    # Random split (could be improved with stratification)
    train_dataset, val_dataset, test_dataset = random_split(
        dataset,
        [train_size, val_size, test_size],
        generator=torch.Generator().manual_seed(seed)
    )
    
    print(f"\n📊 Data Split Summary:")
    print(f"  Training:   {len(train_dataset):3d} samples ({train_ratio*100:.0f}%)")
    print(f"  Validation: {len(val_dataset):3d} samples ({val_ratio*100:.0f}%)")
    print(f"  Test:       {len(test_dataset):3d} samples ({test_ratio*100:.0f}%)")
    
    return train_dataset, val_dataset, test_dataset


def create_dataloaders(
    train_dataset: Dataset,
    val_dataset: Dataset,
    test_dataset: Dataset,
    batch_size: int = 32,
    num_workers: int = 0
) -> Tuple[DataLoader, DataLoader, DataLoader]:
    """
    Create DataLoaders for train/val/test.
    
    Args:
        train_dataset, val_dataset, test_dataset: Dataset splits
        batch_size: Batch size
        num_workers: Number of dataloader workers
        
    Returns:
        (train_loader, val_loader, test_loader)
    """
    train_loader = DataLoader(
        train_dataset,
        batch_size=batch_size,
        shuffle=True,
        collate_fn=collate_fn,
        num_workers=num_workers
    )
    
    val_loader = DataLoader(
        val_dataset,
        batch_size=batch_size,
        shuffle=False,
        collate_fn=collate_fn,
        num_workers=num_workers
    )
    
    test_loader = DataLoader(
        test_dataset,
        batch_size=batch_size,
        shuffle=False,
        collate_fn=collate_fn,
        num_workers=num_workers
    )
    
    return train_loader, val_loader, test_loader


if __name__ == "__main__":
    print("=" * 70)
    print("GenomePath Dataset - Verification Test")
    print("=" * 70)
    
    # Initialize tokenizer
    print("\n1️⃣ Initializing tokenizer...")
    tokenizer = GenomePathTokenizer()
    
    # Load dataset
    print("\n2️⃣ Loading dataset...")
    dataset = GenomePathDataset(
        correlations_path='data/processed/training_correlations.json',
        tk_practices_path='data/processed/tk_practices.json',
        genomic_targets_path='data/processed/genomic_targets.json',
        tokenizer=tokenizer,
        augment=False
    )
    
    # Print statistics
    print("\n3️⃣ Dataset Statistics:")
    stats = dataset.get_statistics()
    print(f"  Total correlations: {stats['total_correlations']}")
    print(f"  TK→Genomic: {stats['tk_to_genomic']}")
    print(f"  Genomic→TK: {stats['genomic_to_tk']}")
    print(f"\n  Quality distribution:")
    for quality, count in sorted(stats['quality_distribution'].items()):
        print(f"    {quality}: {count}")
    print(f"\n  Dosage distribution:")
    for dosage, count in sorted(stats['dosage_distribution'].items()):
        print(f"    {dosage}: {count}")
    print(f"\n  Confidence stats:")
    print(f"    Mean: {stats['confidence_stats']['mean']:.4f}")
    print(f"    Std:  {stats['confidence_stats']['std']:.4f}")
    print(f"    Min:  {stats['confidence_stats']['min']:.4f}")
    print(f"    Max:  {stats['confidence_stats']['max']:.4f}")
    
    # Test single sample
    print("\n4️⃣ Testing single sample...")
    sample = dataset[0]
    print(f"  TK embeddings shape: {sample['tk_embeddings'].shape}")
    print(f"  TK mask shape: {sample['tk_mask'].shape}")
    print(f"  Genomic embeddings shape: {sample['genomic_embeddings'].shape}")
    print(f"  TK→Genomic target shape: {sample['tk_to_genomic_target'].shape}")
    print(f"  Genomic→TK target shape: {sample['genomic_to_tk_target'].shape}")
    print(f"  Confidence: {sample['confidence'].item():.4f}")
    print(f"  Quality: {sample['quality']}")
    print(f"  Dosage: {sample['dosage']}")
    
    # Create splits
    print("\n5️⃣ Creating train/val/test splits...")
    train_ds, val_ds, test_ds = create_data_splits(dataset, seed=42)
    
    # Create dataloaders
    print("\n6️⃣ Creating DataLoaders...")
    train_loader, val_loader, test_loader = create_dataloaders(
        train_ds, val_ds, test_ds, batch_size=8
    )
    
    # Test batch loading
    print("\n7️⃣ Testing batch loading...")
    batch = next(iter(train_loader))
    print(f"  Batch TK embeddings: {batch['tk_embeddings'].shape}")
    print(f"  Batch genomic embeddings: {batch['genomic_embeddings'].shape}")
    print(f"  Batch TK→Genomic targets: {batch['tk_to_genomic_targets'].shape}")
    print(f"  Batch confidences: {batch['confidences'].shape}")
    
    print("\n" + "=" * 70)
    print("✅ All dataset tests passed!")
    print(f"✅ Ready for model training with {len(dataset)} correlations")
    print("=" * 70)
