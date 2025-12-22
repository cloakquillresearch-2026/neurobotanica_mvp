"""
Hybrid Model Training - Chemical + Text + Genomic Fusion
=========================================================

Trains the 3-tower hybrid architecture combining:
1. Text tower: BioBERT embeddings for TK practice descriptions
2. Chemical tower: Morgan fingerprints + molecular descriptors
3. Genomic tower: Multi-modal genomic features

Dataset: 525 correlations, 30 TK practices, 46 genomic targets
Model: 7.4M parameters (vs 30M text-only)
"""

import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import Dataset, DataLoader, Subset
from torch.utils.tensorboard import SummaryWriter
import numpy as np
from pathlib import Path
import json
from tqdm import tqdm
import argparse
from datetime import datetime
import sys

# Add backend to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from backend.services.genomepath.hybrid_model import HybridGenomePathModel, HybridModelConfig
from backend.services.genomepath.tokenizer import GenomePathTokenizer
from backend.services.genomepath.chemical_encoder import ChemicalEncoder


class HybridGenomePathDataset(Dataset):
    """Dataset for hybrid model with chemical features."""
    
    def __init__(
        self,
        correlations_path: Path,
        tk_practices_path: Path,
        genomic_targets_path: Path,
        compound_mappings_path: Path,
        tokenizer: GenomePathTokenizer
    ):
        self.tokenizer = tokenizer
        
        # Load data
        print(f"Loading correlations from {correlations_path}...")
        with open(correlations_path, 'r', encoding='utf-8') as f:
            corr_data = json.load(f)
        self.correlations = corr_data['correlations']
        
        print(f"Loading TK practices from {tk_practices_path}...")
        with open(tk_practices_path, 'r', encoding='utf-8') as f:
            tk_data = json.load(f)
        self.tk_practices = {p['practice_id']: p for p in tk_data['practices']}
        
        print(f"Loading genomic targets from {genomic_targets_path}...")
        with open(genomic_targets_path, 'r', encoding='utf-8') as f:
            genomic_data = json.load(f)
        # Handle both dict and list formats
        if isinstance(genomic_data.get('genes'), dict):
            self.genomic_targets = genomic_data['genes']
        elif isinstance(genomic_data.get('genes'), list):
            self.genomic_targets = {g['gene_id']: g for g in genomic_data['genes']}
        else:
            # Fallback: assume top-level keys are gene IDs
            self.genomic_targets = {k: v for k, v in genomic_data.items() if k != 'metadata'}
        
        print(f"Loading compound mappings from {compound_mappings_path}...")
        with open(compound_mappings_path, 'r', encoding='utf-8') as f:
            compound_data = json.load(f)
        self.compound_mappings = compound_data['tk_compound_mappings']['practices']
        
        # Create ID mappings
        self.tk_ids = sorted(list(self.tk_practices.keys()))
        self.genomic_ids = sorted(list(self.genomic_targets.keys()))
        
        self.tk_id_to_idx = {tid: idx for idx, tid in enumerate(self.tk_ids)}
        self.genomic_id_to_idx = {gid: idx for idx, gid in enumerate(self.genomic_ids)}
        
        print(f"✅ Loaded {len(self.correlations)} correlations")
        print(f"   {len(self.tk_ids)} TK practices, {len(self.genomic_ids)} genomic targets")
        print(f"   {len(self.compound_mappings)} practices with compound mappings")
    
    def __len__(self):
        return len(self.correlations)
    
    def __getitem__(self, idx):
        corr = self.correlations[idx]
        
        # Handle both TK→Genomic and Genomic→TK directions
        if corr['direction'] == 'tk_to_genomic':
            tk_id = corr['tk_practice_id']
            genomic_target_name = corr['genomic_target']
        else:  # genomic_to_tk
            # For genomic→TK, use predicted practice or skip if invalid
            tk_predicted = corr.get('tk_practice_predicted', '')
            # Try to match to a real practice ID
            tk_id = None
            for pid in self.tk_ids:
                if pid in tk_predicted or tk_predicted in str(self.tk_practices[pid].get('practice_name', '')):
                    tk_id = pid
                    break
            
            # If no match, use first TK practice as fallback
            if tk_id is None:
                tk_id = self.tk_ids[0]
            
            genomic_target_name = corr['genomic_target']
        
        tk_practice = self.tk_practices[tk_id]
        
        # Get compound mapping (use practice_id as key)
        compounds = []
        if tk_id in self.compound_mappings:
            compounds = self.compound_mappings[tk_id]['compounds']
        
        # Encode TK practice with chemical features
        tk_encoding = self.tokenizer.encode_tk_practice(
            practice_name=tk_practice['practice_name'],
            indications=tk_practice['indications'],
            preparation_method=tk_practice.get('preparation_method', ''),
            cultural_context=tk_practice.get('cultural_context', ''),
            source_community=tk_practice.get('source_community_id', ''),
            active_compounds=compounds
        )
        
        # Get genomic target
        genomic_target_name = genomic_target_name
        
        # Find genomic ID (may need to search by gene name)
        genomic_id = None
        for gid, gdata in self.genomic_targets.items():
            if gdata.get('gene_symbol') == genomic_target_name or gdata.get('gene_id') == genomic_target_name:
                genomic_id = gid
                break
        
        if genomic_id is None:
            # Fallback: use first genomic target (for robustness)
            genomic_id = self.genomic_ids[0]
        
        genomic_target = self.genomic_targets[genomic_id]
        
        # Encode genomic target
        genomic_encoding = self.tokenizer.encode_genomic_target(
            gene_id=genomic_target.get('gene_symbol', genomic_target.get('gene_id', '')),
            pathways=genomic_target.get('pathways', []),
            tissues=genomic_target.get('tissues', []),
            diseases=genomic_target.get('diseases', [])
        )
        
        # Create targets
        tk_idx = self.tk_id_to_idx[tk_id]
        genomic_idx = self.genomic_id_to_idx[genomic_id]
        
        tk_to_genomic_target = torch.zeros(len(self.genomic_ids))
        tk_to_genomic_target[genomic_idx] = 1.0
        
        genomic_to_tk_target = torch.zeros(len(self.tk_ids))
        genomic_to_tk_target[tk_idx] = 1.0
        
        # Get confidence and quality
        confidence = torch.tensor(corr.get('confidence', 0.5), dtype=torch.float32)
        
        return {
            # TK inputs
            'tk_text_embeddings': tk_encoding['text_embeddings'],  # [seq_len, 768]
            'tk_attention_mask': tk_encoding['attention_mask'],  # [seq_len]
            'tk_chemical_fingerprint': tk_encoding['chemical_fingerprint'],  # [2048]
            'tk_chemical_descriptors': tk_encoding['chemical_descriptors'],  # [7]
            'tk_chemical_valid': tk_encoding['chemical_valid'],  # [1]
            
            # Genomic inputs
            'genomic_embeddings': genomic_encoding['embeddings'].squeeze(0),  # [3072]
            
            # Targets
            'tk_to_genomic_target': tk_to_genomic_target,  # [46]
            'genomic_to_tk_target': genomic_to_tk_target,  # [30]
            'confidence': confidence
        }


def collate_fn(batch):
    """Collate batch with variable-length sequences."""
    # Find max sequence length
    max_seq_len = max(item['tk_text_embeddings'].size(0) for item in batch)
    
    batch_size = len(batch)
    text_dim = batch[0]['tk_text_embeddings'].size(1)
    genomic_dim = batch[0]['genomic_embeddings'].size(0)
    num_genomic = batch[0]['tk_to_genomic_target'].size(0)
    num_tk = batch[0]['genomic_to_tk_target'].size(0)
    
    # Prepare batched tensors
    tk_text_embeddings = torch.zeros(batch_size, max_seq_len, text_dim)
    tk_attention_mask = torch.zeros(batch_size, max_seq_len)
    tk_chemical_fingerprint = torch.stack([item['tk_chemical_fingerprint'] for item in batch])
    tk_chemical_descriptors = torch.stack([item['tk_chemical_descriptors'] for item in batch])
    tk_chemical_valid = torch.stack([item['tk_chemical_valid'] for item in batch])
    genomic_embeddings = torch.stack([item['genomic_embeddings'] for item in batch])
    tk_to_genomic_target = torch.stack([item['tk_to_genomic_target'] for item in batch])
    genomic_to_tk_target = torch.stack([item['genomic_to_tk_target'] for item in batch])
    confidence = torch.stack([item['confidence'] for item in batch])
    
    # Fill in variable-length sequences
    for i, item in enumerate(batch):
        seq_len = item['tk_text_embeddings'].size(0)
        tk_text_embeddings[i, :seq_len] = item['tk_text_embeddings']
        tk_attention_mask[i, :seq_len] = item['tk_attention_mask']
    
    return {
        'tk_text_embeddings': tk_text_embeddings,
        'tk_attention_mask': tk_attention_mask,
        'tk_chemical_fingerprint': tk_chemical_fingerprint,
        'tk_chemical_descriptors': tk_chemical_descriptors,
        'tk_chemical_valid': tk_chemical_valid,
        'genomic_embeddings': genomic_embeddings,
        'tk_to_genomic_target': tk_to_genomic_target,
        'genomic_to_tk_target': genomic_to_tk_target,
        'confidence': confidence
    }


def train_epoch(model, dataloader, optimizer, device, epoch):
    """Train for one epoch."""
    model.train()
    total_loss = 0
    tk_correct = 0
    g_correct = 0
    total_samples = 0
    
    pbar = tqdm(dataloader, desc=f"Epoch {epoch}")
    
    for batch in pbar:
        # Move to device
        for k in batch:
            if isinstance(batch[k], torch.Tensor):
                batch[k] = batch[k].to(device)
        
        # Forward pass
        outputs = model(
            tk_text_embeddings=batch['tk_text_embeddings'],
            tk_attention_mask=batch['tk_attention_mask'],
            tk_chemical_fingerprint=batch['tk_chemical_fingerprint'],
            tk_chemical_descriptors=batch['tk_chemical_descriptors'],
            tk_chemical_valid=batch['tk_chemical_valid'],
            genomic_embeddings=batch['genomic_embeddings']
        )
        
        # Calculate loss (binary cross-entropy with confidence weighting)
        tk_to_g_loss = nn.functional.binary_cross_entropy_with_logits(
            outputs['tk_to_genomic_logits'],
            batch['tk_to_genomic_target'],
            reduction='none'
        )
        tk_to_g_loss = (tk_to_g_loss.mean(dim=1) * batch['confidence']).mean()
        
        g_to_tk_loss = nn.functional.binary_cross_entropy_with_logits(
            outputs['genomic_to_tk_logits'],
            batch['genomic_to_tk_target'],
            reduction='none'
        )
        g_to_tk_loss = (g_to_tk_loss.mean(dim=1) * batch['confidence']).mean()
        
        # Combined loss (60% TK→G, 40% G→TK)
        loss = 0.6 * tk_to_g_loss + 0.4 * g_to_tk_loss
        
        # Backward pass
        optimizer.zero_grad()
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
        optimizer.step()
        
        # Calculate accuracy@5 for TK→G
        _, tk_pred = torch.topk(outputs['tk_to_genomic_logits'], k=5, dim=1)
        tk_targets = batch['tk_to_genomic_target'].argmax(dim=1, keepdim=True)
        tk_correct += (tk_pred == tk_targets).any(dim=1).sum().item()
        
        # Calculate accuracy@3 for G→TK
        _, g_pred = torch.topk(outputs['genomic_to_tk_logits'], k=3, dim=1)
        g_targets = batch['genomic_to_tk_target'].argmax(dim=1, keepdim=True)
        g_correct += (g_pred == g_targets).any(dim=1).sum().item()
        
        total_samples += batch['tk_text_embeddings'].size(0)
        total_loss += loss.item()
        
        # Update progress bar
        pbar.set_postfix({
            'loss': f"{loss.item():.4f}",
            'tk_acc': f"{tk_correct/total_samples*100:.1f}%",
            'g_acc': f"{g_correct/total_samples*100:.1f}%"
        })
    
    avg_loss = total_loss / len(dataloader)
    tk_acc = tk_correct / total_samples
    g_acc = g_correct / total_samples
    
    return avg_loss, tk_acc, g_acc


@torch.no_grad()
def evaluate(model, dataloader, device):
    """Evaluate model."""
    model.eval()
    total_loss = 0
    tk_correct = 0
    g_correct = 0
    total_samples = 0
    
    for batch in tqdm(dataloader, desc="Evaluating"):
        # Move to device
        for k in batch:
            if isinstance(batch[k], torch.Tensor):
                batch[k] = batch[k].to(device)
        
        # Forward pass
        outputs = model(
            tk_text_embeddings=batch['tk_text_embeddings'],
            tk_attention_mask=batch['tk_attention_mask'],
            tk_chemical_fingerprint=batch['tk_chemical_fingerprint'],
            tk_chemical_descriptors=batch['tk_chemical_descriptors'],
            tk_chemical_valid=batch['tk_chemical_valid'],
            genomic_embeddings=batch['genomic_embeddings']
        )
        
        # Calculate loss
        tk_to_g_loss = nn.functional.binary_cross_entropy_with_logits(
            outputs['tk_to_genomic_logits'],
            batch['tk_to_genomic_target'],
            reduction='none'
        )
        tk_to_g_loss = (tk_to_g_loss.mean(dim=1) * batch['confidence']).mean()
        
        g_to_tk_loss = nn.functional.binary_cross_entropy_with_logits(
            outputs['genomic_to_tk_logits'],
            batch['genomic_to_tk_target'],
            reduction='none'
        )
        g_to_tk_loss = (g_to_tk_loss.mean(dim=1) * batch['confidence']).mean()
        
        loss = 0.6 * tk_to_g_loss + 0.4 * g_to_tk_loss
        
        # Calculate accuracy
        _, tk_pred = torch.topk(outputs['tk_to_genomic_logits'], k=5, dim=1)
        tk_targets = batch['tk_to_genomic_target'].argmax(dim=1, keepdim=True)
        tk_correct += (tk_pred == tk_targets).any(dim=1).sum().item()
        
        _, g_pred = torch.topk(outputs['genomic_to_tk_logits'], k=3, dim=1)
        g_targets = batch['genomic_to_tk_target'].argmax(dim=1, keepdim=True)
        g_correct += (g_pred == g_targets).any(dim=1).sum().item()
        
        total_samples += batch['tk_text_embeddings'].size(0)
        total_loss += loss.item()
    
    avg_loss = total_loss / len(dataloader)
    tk_acc = tk_correct / total_samples
    g_acc = g_correct / total_samples
    
    return avg_loss, tk_acc, g_acc


def main():
    parser = argparse.ArgumentParser(description='Train Hybrid GenomePath Model')
    parser.add_argument('--epochs', type=int, default=20, help='Number of epochs')
    parser.add_argument('--batch-size', type=int, default=16, help='Batch size')
    parser.add_argument('--lr', type=float, default=1e-4, help='Learning rate')
    parser.add_argument('--device', type=str, default='cuda' if torch.cuda.is_available() else 'cpu')
    parser.add_argument('--quick-test', action='store_true', help='Quick test with 2 epochs')
    args = parser.parse_args()
    
    if args.quick_test:
        args.epochs = 2
        args.batch_size = 8
    
    print("\n" + "=" * 70)
    print("Hybrid GenomePath Model Training")
    print("Chemical + Text + Genomic Fusion")
    print("=" * 70)
    
    # Paths
    base_dir = Path(__file__).parent.parent
    data_dir = base_dir / 'data'
    
    # Initialize tokenizer (includes chemical encoder)
    print("\nInitializing tokenizer with chemical encoder...")
    tokenizer = GenomePathTokenizer()
    
    # Load dataset
    print("\nLoading dataset...")
    dataset = HybridGenomePathDataset(
        correlations_path=data_dir / 'processed' / 'training_correlations.json',
        tk_practices_path=data_dir / 'processed' / 'tk_practices.json',
        genomic_targets_path=data_dir / 'processed' / 'genomic_targets.json',
        compound_mappings_path=data_dir / 'training' / 'tk_compound_mappings.json',
        tokenizer=tokenizer
    )
    
    # Create splits
    total_size = len(dataset)
    train_size = int(0.8 * total_size)
    val_size = int(0.1 * total_size)
    test_size = total_size - train_size - val_size
    
    indices = list(range(total_size))
    np.random.seed(42)
    np.random.shuffle(indices)
    
    train_indices = indices[:train_size]
    val_indices = indices[train_size:train_size+val_size]
    test_indices = indices[train_size+val_size:]
    
    train_dataset = Subset(dataset, train_indices)
    val_dataset = Subset(dataset, val_indices)
    test_dataset = Subset(dataset, test_indices)
    
    print(f"Train: {len(train_dataset)}, Val: {len(val_dataset)}, Test: {len(test_dataset)}")
    
    # Create dataloaders
    train_loader = DataLoader(
        train_dataset,
        batch_size=args.batch_size,
        shuffle=True,
        collate_fn=collate_fn,
        num_workers=0  # Windows compatibility
    )
    val_loader = DataLoader(
        val_dataset,
        batch_size=args.batch_size,
        shuffle=False,
        collate_fn=collate_fn,
        num_workers=0
    )
    test_loader = DataLoader(
        test_dataset,
        batch_size=args.batch_size,
        shuffle=False,
        collate_fn=collate_fn,
        num_workers=0
    )
    
    # Initialize model
    print("\nInitializing hybrid model...")
    config = HybridModelConfig(
        num_genomic_targets=len(dataset.genomic_ids),
        num_tk_practices=len(dataset.tk_ids)
    )
    model = HybridGenomePathModel(config).to(args.device)
    
    print(f"Model parameters: {model.count_parameters():,}")
    print(f"Device: {args.device}")
    
    # Optimizer
    optimizer = optim.AdamW(model.parameters(), lr=args.lr, weight_decay=0.01)
    
    # TensorBoard
    log_dir = base_dir / 'runs' / f'hybrid_{datetime.now().strftime("%Y%m%d_%H%M%S")}'
    writer = SummaryWriter(log_dir)
    
    # Training loop
    print(f"\nTraining for {args.epochs} epochs...")
    best_val_loss = float('inf')
    best_model_path = base_dir / 'models' / 'genomepath' / 'hybrid_best_model.pt'
    best_model_path.parent.mkdir(parents=True, exist_ok=True)
    
    for epoch in range(1, args.epochs + 1):
        print(f"\n{'='*70}")
        print(f"Epoch {epoch}/{args.epochs}")
        print('='*70)
        
        # Train
        train_loss, train_tk_acc, train_g_acc = train_epoch(
            model, train_loader, optimizer, args.device, epoch
        )
        
        # Validate
        val_loss, val_tk_acc, val_g_acc = evaluate(model, val_loader, args.device)
        
        # Log
        writer.add_scalar('Loss/train', train_loss, epoch)
        writer.add_scalar('Loss/val', val_loss, epoch)
        writer.add_scalar('Accuracy/train_tk_to_g', train_tk_acc, epoch)
        writer.add_scalar('Accuracy/val_tk_to_g', val_tk_acc, epoch)
        writer.add_scalar('Accuracy/train_g_to_tk', train_g_acc, epoch)
        writer.add_scalar('Accuracy/val_g_to_tk', val_g_acc, epoch)
        
        print(f"\nResults:")
        print(f"  Train Loss: {train_loss:.4f}, TK→G Acc@5: {train_tk_acc*100:.2f}%, G→TK Acc@3: {train_g_acc*100:.2f}%")
        print(f"  Val Loss:   {val_loss:.4f}, TK→G Acc@5: {val_tk_acc*100:.2f}%, G→TK Acc@3: {val_g_acc*100:.2f}%")
        
        # Save best model
        if val_loss < best_val_loss:
            best_val_loss = val_loss
            torch.save({
                'epoch': epoch,
                'model_state_dict': model.state_dict(),
                'optimizer_state_dict': optimizer.state_dict(),
                'val_loss': val_loss,
                'val_tk_acc': val_tk_acc,
                'val_g_acc': val_g_acc,
                'config': config
            }, best_model_path)
            print(f"  ✅ Saved best model (val_loss: {val_loss:.4f})")
    
    # Final test evaluation
    print(f"\n{'='*70}")
    print("Final Test Evaluation")
    print('='*70)
    
    # Load best model
    checkpoint = torch.load(best_model_path, weights_only=False)
    model.load_state_dict(checkpoint['model_state_dict'])
    
    test_loss, test_tk_acc, test_g_acc = evaluate(model, test_loader, args.device)
    
    print(f"\nTest Results:")
    print(f"  Loss: {test_loss:.4f}")
    print(f"  TK→Genomic Acc@5: {test_tk_acc*100:.2f}%")
    print(f"  Genomic→TK Acc@3: {test_g_acc*100:.2f}%")
    
    # Save results
    results = {
        'test_loss': float(test_loss),
        'test_tk_to_g_acc@5': float(test_tk_acc),
        'test_g_to_tk_acc@3': float(test_g_acc),
        'best_epoch': checkpoint['epoch'],
        'best_val_loss': float(checkpoint['val_loss']),
        'model_parameters': model.count_parameters(),
        'training_date': datetime.now().isoformat()
    }
    
    results_path = base_dir / 'models' / 'genomepath' / 'hybrid_test_results.json'
    with open(results_path, 'w') as f:
        json.dump(results, f, indent=2)
    
    print(f"\n✅ Results saved to {results_path}")
    print(f"✅ Best model saved to {best_model_path}")
    print(f"✅ TensorBoard logs: {log_dir}")
    
    writer.close()
    
    print("\n" + "=" * 70)
    print("Training complete!")
    print("=" * 70)


if __name__ == "__main__":
    main()
