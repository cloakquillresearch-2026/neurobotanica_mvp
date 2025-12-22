"""
GenomePath Model Training Script
=================================

Trains the GenomePath transformer model on 525 TK↔Genomic correlations.

Training configuration:
- 420 training samples (80%)
- 52 validation samples (10%)
- 53 test samples (10%)
- Batch size: 32
- Epochs: 50 (with early stopping)
- Optimizer: AdamW (lr=1e-4, weight_decay=0.01)
- Loss: Weighted binary cross-entropy (bidirectional)

Target metrics:
- TK→Genomic Accuracy@5: ≥70%
- Genomic→TK Accuracy@3: ≥60%
- Bidirectional Consistency: ≥75%
- Inference Speed: <200ms

Trade Secret: TS-GP-001 ($6.2B value)
"""

import sys
from pathlib import Path

# Add project root to Python path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.tensorboard import SummaryWriter
import json
from datetime import datetime
from tqdm import tqdm
import numpy as np

from backend.services.genomepath.model import GenomePathModel, ModelConfig
from backend.services.genomepath.tokenizer import GenomePathTokenizer
from backend.services.genomepath.dataset import (
    GenomePathDataset,
    create_data_splits,
    create_dataloaders
)


class BidirectionalLoss(nn.Module):
    """
    Bidirectional loss for TK↔Genomic prediction.
    
    Combines TK→Genomic and Genomic→TK losses with confidence weighting.
    """
    
    def __init__(self, tk_to_genomic_weight=0.6, genomic_to_tk_weight=0.4):
        super().__init__()
        self.tk_weight = tk_to_genomic_weight
        self.genomic_weight = genomic_to_tk_weight
        self.bce = nn.BCELoss(reduction='none')  # Per-sample loss
        
    def forward(self, outputs, targets, confidences):
        """
        Calculate bidirectional loss.
        
        Args:
            outputs: Model outputs dict
            targets: Dict with 'tk_to_genomic' and 'genomic_to_tk' targets
            confidences: Confidence scores for weighting
            
        Returns:
            Total loss, TK→G loss, G→TK loss
        """
        # TK → Genomic loss
        tk_to_g_loss = self.bce(
            outputs['tk_to_genomic_scores'],
            targets['tk_to_genomic']
        )  # [batch, 46]
        
        # Weight by confidence (higher confidence = higher weight)
        confidence_weights = confidences.unsqueeze(1)  # [batch, 1]
        tk_to_g_loss = (tk_to_g_loss * confidence_weights).mean()
        
        # Genomic → TK loss
        genomic_to_tk_loss = self.bce(
            outputs['genomic_to_tk_scores'],
            targets['genomic_to_tk']
        )  # [batch, 30]
        genomic_to_tk_loss = (genomic_to_tk_loss * confidence_weights).mean()
        
        # Combined loss
        total_loss = (
            self.tk_weight * tk_to_g_loss +
            self.genomic_weight * genomic_to_tk_loss
        )
        
        return total_loss, tk_to_g_loss, genomic_to_tk_loss


def accuracy_at_k(predictions, targets, k=5):
    """
    Calculate accuracy@k metric.
    
    Args:
        predictions: [batch, num_classes] confidence scores
        targets: [batch, num_classes] one-hot targets
        k: Number of top predictions to consider
        
    Returns:
        Accuracy@k as float
    """
    batch_size = predictions.size(0)
    
    # Get top-k predictions
    _, top_k_indices = torch.topk(predictions, k=min(k, predictions.size(1)), dim=1)
    
    # Check if true target is in top-k
    correct = 0
    for i in range(batch_size):
        true_indices = torch.where(targets[i] > 0.5)[0]
        if len(true_indices) > 0:
            # Check if any true target is in top-k
            if torch.any(torch.isin(true_indices, top_k_indices[i])):
                correct += 1
    
    return correct / batch_size if batch_size > 0 else 0.0


def mean_reciprocal_rank(predictions, targets):
    """
    Calculate Mean Reciprocal Rank (MRR).
    
    Args:
        predictions: [batch, num_classes] confidence scores
        targets: [batch, num_classes] one-hot targets
        
    Returns:
        MRR as float
    """
    batch_size = predictions.size(0)
    reciprocal_ranks = []
    
    for i in range(batch_size):
        # Get sorted indices (highest confidence first)
        sorted_indices = torch.argsort(predictions[i], descending=True)
        
        # Find true targets
        true_indices = torch.where(targets[i] > 0.5)[0]
        
        if len(true_indices) > 0:
            # Find rank of first true target
            for rank, idx in enumerate(sorted_indices, start=1):
                if idx in true_indices:
                    reciprocal_ranks.append(1.0 / rank)
                    break
        else:
            reciprocal_ranks.append(0.0)
    
    return np.mean(reciprocal_ranks) if reciprocal_ranks else 0.0


def bidirectional_consistency(tk_to_g_scores, g_to_tk_scores, tk_targets, g_targets, threshold=0.75):
    """
    Calculate bidirectional consistency.
    
    Checks if TK→Genomic and Genomic→TK predictions agree.
    
    Args:
        tk_to_g_scores: [batch, 46] TK→Genomic predictions
        g_to_tk_scores: [batch, 30] Genomic→TK predictions
        tk_targets, g_targets: Ground truth targets
        threshold: Confidence threshold for consistency
        
    Returns:
        Consistency ratio (0-1)
    """
    batch_size = tk_to_g_scores.size(0)
    consistent = 0
    
    for i in range(batch_size):
        # Get predicted genomic target (highest confidence)
        pred_genomic_idx = torch.argmax(tk_to_g_scores[i])
        pred_genomic_conf = tk_to_g_scores[i, pred_genomic_idx]
        
        # Get predicted TK practice
        pred_tk_idx = torch.argmax(g_to_tk_scores[i])
        pred_tk_conf = g_to_tk_scores[i, pred_tk_idx]
        
        # Check if both directions agree with ground truth
        true_genomic_idx = torch.argmax(tk_targets[i])
        true_tk_idx = torch.argmax(g_targets[i])
        
        if (pred_genomic_conf >= threshold and pred_tk_conf >= threshold and
            pred_genomic_idx == true_genomic_idx and pred_tk_idx == true_tk_idx):
            consistent += 1
    
    return consistent / batch_size if batch_size > 0 else 0.0


def train_epoch(model, dataloader, criterion, optimizer, device, epoch):
    """Train for one epoch."""
    model.train()
    total_loss = 0.0
    tk_to_g_loss = 0.0
    g_to_tk_loss = 0.0
    
    progress_bar = tqdm(dataloader, desc=f"Epoch {epoch}")
    
    for batch in progress_bar:
        # Move to device
        tk_emb = batch['tk_embeddings'].to(device)
        tk_mask = batch['tk_mask'].to(device)
        gen_emb = batch['genomic_embeddings'].to(device)
        gen_mask = batch['genomic_mask'].to(device)
        tk_targets = batch['tk_to_genomic_targets'].to(device)
        g_targets = batch['genomic_to_tk_targets'].to(device)
        confidences = batch['confidences'].to(device)
        
        # Forward pass
        optimizer.zero_grad()
        outputs = model(tk_emb, gen_emb, tk_mask, gen_mask)
        
        # Calculate loss
        loss, tk_loss, g_loss = criterion(
            outputs,
            {'tk_to_genomic': tk_targets, 'genomic_to_tk': g_targets},
            confidences
        )
        
        # Backward pass
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
        optimizer.step()
        
        # Accumulate
        total_loss += loss.item()
        tk_to_g_loss += tk_loss.item()
        g_to_tk_loss += g_loss.item()
        
        progress_bar.set_postfix({
            'loss': f'{loss.item():.4f}',
            'tk→g': f'{tk_loss.item():.4f}',
            'g→tk': f'{g_loss.item():.4f}'
        })
    
    num_batches = len(dataloader)
    return {
        'loss': total_loss / num_batches,
        'tk_to_genomic_loss': tk_to_g_loss / num_batches,
        'genomic_to_tk_loss': g_to_tk_loss / num_batches
    }


def evaluate(model, dataloader, criterion, device):
    """Evaluate model on validation/test set."""
    model.eval()
    total_loss = 0.0
    tk_to_g_loss = 0.0
    g_to_tk_loss = 0.0
    
    # Metrics
    tk_acc_at_1 = []
    tk_acc_at_3 = []
    tk_acc_at_5 = []
    tk_acc_at_10 = []
    g_acc_at_1 = []
    g_acc_at_3 = []
    tk_mrr = []
    g_mrr = []
    consistency = []
    
    with torch.no_grad():
        for batch in dataloader:
            # Move to device
            tk_emb = batch['tk_embeddings'].to(device)
            tk_mask = batch['tk_mask'].to(device)
            gen_emb = batch['genomic_embeddings'].to(device)
            gen_mask = batch['genomic_mask'].to(device)
            tk_targets = batch['tk_to_genomic_targets'].to(device)
            g_targets = batch['genomic_to_tk_targets'].to(device)
            confidences = batch['confidences'].to(device)
            
            # Forward pass
            outputs = model(tk_emb, gen_emb, tk_mask, gen_mask)
            
            # Calculate loss
            loss, tk_loss, g_loss = criterion(
                outputs,
                {'tk_to_genomic': tk_targets, 'genomic_to_tk': g_targets},
                confidences
            )
            
            total_loss += loss.item()
            tk_to_g_loss += tk_loss.item()
            g_to_tk_loss += g_loss.item()
            
            # Calculate metrics
            tk_scores = outputs['tk_to_genomic_scores']
            g_scores = outputs['genomic_to_tk_scores']
            
            tk_acc_at_1.append(accuracy_at_k(tk_scores, tk_targets, k=1))
            tk_acc_at_3.append(accuracy_at_k(tk_scores, tk_targets, k=3))
            tk_acc_at_5.append(accuracy_at_k(tk_scores, tk_targets, k=5))
            tk_acc_at_10.append(accuracy_at_k(tk_scores, tk_targets, k=10))
            g_acc_at_1.append(accuracy_at_k(g_scores, g_targets, k=1))
            g_acc_at_3.append(accuracy_at_k(g_scores, g_targets, k=3))
            tk_mrr.append(mean_reciprocal_rank(tk_scores, tk_targets))
            g_mrr.append(mean_reciprocal_rank(g_scores, g_targets))
            consistency.append(bidirectional_consistency(tk_scores, g_scores, tk_targets, g_targets))
    
    num_batches = len(dataloader)
    return {
        'loss': total_loss / num_batches,
        'tk_to_genomic_loss': tk_to_g_loss / num_batches,
        'genomic_to_tk_loss': g_to_tk_loss / num_batches,
        'tk_acc@1': np.mean(tk_acc_at_1),
        'tk_acc@3': np.mean(tk_acc_at_3),
        'tk_acc@5': np.mean(tk_acc_at_5),
        'tk_acc@10': np.mean(tk_acc_at_10),
        'g_acc@1': np.mean(g_acc_at_1),
        'g_acc@3': np.mean(g_acc_at_3),
        'tk_mrr': np.mean(tk_mrr),
        'g_mrr': np.mean(g_mrr),
        'bidirectional_consistency': np.mean(consistency)
    }


def main(epochs=50, batch_size=32, quick_test=False):
    """
    Main training function.
    
    Args:
        epochs: Number of training epochs
        batch_size: Batch size for training
        quick_test: If True, run 2-epoch POC test
    """
    print("=" * 70)
    print("GenomePath Transformer Training")
    if quick_test:
        print("(Quick Test Mode - 2 epochs)")
    print("=" * 70)
    
    # Configuration
    config = ModelConfig()
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"\nDevice: {device}")
    
    # Hyperparameters
    BATCH_SIZE = batch_size
    EPOCHS = 2 if quick_test else epochs
    LEARNING_RATE = 1e-4
    WEIGHT_DECAY = 0.01
    PATIENCE = 10 if not quick_test else 5
    
    print(f"Batch size: {BATCH_SIZE}")
    print(f"Epochs: {EPOCHS}")
    print(f"Learning rate: {LEARNING_RATE}")
    print(f"Weight decay: {WEIGHT_DECAY}")
    print(f"Early stopping patience: {PATIENCE}")
    
    # Create directories
    checkpoint_dir = Path('models/genomepath')
    checkpoint_dir.mkdir(parents=True, exist_ok=True)
    log_dir = Path('runs/genomepath') / datetime.now().strftime('%Y%m%d_%H%M%S')
    log_dir.mkdir(parents=True, exist_ok=True)
    
    # Initialize tensorboard
    writer = SummaryWriter(log_dir)
    
    # Load data
    print("\n" + "=" * 70)
    print("Loading Data...")
    print("=" * 70)
    
    tokenizer = GenomePathTokenizer()
    dataset = GenomePathDataset(
        correlations_path='data/processed/training_correlations.json',
        tk_practices_path='data/processed/tk_practices.json',
        genomic_targets_path='data/processed/genomic_targets.json',
        tokenizer=tokenizer,
        augment=True  # Enable data augmentation
    )
    
    # Create splits
    train_ds, val_ds, test_ds = create_data_splits(dataset, seed=42)
    train_loader, val_loader, test_loader = create_dataloaders(
        train_ds, val_ds, test_ds, batch_size=BATCH_SIZE
    )
    
    # Initialize model
    print("\n" + "=" * 70)
    print("Initializing Model...")
    print("=" * 70)
    
    model = GenomePathModel(config).to(device)
    print(f"Parameters: {model.count_parameters():,}")
    
    # Loss and optimizer
    criterion = BidirectionalLoss()
    optimizer = optim.AdamW(model.parameters(), lr=LEARNING_RATE, weight_decay=WEIGHT_DECAY)
    scheduler = optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=EPOCHS)
    
    # Training loop
    print("\n" + "=" * 70)
    print("Training...")
    print("=" * 70)
    
    best_val_loss = float('inf')
    patience_counter = 0
    
    for epoch in range(1, EPOCHS + 1):
        print(f"\nEpoch {epoch}/{EPOCHS}")
        print("-" * 70)
        
        # Train
        train_metrics = train_epoch(model, train_loader, criterion, optimizer, device, epoch)
        
        # Validate every epoch
        val_metrics = evaluate(model, val_loader, criterion, device)
        
        # Update learning rate
        scheduler.step()
        
        # Log metrics
        writer.add_scalar('Loss/train', train_metrics['loss'], epoch)
        writer.add_scalar('Loss/val', val_metrics['loss'], epoch)
        writer.add_scalar('Accuracy/tk_acc@5', val_metrics['tk_acc@5'], epoch)
        writer.add_scalar('Accuracy/g_acc@3', val_metrics['g_acc@3'], epoch)
        writer.add_scalar('Metrics/MRR_tk', val_metrics['tk_mrr'], epoch)
        writer.add_scalar('Metrics/MRR_g', val_metrics['g_mrr'], epoch)
        writer.add_scalar('Metrics/bidirectional_consistency', val_metrics['bidirectional_consistency'], epoch)
        writer.add_scalar('LR', optimizer.param_groups[0]['lr'], epoch)
        
        # Print summary
        print(f"\nTrain Loss: {train_metrics['loss']:.4f}")
        print(f"Val Loss:   {val_metrics['loss']:.4f}")
        print(f"TK→G Acc@5: {val_metrics['tk_acc@5']:.2%} (target ≥70%)")
        print(f"G→TK Acc@3: {val_metrics['g_acc@3']:.2%} (target ≥60%)")
        print(f"Consistency: {val_metrics['bidirectional_consistency']:.2%} (target ≥75%)")
        print(f"TK MRR: {val_metrics['tk_mrr']:.4f}")
        print(f"G MRR: {val_metrics['g_mrr']:.4f}")
        
        # Save best model
        if val_metrics['loss'] < best_val_loss:
            best_val_loss = val_metrics['loss']
            patience_counter = 0
            
            checkpoint = {
                'epoch': epoch,
                'model_state_dict': model.state_dict(),
                'optimizer_state_dict': optimizer.state_dict(),
                'val_loss': val_metrics['loss'],
                'val_metrics': val_metrics,
                'config': config
            }
            torch.save(checkpoint, checkpoint_dir / 'best_model.pt')
            print(f"✅ Best model saved (val_loss: {val_metrics['loss']:.4f})")
        else:
            patience_counter += 1
            print(f"Patience: {patience_counter}/{PATIENCE}")
        
        # Early stopping
        if patience_counter >= PATIENCE:
            print(f"\n⚠️ Early stopping triggered after {epoch} epochs")
            break
    
    # Final evaluation on test set
    print("\n" + "=" * 70)
    print("Final Evaluation on Test Set...")
    print("=" * 70)
    
    # Load best model (weights_only=False for compatibility with numpy in checkpoint)
    checkpoint = torch.load(checkpoint_dir / 'best_model.pt', weights_only=False)
    model.load_state_dict(checkpoint['model_state_dict'])
    
    test_metrics = evaluate(model, test_loader, criterion, device)
    
    print(f"\nTest Results:")
    print(f"  Loss: {test_metrics['loss']:.4f}")
    print(f"  TK→Genomic Metrics:")
    print(f"    Acc@1:  {test_metrics['tk_acc@1']:.2%}")
    print(f"    Acc@3:  {test_metrics['tk_acc@3']:.2%}")
    print(f"    Acc@5:  {test_metrics['tk_acc@5']:.2%} {'✅' if test_metrics['tk_acc@5'] >= 0.70 else '❌'} (target ≥70%)")
    print(f"    Acc@10: {test_metrics['tk_acc@10']:.2%}")
    print(f"    MRR:    {test_metrics['tk_mrr']:.4f}")
    print(f"  Genomic→TK Metrics:")
    print(f"    Acc@1:  {test_metrics['g_acc@1']:.2%}")
    print(f"    Acc@3:  {test_metrics['g_acc@3']:.2%} {'✅' if test_metrics['g_acc@3'] >= 0.60 else '❌'} (target ≥60%)")
    print(f"    MRR:    {test_metrics['g_mrr']:.4f}")
    print(f"  Bidirectional Consistency: {test_metrics['bidirectional_consistency']:.2%} {'✅' if test_metrics['bidirectional_consistency'] >= 0.75 else '❌'} (target ≥75%)")
    
    # Save final report
    report = {
        'training_completed': datetime.now().isoformat(),
        'best_epoch': checkpoint['epoch'],
        'test_metrics': {k: float(v) for k, v in test_metrics.items()},
        'hyperparameters': {
            'batch_size': BATCH_SIZE,
            'epochs': EPOCHS,
            'learning_rate': LEARNING_RATE,
            'weight_decay': WEIGHT_DECAY
        },
        'model_config': {
            'tk_embedding_dim': config.tk_embedding_dim,
            'genomic_embedding_dim': config.genomic_embedding_dim,
            'bridge_dim': config.bridge_dim,
            'num_transformer_layers': config.num_transformer_layers,
            'num_attention_heads': config.num_attention_heads,
            'total_parameters': model.count_parameters()
        }
    }
    
    with open(checkpoint_dir / 'training_report.json', 'w') as f:
        json.dump(report, f, indent=2)
    
    writer.close()
    
    print("\n" + "=" * 70)
    print("✅ Training Complete!")
    print(f"✅ Best model saved to: {checkpoint_dir / 'best_model.pt'}")
    print(f"✅ Training report: {checkpoint_dir / 'training_report.json'}")
    print(f"✅ TensorBoard logs: {log_dir}")
    print("=" * 70)


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description='Train GenomePath Model')
    parser.add_argument('--epochs', type=int, default=50, help='Number of epochs')
    parser.add_argument('--batch-size', type=int, default=32, help='Batch size')
    parser.add_argument('--quick-test', action='store_true', help='Run 2-epoch POC test')
    
    args = parser.parse_args()
    
    main(epochs=args.epochs, batch_size=args.batch_size, quick_test=args.quick_test)
