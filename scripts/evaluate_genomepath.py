"""Evaluate trained GenomePath model on test set."""

import sys
from pathlib import Path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

import torch
import json

from backend.services.genomepath.model import GenomePathModel
from backend.services.genomepath.tokenizer import GenomePathTokenizer
from backend.services.genomepath.dataset import GenomePathDataset, create_data_splits, create_dataloaders
from scripts.train_genomepath import BidirectionalLoss, evaluate

# Load model
checkpoint_path = Path('models/genomepath/best_model.pt')
print(f"Loading checkpoint from {checkpoint_path}...")
checkpoint = torch.load(checkpoint_path, weights_only=False)

# Initialize model
config = checkpoint['config']
model = GenomePathModel(config)
model.load_state_dict(checkpoint['model_state_dict'])
model.eval()

print(f"✅ Model loaded from epoch {checkpoint['epoch']}")
print(f"✅ Val loss: {checkpoint['val_loss']:.4f}")

# Load dataset
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
_, _, test_loader = create_dataloaders(train_ds, val_ds, test_ds, batch_size=32)

# Evaluate
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
model = model.to(device)
criterion = BidirectionalLoss()

print("\n" + "=" * 70)
print("Test Set Evaluation")
print("=" * 70)

test_metrics = evaluate(model, test_loader, criterion, device)

print(f"\n📊 Test Results:")
print(f"  Loss: {test_metrics['loss']:.4f}")
print(f"\n  TK→Genomic Metrics:")
print(f"    Acc@1:  {test_metrics['tk_acc@1']:.2%}")
print(f"    Acc@3:  {test_metrics['tk_acc@3']:.2%}")
print(f"    Acc@5:  {test_metrics['tk_acc@5']:.2%} {'✅' if test_metrics['tk_acc@5'] >= 0.70 else '❌'} (target ≥70%)")
print(f"    Acc@10: {test_metrics['tk_acc@10']:.2%}")
print(f"    MRR:    {test_metrics['tk_mrr']:.4f}")
print(f"\n  Genomic→TK Metrics:")
print(f"    Acc@1:  {test_metrics['g_acc@1']:.2%}")
print(f"    Acc@3:  {test_metrics['g_acc@3']:.2%} {'✅' if test_metrics['g_acc@3'] >= 0.60 else '❌'} (target ≥60%)")
print(f"    MRR:    {test_metrics['g_mrr']:.4f}")
print(f"\n  Bidirectional Consistency: {test_metrics['bidirectional_consistency']:.2%} {'✅' if test_metrics['bidirectional_consistency'] >= 0.75 else '❌'} (target ≥75%)")

# Save results
results = {
    'test_metrics': {k: float(v) for k, v in test_metrics.items()},
    'checkpoint_epoch': checkpoint['epoch'],
    'val_loss': float(checkpoint['val_loss'])
}

results_path = Path('models/genomepath/test_results.json')
with open(results_path, 'w') as f:
    json.dump(results, f, indent=2)

print(f"\n✅ Results saved to {results_path}")
