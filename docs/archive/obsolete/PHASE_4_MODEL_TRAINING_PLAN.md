# GenomePath Phase 4: Transformer Model Training
**Status**: 🟢 READY TO BEGIN  
**Prerequisites**: ✅ ALL COMPLETE (525 correlations, 30 TK practices, 46 genomic targets)  
**Timeline**: 3-5 days (AI-accelerated development)  
**Trade Secret**: TS-GP-001 ($6.2B bidirectional semantic bridge)

---

## 🎯 Phase 4 Objectives

Transform the 525-correlation dataset into a production-ready transformer model capable of:

1. **TK → Genomic Predictions**: Given traditional practice, predict relevant genomic targets with confidence scores
2. **Genomic → TK Suggestions**: Given genomic findings, suggest validated traditional practices
3. **Dosage-Specific Recommendations**: Provide low/medium/high dose guidance
4. **Multi-Target Synergy Recognition**: Identify synergistic target combinations
5. **<2-Second Inference**: Meet Nevada dispensary real-time requirements

---

## 📊 Training Dataset Summary

**Inputs** (from Phase 3 completion):
- **525 total correlations** (103% of target)
  - 444 TK → Genomic (84.6%)
  - 81 Genomic → TK (15.4%)
- **30 TK practices** from 30 global traditions
- **46 genomic targets** from 197 NORML studies
- **3 dosage profiles** per practice (low/medium/high)
- **11 synergy pairs** for multi-target recognition

**Quality Metrics**:
- Average confidence: 0.5983
- Quality distribution: 3.4% good, 52.8% moderate, 43.8% poor
- 100% ethical compliance (0 sacred violations)

---

## 🏗️ Architecture Design

### Model Type: **Dual-Tower Transformer**

**Rationale**: Bidirectional correlation requires separate encoders for TK and Genomic spaces, then cross-attention for alignment.

```
┌─────────────────────────────────────────────────────────────┐
│                    GenomePath Transformer                   │
├─────────────────────────────────────────────────────────────┤
│                                                             │
│  ┌──────────────────┐         ┌──────────────────┐        │
│  │   TK Encoder     │         │ Genomic Encoder  │        │
│  │  (512-d output)  │         │  (512-d output)  │        │
│  │                  │         │                  │        │
│  │ - Practice name  │         │ - Gene ID        │        │
│  │ - Indications    │         │ - Pathway        │        │
│  │ - Preparation    │         │ - Tissue         │        │
│  │ - Cultural ctx   │         │ - Disease assoc  │        │
│  └────────┬─────────┘         └────────┬─────────┘        │
│           │                            │                   │
│           └──────────┬─────────────────┘                   │
│                      │                                     │
│            ┌─────────▼──────────┐                         │
│            │  Cross-Attention   │                         │
│            │   Bridge Layer     │                         │
│            │   (256-d shared)   │                         │
│            └─────────┬──────────┘                         │
│                      │                                     │
│         ┌────────────┴────────────┐                       │
│         │                         │                       │
│   ┌─────▼──────┐         ┌───────▼──────┐               │
│   │ TK→Genomic │         │ Genomic→TK   │               │
│   │  Predictor │         │  Predictor   │               │
│   │ (10 targets)│        │ (8 practices)│               │
│   └────────────┘         └──────────────┘               │
│                                                             │
└─────────────────────────────────────────────────────────────┘
```

### Components

**1. TK Encoder** (PyTorch Transformer Encoder):
- **Input**: Practice name + indications + preparation + cultural context
- **Tokenization**: SentenceTransformer embeddings (all-MiniLM-L6-v2)
- **Architecture**: 4 transformer layers, 8 attention heads
- **Output**: 512-dimensional semantic vector

**2. Genomic Encoder** (PyTorch Transformer Encoder):
- **Input**: Gene ID + pathway + tissue + disease associations
- **Tokenization**: Scientific term embeddings
- **Architecture**: 4 transformer layers, 8 attention heads
- **Output**: 512-dimensional semantic vector

**3. Cross-Attention Bridge**:
- **Input**: Concatenated TK + Genomic embeddings (1024-d)
- **Architecture**: Multi-head cross-attention (8 heads)
- **Output**: 256-dimensional shared semantic space
- **Trade Secret**: Alignment matrix learns TK↔Genomic correspondences

**4. Prediction Heads**:

**TK→Genomic Predictor**:
- Input: Bridge output (256-d) + TK embedding (512-d)
- Architecture: 2-layer MLP [768-d → 384-d → 46-d]
- Output: Confidence scores for 46 genomic targets
- Loss: Binary cross-entropy (multi-label classification)

**Genomic→TK Predictor**:
- Input: Bridge output (256-d) + Genomic embedding (512-d)
- Architecture: 2-layer MLP [768-d → 384-d → 30-d]
- Output: Confidence scores for 30 TK practices
- Loss: Binary cross-entropy (multi-label classification)

---

## 🔢 Training Strategy

### Data Preparation

**Train/Val/Test Split**:
- Training: 420 correlations (80%)
- Validation: 53 correlations (10%)
- Test: 52 correlations (10%)

**Stratification**: By correlation quality (good/moderate/poor) to ensure balanced representation

**Data Augmentation**:
1. **Dosage Mixing**: Use all 3 dosage variations per practice
2. **Synonym Expansion**: Vary practice names (e.g., "Cannabis for pain" ↔ "Hemp for analgesia")
3. **Indication Shuffling**: Reorder indications to prevent positional bias
4. **Target Sampling**: Random subset of genomic targets per batch

### Training Configuration

**Hyperparameters**:
```python
BATCH_SIZE = 32
LEARNING_RATE = 1e-4
EPOCHS = 50
OPTIMIZER = AdamW(weight_decay=0.01)
SCHEDULER = CosineAnnealingLR(T_max=50)
DROPOUT = 0.1
ATTENTION_HEADS = 8
TRANSFORMER_LAYERS = 4
EMBEDDING_DIM = 512
BRIDGE_DIM = 256
```

**Loss Function**:
```python
# Weighted multi-label binary cross-entropy
loss_tk2g = BCE(predicted_targets, true_targets, weight=confidence_scores)
loss_g2tk = BCE(predicted_practices, true_practices, weight=confidence_scores)
total_loss = 0.6 * loss_tk2g + 0.4 * loss_g2tk  # TK→G weighted higher
```

**Regularization**:
- Dropout: 0.1 in all transformer layers
- Weight decay: 0.01 in AdamW
- Early stopping: Patience 10 epochs on validation loss
- Gradient clipping: Max norm 1.0

### Training Process

**Epoch Structure**:
1. Shuffle training data
2. For each batch:
   - Forward pass through both encoders
   - Cross-attention bridge alignment
   - Prediction head outputs
   - Calculate bidirectional loss
   - Backpropagate and update weights
3. Validation evaluation every 5 epochs
4. Save best model checkpoint (lowest val loss)
5. Log metrics: loss, accuracy@k, bidirectional consistency

**Convergence Criteria**:
- Validation loss < 0.15 (target)
- TK→Genomic accuracy@5 > 70%
- Genomic→TK accuracy@3 > 60%
- Bidirectional consistency > 75%

---

## 📈 Evaluation Metrics

### Primary Metrics

**1. Accuracy@K** (K=1,3,5,10):
- **TK→Genomic**: % of times true targets appear in top-K predictions
- **Genomic→TK**: % of times true practices appear in top-K predictions
- **Target**: Accuracy@5 ≥ 70% for TK→G, Accuracy@3 ≥ 60% for G→TK

**2. Bidirectional Consistency**:
- For each TK↔Genomic pair, check if both directions predict each other
- **Target**: ≥75% consistency (matches current correlator threshold)

**3. Mean Reciprocal Rank (MRR)**:
- Average of 1/rank for first correct prediction
- Measures ranking quality
- **Target**: MRR ≥ 0.40

**4. Confidence Calibration**:
- Compare predicted confidence vs actual accuracy
- Expected Calibration Error (ECE) < 0.10
- Ensures confidence scores are reliable

### Secondary Metrics

**5. Dosage-Specific Accuracy**:
- Low dose: Accuracy for top 5 targets
- Medium dose: Accuracy for top 8 targets
- High dose: Accuracy for all 10 targets
- **Target**: Low dose ≥75%, Medium ≥65%, High ≥55%

**6. Synergy Recognition**:
- % of predictions containing known synergy pairs
- **Target**: ≥30% of top-5 predictions include synergy pairs

**7. Inference Speed**:
- Time to generate predictions for single practice
- **Target**: <200ms (10x safety margin for 2-second Nevada requirement)

**8. Ethical Compliance**:
- 0 sacred knowledge violations in predictions
- 100% community attribution preserved
- **Target**: Maintain 100% compliance

---

## 🛠️ Implementation Plan

### Day 1: Model Architecture (6-8 hours)

**Tasks**:
1. Install PyTorch + transformers library
2. Create `backend/services/genomepath/model.py`:
   - `TKEncoder` class (transformer-based)
   - `GenomicEncoder` class (transformer-based)
   - `CrossAttentionBridge` class
   - `GenomePathModel` class (dual-tower architecture)
3. Create `backend/services/genomepath/tokenizer.py`:
   - TK practice text preprocessing
   - Genomic feature encoding
   - SentenceTransformer integration
4. Unit tests for model components

**Deliverables**:
- `model.py` (~400 lines)
- `tokenizer.py` (~200 lines)
- `tests/test_genomepath_model.py` (~300 lines)
- Model instantiation verified

### Day 2: Data Pipeline (5-7 hours)

**Tasks**:
1. Create `scripts/prepare_training_data.py`:
   - Load 525 correlations from JSON
   - Create train/val/test splits (80/10/10)
   - Implement data augmentation
   - Export PyTorch datasets
2. Create `backend/services/genomepath/dataset.py`:
   - `GenomePathDataset` class (PyTorch Dataset)
   - Batch collation with padding
   - Dosage-specific sampling
3. Verify data loading with DataLoader
4. Statistics on dataset characteristics

**Deliverables**:
- `prepare_training_data.py` (~300 lines)
- `dataset.py` (~250 lines)
- `data/processed/train_data.pt`, `val_data.pt`, `test_data.pt`
- Data loading tests passing

### Day 3: Training Loop (6-8 hours)

**Tasks**:
1. Create `scripts/train_genomepath.py`:
   - Training loop with progress bars
   - Validation evaluation
   - Checkpoint saving (best model)
   - TensorBoard logging
2. Implement loss functions:
   - Weighted binary cross-entropy
   - Bidirectional loss combination
3. Implement evaluation metrics:
   - Accuracy@K calculation
   - Bidirectional consistency check
   - MRR calculation
4. Run initial training (10 epochs proof-of-concept)

**Deliverables**:
- `train_genomepath.py` (~500 lines)
- `backend/services/genomepath/metrics.py` (~200 lines)
- Initial model checkpoint
- Training logs and TensorBoard visualizations

### Day 4: Hyperparameter Tuning (5-7 hours)

**Tasks**:
1. Grid search over key hyperparameters:
   - Learning rate: [1e-5, 5e-5, 1e-4, 5e-4]
   - Batch size: [16, 32, 64]
   - Dropout: [0.05, 0.1, 0.2]
   - Bridge dimension: [128, 256, 512]
2. Run full 50-epoch training with best config
3. Evaluate on test set
4. Generate performance report

**Deliverables**:
- Hyperparameter search results
- Best model checkpoint (50 epochs)
- Test set evaluation metrics
- `docs/MODEL_TRAINING_REPORT.md`

### Day 5: Inference & Deployment (6-8 hours)

**Tasks**:
1. Create `backend/services/genomepath/inference.py`:
   - Load trained model
   - Predict TK → Genomic targets
   - Predict Genomic → TK practices
   - Batch inference optimization
2. Integrate into existing API (`genomepath_api.py`)
3. Add model caching (load once, reuse)
4. Benchmark inference speed (<200ms target)
5. Create prediction examples for documentation

**Deliverables**:
- `inference.py` (~300 lines)
- Updated `genomepath_api.py` with model predictions
- Inference speed benchmarks
- API documentation with examples
- Production-ready model file

---

## 📦 Dependencies

**New Python Packages** (add to requirements.txt):
```
torch>=2.0.0
torchvision>=0.15.0
sentence-transformers>=2.2.0
tensorboard>=2.13.0
scikit-learn>=1.3.0
tqdm>=4.65.0
```

**Estimated Installation**: 5-10 minutes (PyTorch is 2GB)

---

## 🎯 Success Criteria

### Minimum Viable Model (MVP)

| Metric | Target | Rationale |
|--------|--------|-----------|
| TK→Genomic Accuracy@5 | ≥70% | 7/10 top predictions contain true targets |
| Genomic→TK Accuracy@3 | ≥60% | 2/3 top predictions contain true practices |
| Bidirectional Consistency | ≥75% | Matches current correlator threshold |
| Inference Speed | <200ms | 10x safety margin for Nevada (<2s) |
| Ethical Compliance | 100% | 0 sacred violations, full attribution |

### Production-Ready Model

| Metric | Target | Rationale |
|--------|--------|-----------|
| TK→Genomic Accuracy@5 | ≥80% | 8/10 top predictions correct |
| Genomic→TK Accuracy@3 | ≥70% | 7/10 top predictions correct |
| Bidirectional Consistency | ≥85% | High confidence in both directions |
| Confidence Calibration (ECE) | <0.10 | Reliable confidence scores |
| Dosage-Specific Accuracy | Low ≥75%, Med ≥65%, High ≥55% | Dose-dependent targeting |
| Synergy Recognition | ≥30% | Multi-target optimization |

---

## 💰 Resource Estimates

**Compute Requirements**:
- **Training**: ~2-4 hours on modern GPU (RTX 3060 or better)
- **Inference**: CPU-only (100-200ms per prediction)
- **Storage**: ~500MB for model weights + embeddings

**Development Time**:
- Day 1: Model architecture (6-8 hours)
- Day 2: Data pipeline (5-7 hours)
- Day 3: Training loop (6-8 hours)
- Day 4: Hyperparameter tuning (5-7 hours)
- Day 5: Inference & deployment (6-8 hours)
- **Total**: 28-38 hours over 5 days

**Cost**:
- **Free** (using local GPU or Google Colab free tier)
- Alternative: AWS g4dn.xlarge ~$0.526/hr × 4 hrs = ~$2

---

## 🚀 Next Steps

### Immediate (Today)
1. Install PyTorch and transformers library
2. Create model architecture skeleton
3. Verify SentenceTransformer embeddings work
4. Load 525-correlation dataset

### Tomorrow
1. Implement data preparation pipeline
2. Create PyTorch Dataset classes
3. Test data loading with batching
4. Verify train/val/test splits

### This Week
1. Build complete training loop
2. Run 10-epoch proof-of-concept
3. Evaluate initial results
4. Tune hyperparameters
5. Achieve MVP success criteria

---

## 📝 Deliverables Checklist

- [ ] `backend/services/genomepath/model.py` - Dual-tower transformer
- [ ] `backend/services/genomepath/tokenizer.py` - Text preprocessing
- [ ] `backend/services/genomepath/dataset.py` - PyTorch dataset
- [ ] `backend/services/genomepath/metrics.py` - Evaluation metrics
- [ ] `backend/services/genomepath/inference.py` - Prediction interface
- [ ] `scripts/prepare_training_data.py` - Data preparation
- [ ] `scripts/train_genomepath.py` - Training loop
- [ ] `data/processed/train_data.pt` - Training dataset
- [ ] `data/processed/val_data.pt` - Validation dataset
- [ ] `data/processed/test_data.pt` - Test dataset
- [ ] `models/genomepath_best.pt` - Trained model checkpoint
- [ ] `docs/MODEL_TRAINING_REPORT.md` - Performance analysis
- [ ] Updated `tests/test_genomepath_model.py` - Model tests
- [ ] Updated API endpoints with model inference
- [ ] TensorBoard training visualizations

---

## 🏆 Expected Outcomes

**Week 1 Completion**:
- ✅ Trained GenomePath transformer model (50 epochs)
- ✅ TK→Genomic accuracy ≥70% @top-5
- ✅ Genomic→TK accuracy ≥60% @top-3
- ✅ Inference speed <200ms (10x Nevada requirement)
- ✅ Production-ready model checkpoint
- ✅ API integration complete
- ✅ Ready for Nevada dispensary pilot testing

**Strategic Value**:
- First AI model for TK↔Genomic bidirectional prediction
- Dosage-specific cannabis recommendations
- Multi-target synergy recognition
- 100% ethical compliance maintained
- Foundation for $6,000 MRR Nevada pilot

---

**Status**: 🟢 READY TO BEGIN  
**Prerequisites**: ✅ ALL COMPLETE  
**Next Action**: Install PyTorch and create model architecture skeleton

---

*NeuroBotanica MVP Development - GenomePath Phase 4*  
*Cloak and Quill Research 501(c)(3) | Henderson, Nevada*
