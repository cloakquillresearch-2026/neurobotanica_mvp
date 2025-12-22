"""
Hybrid Chemical-Text GenomePath Model
=====================================

3-tower architecture:
1. Text encoder (BioBERT): TK practice descriptions
2. Chemical encoder (Morgan fingerprints): Active compounds
3. Genomic encoder (BioBERT): Gene/pathway features

Fusion strategy: Early fusion of text + chemical, then cross-attention with genomic
"""

import torch
import torch.nn as nn
import torch.nn.functional as F
from typing import Dict, Optional
from dataclasses import dataclass


@dataclass
class HybridModelConfig:
    """Hybrid chemical-text model configuration."""
    
    # Input dimensions
    text_dim: int = 768  # BioBERT
    chemical_fingerprint_dim: int = 2048  # Morgan ECFP4
    chemical_descriptors_dim: int = 7  # Molecular properties
    genomic_dim: int = 3072  # 4 × 768-d
    
    # Encoded dimensions
    text_hidden_dim: int = 256
    chemical_hidden_dim: int = 256
    fusion_dim: int = 512  # text + chemical combined
    genomic_hidden_dim: int = 512
    
    # Transformer
    num_layers: int = 4
    num_heads: int = 8
    ffn_dim: int = 1024
    dropout: float = 0.1
    
    # Prediction
    num_genomic_targets: int = 46
    num_tk_practices: int = 30


class ChemicalFingerprintEncoder(nn.Module):
    """Encodes Morgan fingerprints + molecular descriptors"""
    
    def __init__(self, config: HybridModelConfig):
        super().__init__()
        
        # Fingerprint projection
        self.fingerprint_proj = nn.Sequential(
            nn.Linear(config.chemical_fingerprint_dim, 512),
            nn.GELU(),
            nn.Dropout(config.dropout),
            nn.Linear(512, config.chemical_hidden_dim)
        )
        
        # Descriptor projection
        self.descriptor_proj = nn.Sequential(
            nn.Linear(config.chemical_descriptors_dim, 64),
            nn.GELU(),
            nn.Linear(64, config.chemical_hidden_dim)
        )
        
        # Fusion
        self.fusion = nn.Sequential(
            nn.Linear(config.chemical_hidden_dim * 2, config.chemical_hidden_dim),
            nn.LayerNorm(config.chemical_hidden_dim),
            nn.GELU(),
            nn.Dropout(config.dropout)
        )
    
    def forward(
        self,
        fingerprint: torch.Tensor,  # [batch, 2048]
        descriptors: torch.Tensor,  # [batch, 7]
        is_valid: Optional[torch.Tensor] = None  # [batch]
    ) -> torch.Tensor:
        """
        Args:
            fingerprint: Morgan fingerprint [batch, 2048]
            descriptors: Molecular descriptors [batch, 7]
            is_valid: Validity mask [batch] (1.0 if valid, 0.0 if placeholder)
            
        Returns:
            Chemical encoding [batch, chemical_hidden_dim]
        """
        fp_encoded = self.fingerprint_proj(fingerprint)  # [batch, 256]
        desc_encoded = self.descriptor_proj(descriptors)  # [batch, 256]
        
        # Concatenate and fuse
        combined = torch.cat([fp_encoded, desc_encoded], dim=1)  # [batch, 512]
        fused = self.fusion(combined)  # [batch, 256]
        
        # Zero out invalid compounds
        if is_valid is not None:
            fused = fused * is_valid.unsqueeze(-1)
        
        return fused


class TextEncoder(nn.Module):
    """Encodes TK practice text (BioBERT embeddings)"""
    
    def __init__(self, config: HybridModelConfig):
        super().__init__()
        
        # Project text to hidden dim
        self.proj = nn.Sequential(
            nn.Linear(config.text_dim, config.text_hidden_dim),
            nn.LayerNorm(config.text_hidden_dim),
            nn.GELU(),
            nn.Dropout(config.dropout)
        )
        
        # Simple transformer for multi-component text
        encoder_layer = nn.TransformerEncoderLayer(
            d_model=config.text_hidden_dim,
            nhead=config.num_heads // 2,  # Smaller for text-only
            dim_feedforward=config.ffn_dim // 2,
            dropout=config.dropout,
            activation='gelu',
            batch_first=True,
            norm_first=True
        )
        self.transformer = nn.TransformerEncoder(encoder_layer, num_layers=2)
    
    def forward(
        self,
        text_embeddings: torch.Tensor,  # [batch, seq_len, 768]
        attention_mask: Optional[torch.Tensor] = None  # [batch, seq_len]
    ) -> torch.Tensor:
        """
        Returns:
            [batch, text_hidden_dim] pooled representation
        """
        x = self.proj(text_embeddings)  # [batch, seq_len, 256]
        
        # Create padding mask for transformer
        if attention_mask is not None:
            src_key_padding_mask = (attention_mask == 0)
        else:
            src_key_padding_mask = None
        
        x = self.transformer(x, src_key_padding_mask=src_key_padding_mask)
        
        # Mean pool
        if attention_mask is not None:
            mask_expanded = attention_mask.unsqueeze(-1).float()
            x = (x * mask_expanded).sum(1) / mask_expanded.sum(1).clamp(min=1)
        else:
            x = x.mean(1)
        
        return x  # [batch, 256]


class TKFusionEncoder(nn.Module):
    """Fuses text + chemical representations"""
    
    def __init__(self, config: HybridModelConfig):
        super().__init__()
        
        self.fusion = nn.Sequential(
            nn.Linear(config.text_hidden_dim + config.chemical_hidden_dim, config.fusion_dim),
            nn.LayerNorm(config.fusion_dim),
            nn.GELU(),
            nn.Dropout(config.dropout),
            nn.Linear(config.fusion_dim, config.fusion_dim),
            nn.LayerNorm(config.fusion_dim)
        )
    
    def forward(
        self,
        text_encoding: torch.Tensor,  # [batch, 256]
        chemical_encoding: torch.Tensor  # [batch, 256]
    ) -> torch.Tensor:
        combined = torch.cat([text_encoding, chemical_encoding], dim=1)  # [batch, 512]
        return self.fusion(combined)  # [batch, 512]


class GenomicEncoder(nn.Module):
    """Encodes genomic features (gene + pathways + tissues + diseases)"""
    
    def __init__(self, config: HybridModelConfig):
        super().__init__()
        
        self.proj = nn.Sequential(
            nn.Linear(config.genomic_dim, config.genomic_hidden_dim),
            nn.LayerNorm(config.genomic_hidden_dim),
            nn.GELU(),
            nn.Dropout(config.dropout)
        )
    
    def forward(self, genomic_embeddings: torch.Tensor) -> torch.Tensor:
        return self.proj(genomic_embeddings)  # [batch, 512]


class HybridGenomePathModel(nn.Module):
    """
    3-tower chemical-text-genomic fusion model
    
    Architecture:
    1. Text tower: BioBERT → Transformer → Pool
    2. Chemical tower: Fingerprint + Descriptors → MLP
    3. Fusion: Text + Chemical → Combined TK representation
    4. Genomic tower: Multi-modal genomic features → MLP
    5. Cross-attention: TK ↔ Genomic
    6. Prediction heads: Bidirectional TK→G and G→TK
    """
    
    def __init__(self, config: Optional[HybridModelConfig] = None):
        super().__init__()
        
        self.config = config or HybridModelConfig()
        
        # Encoders
        self.text_encoder = TextEncoder(self.config)
        self.chemical_encoder = ChemicalFingerprintEncoder(self.config)
        self.tk_fusion = TKFusionEncoder(self.config)
        self.genomic_encoder = GenomicEncoder(self.config)
        
        # Cross-attention (TK ↔ Genomic)
        self.tk_to_genomic_attn = nn.MultiheadAttention(
            embed_dim=self.config.fusion_dim,
            num_heads=self.config.num_heads,
            dropout=self.config.dropout,
            batch_first=True
        )
        self.genomic_to_tk_attn = nn.MultiheadAttention(
            embed_dim=self.config.genomic_hidden_dim,
            num_heads=self.config.num_heads,
            dropout=self.config.dropout,
            batch_first=True
        )
        
        # Prediction heads
        self.tk_to_genomic_head = nn.Sequential(
            nn.Linear(self.config.fusion_dim, self.config.fusion_dim),
            nn.LayerNorm(self.config.fusion_dim),
            nn.GELU(),
            nn.Dropout(self.config.dropout),
            nn.Linear(self.config.fusion_dim, self.config.num_genomic_targets)
        )
        
        self.genomic_to_tk_head = nn.Sequential(
            nn.Linear(self.config.genomic_hidden_dim, self.config.genomic_hidden_dim),
            nn.LayerNorm(self.config.genomic_hidden_dim),
            nn.GELU(),
            nn.Dropout(self.config.dropout),
            nn.Linear(self.config.genomic_hidden_dim, self.config.num_tk_practices)
        )
    
    def forward(
        self,
        # TK inputs
        tk_text_embeddings: torch.Tensor,  # [batch, seq_len, 768]
        tk_attention_mask: torch.Tensor,  # [batch, seq_len]
        tk_chemical_fingerprint: torch.Tensor,  # [batch, 2048]
        tk_chemical_descriptors: torch.Tensor,  # [batch, 7]
        tk_chemical_valid: torch.Tensor,  # [batch]
        # Genomic inputs
        genomic_embeddings: torch.Tensor,  # [batch, 3072]
        genomic_attention_mask: Optional[torch.Tensor] = None
    ) -> Dict[str, torch.Tensor]:
        """
        Forward pass through hybrid model.
        
        Returns:
            Dictionary with:
            - tk_to_genomic_logits: [batch, 46] predictions
            - genomic_to_tk_logits: [batch, 30] predictions
            - tk_embedding: [batch, 512] TK representation
            - genomic_embedding: [batch, 512] Genomic representation
        """
        # Encode TK components
        text_encoded = self.text_encoder(tk_text_embeddings, tk_attention_mask)  # [batch, 256]
        chemical_encoded = self.chemical_encoder(
            tk_chemical_fingerprint,
            tk_chemical_descriptors,
            tk_chemical_valid
        )  # [batch, 256]
        
        # Fuse TK representations
        tk_fused = self.tk_fusion(text_encoded, chemical_encoded)  # [batch, 512]
        
        # Encode genomic features
        genomic_encoded = self.genomic_encoder(genomic_embeddings)  # [batch, 512]
        
        # Cross-attention (expand to sequence dimension for attention)
        tk_seq = tk_fused.unsqueeze(1)  # [batch, 1, 512]
        genomic_seq = genomic_encoded.unsqueeze(1)  # [batch, 1, 512]
        
        # TK attends to genomic
        tk_attended, _ = self.tk_to_genomic_attn(
            tk_seq, genomic_seq, genomic_seq
        )  # [batch, 1, 512]
        tk_attended = tk_attended.squeeze(1)  # [batch, 512]
        
        # Genomic attends to TK
        genomic_attended, _ = self.genomic_to_tk_attn(
            genomic_seq, tk_seq, tk_seq
        )  # [batch, 1, 512]
        genomic_attended = genomic_attended.squeeze(1)  # [batch, 512]
        
        # Predictions
        tk_to_genomic_logits = self.tk_to_genomic_head(tk_attended)  # [batch, 46]
        genomic_to_tk_logits = self.genomic_to_tk_head(genomic_attended)  # [batch, 30]
        
        return {
            'tk_to_genomic_logits': tk_to_genomic_logits,
            'genomic_to_tk_logits': genomic_to_tk_logits,
            'tk_embedding': tk_fused,
            'genomic_embedding': genomic_encoded
        }
    
    def count_parameters(self) -> int:
        """Count total trainable parameters"""
        return sum(p.numel() for p in self.parameters() if p.requires_grad)


if __name__ == "__main__":
    # Test model
    config = HybridModelConfig()
    model = HybridGenomePathModel(config)
    
    print(f"\nHybrid GenomePath Model")
    print(f"=" * 70)
    print(f"Total parameters: {model.count_parameters():,}")
    print(f"\nConfiguration:")
    print(f"  Text dim: {config.text_dim}")
    print(f"  Chemical fingerprint: {config.chemical_fingerprint_dim}")
    print(f"  Chemical descriptors: {config.chemical_descriptors_dim}")
    print(f"  Genomic dim: {config.genomic_dim}")
    print(f"  Fusion dim: {config.fusion_dim}")
    
    # Test forward pass
    batch_size = 4
    seq_len = 5
    
    dummy_inputs = {
        'tk_text_embeddings': torch.randn(batch_size, seq_len, config.text_dim),
        'tk_attention_mask': torch.ones(batch_size, seq_len),
        'tk_chemical_fingerprint': torch.randn(batch_size, config.chemical_fingerprint_dim),
        'tk_chemical_descriptors': torch.randn(batch_size, config.chemical_descriptors_dim),
        'tk_chemical_valid': torch.ones(batch_size),
        'genomic_embeddings': torch.randn(batch_size, config.genomic_dim)
    }
    
    outputs = model(**dummy_inputs)
    
    print(f"\nOutput shapes:")
    for k, v in outputs.items():
        print(f"  {k}: {v.shape}")
    
    print(f"\nModel ready for training!")
