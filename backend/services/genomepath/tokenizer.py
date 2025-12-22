"""
GenomePath Tokenizer - Text Preprocessing for TK and Genomic Features
======================================================================

Converts TK practices and genomic targets into model-ready embeddings.

Components:
- SentenceTransformer for semantic text encoding (384-d)
- TK practice text preprocessing
- Genomic feature concatenation (gene + pathway + tissue + disease = 1536-d)
- Batch processing and caching

License: Proprietary - Trade Secret Protection Required
"""

import torch
import numpy as np
from sentence_transformers import SentenceTransformer
from typing import Dict, List, Tuple, Optional, Union
import json
from pathlib import Path


class GenomePathTokenizer:
    """
    Tokenizer for GenomePath model inputs.
    
    Converts text descriptions into SentenceTransformer embeddings and
    prepares genomic features for model consumption.
    """
    
    def __init__(
        self,
        model_name: str = 'pritamdeka/BioBERT-mnli-snli-scinli-scitail-mednli-stsb',
        cache_dir: Optional[str] = None,
        device: Optional[str] = 'cpu'  # Force CPU for DataLoader compatibility
    ):
        """
        Initialize tokenizer with SentenceTransformer model.
        
        Args:
            model_name: HuggingFace model name
                - 'pritamdeka/BioBERT-mnli-snli-scinli-scitail-mednli-stsb' (default, BioBERT fine-tuned, 768-d)
                - 'microsoft/BiomedNLP-PubMedBERT-base-uncased-abstract-fulltext' (PubMedBERT, 768-d)
                - 'dmis-lab/biobert-v1.1' (BioBERT base, 768-d)
                - 'allenai/scibert_scivocab_uncased' (SciBERT, 768-d)
                - 'sentence-transformers/all-MiniLM-L6-v2' (general, 384-d)
            cache_dir: Directory to cache embeddings
            device: Device for model (default: 'cpu' for DataLoader compatibility)
        """
        # Always use CPU for embedding generation in DataLoader
        # Model tensors will be moved to GPU separately
        device = 'cpu'
        
        print(f"Loading SentenceTransformer model: {model_name}")
        print(f"Device: {device}")
        
        self.model = SentenceTransformer(model_name, device=device)
        self.embedding_dim = self.model.get_sentence_embedding_dimension()
        
        # BioBERT models output 768-d embeddings
        print(f"✅ Tokenizer initialized with {self.embedding_dim}-d embeddings")
        
        self.cache_dir = Path(cache_dir) if cache_dir else None
        if self.cache_dir:
            self.cache_dir.mkdir(parents=True, exist_ok=True)
    
    def encode_tk_practice(
        self,
        practice_name: str,
        indications: List[str],
        preparation_method: str = "",
        cultural_context: str = "",
        source_community: str = ""
    ) -> Dict[str, torch.Tensor]:
        """
        Encode TK practice into model-ready format.
        
        Args:
            practice_name: Name of traditional practice
            indications: List of therapeutic indications
            preparation_method: How it's prepared
            cultural_context: Cultural significance
            source_community: Community attribution
            
        Returns:
            Dictionary with:
            - embeddings: [seq_len, 384] tensor
            - attention_mask: [seq_len] tensor (1=real, 0=padding)
            - metadata: Original text components
        """
        # Build text sequence for encoding
        text_components = []
        
        # Practice name (always include)
        text_components.append(practice_name)
        
        # Indications (core therapeutic information)
        if indications:
            indications_text = "Used for: " + ", ".join(indications)
            text_components.append(indications_text)
        
        # Preparation method
        if preparation_method:
            text_components.append(f"Preparation: {preparation_method}")
        
        # Cultural context (for cultural sensitivity preservation)
        if cultural_context:
            text_components.append(f"Cultural context: {cultural_context}")
        
        # Source community (for attribution)
        if source_community:
            text_components.append(f"Source: {source_community}")
        
        # Encode each component separately (preserves structure)
        embeddings = self.model.encode(
            text_components,
            convert_to_tensor=True,
            show_progress_bar=False
        )  # [seq_len, 384]
        
        # Create attention mask (all real tokens)
        attention_mask = torch.ones(len(text_components), dtype=torch.long)
        
        return {
            'embeddings': embeddings,
            'attention_mask': attention_mask,
            'metadata': {
                'practice_name': practice_name,
                'indications': indications,
                'preparation_method': preparation_method,
                'cultural_context': cultural_context,
                'source_community': source_community,
                'num_components': len(text_components)
            }
        }
    
    def encode_genomic_target(
        self,
        gene_id: str,
        pathways: List[str],
        tissues: List[str],
        diseases: List[str]
    ) -> Dict[str, torch.Tensor]:
        """
        Encode genomic target into model-ready format.
        
        Creates 1536-d embedding from 4 modalities:
        - Gene ID: 384-d
        - Pathways: 384-d (concatenated)
        - Tissues: 384-d (concatenated)
        - Diseases: 384-d (concatenated)
        
        Args:
            gene_id: Gene identifier (e.g., "CNR1", "CB2")
            pathways: List of pathways involved
            tissues: List of tissue expression
            diseases: List of disease associations
            
        Returns:
            Dictionary with:
            - embeddings: [seq_len, 1536] tensor (4 × 384-d)
            - attention_mask: [seq_len] tensor
            - metadata: Original components
        """
        # Encode each modality
        modality_texts = []
        
        # 1. Gene ID
        modality_texts.append(f"Gene: {gene_id}")
        
        # 2. Pathways
        if pathways:
            pathway_text = "Pathways: " + ", ".join(pathways)
        else:
            pathway_text = "Pathways: unknown"
        modality_texts.append(pathway_text)
        
        # 3. Tissues
        if tissues:
            tissue_text = "Expressed in: " + ", ".join(tissues)
        else:
            tissue_text = "Expressed in: unknown"
        modality_texts.append(tissue_text)
        
        # 4. Diseases
        if diseases:
            disease_text = "Associated with: " + ", ".join(diseases)
        else:
            disease_text = "Associated with: unknown"
        modality_texts.append(disease_text)
        
        # Encode all modalities
        modality_embeddings = self.model.encode(
            modality_texts,
            convert_to_tensor=True,
            show_progress_bar=False
        )  # [4, 384]
        
        # Concatenate to create 1536-d representation
        concatenated = modality_embeddings.flatten().unsqueeze(0)  # [1, 1536]
        
        # Attention mask (single token after concatenation)
        attention_mask = torch.ones(1, dtype=torch.long)
        
        return {
            'embeddings': concatenated,
            'attention_mask': attention_mask,
            'metadata': {
                'gene_id': gene_id,
                'pathways': pathways,
                'tissues': tissues,
                'diseases': diseases,
                'num_modalities': 4
            }
        }
    
    def batch_encode_tk_practices(
        self,
        practices: List[Dict[str, any]],
        max_length: int = 128
    ) -> Dict[str, torch.Tensor]:
        """
        Batch encode multiple TK practices with padding.
        
        Args:
            practices: List of practice dictionaries
            max_length: Maximum sequence length (for padding)
            
        Returns:
            Batched tensors with padding
        """
        all_embeddings = []
        all_masks = []
        
        for practice in practices:
            encoded = self.encode_tk_practice(
                practice_name=practice.get('practice_name', ''),
                indications=practice.get('indications', []),
                preparation_method=practice.get('preparation_method', ''),
                cultural_context=practice.get('cultural_context', ''),
                source_community=practice.get('source_community', '')
            )
            all_embeddings.append(encoded['embeddings'])
            all_masks.append(encoded['attention_mask'])
        
        # Pad to max_length
        padded_embeddings = self._pad_embeddings(all_embeddings, max_length)
        padded_masks = self._pad_masks(all_masks, max_length)
        
        return {
            'embeddings': padded_embeddings,
            'attention_mask': padded_masks
        }
    
    def batch_encode_genomic_targets(
        self,
        targets: List[Dict[str, any]]
    ) -> Dict[str, torch.Tensor]:
        """
        Batch encode multiple genomic targets.
        
        Args:
            targets: List of genomic target dictionaries
            
        Returns:
            Batched tensors (no padding needed, all 1536-d)
        """
        all_embeddings = []
        all_masks = []
        
        for target in targets:
            encoded = self.encode_genomic_target(
                gene_id=target.get('gene_id', ''),
                pathways=target.get('pathways', []),
                tissues=target.get('tissues', []),
                diseases=target.get('diseases', [])
            )
            all_embeddings.append(encoded['embeddings'])
            all_masks.append(encoded['attention_mask'])
        
        # Stack (all same length for genomic)
        batched_embeddings = torch.cat(all_embeddings, dim=0).unsqueeze(1)  # [batch, 1, 1536]
        batched_masks = torch.stack(all_masks, dim=0)  # [batch, 1]
        
        return {
            'embeddings': batched_embeddings,
            'attention_mask': batched_masks
        }
    
    def _pad_embeddings(
        self,
        embeddings_list: List[torch.Tensor],
        max_length: int
    ) -> torch.Tensor:
        """Pad embeddings to max_length."""
        batch_size = len(embeddings_list)
        embedding_dim = embeddings_list[0].size(-1)
        
        # Create padded tensor
        padded = torch.zeros(batch_size, max_length, embedding_dim)
        
        for i, emb in enumerate(embeddings_list):
            seq_len = min(emb.size(0), max_length)
            padded[i, :seq_len, :] = emb[:seq_len, :]
        
        return padded
    
    def _pad_masks(
        self,
        masks_list: List[torch.Tensor],
        max_length: int
    ) -> torch.Tensor:
        """Pad attention masks to max_length."""
        batch_size = len(masks_list)
        
        # Create padded masks (0 = padding)
        padded = torch.zeros(batch_size, max_length, dtype=torch.long)
        
        for i, mask in enumerate(masks_list):
            seq_len = min(mask.size(0), max_length)
            padded[i, :seq_len] = mask[:seq_len]
        
        return padded
    
    def get_embedding_dim(self) -> int:
        """Get embedding dimension (384)."""
        return self.embedding_dim
    
    def get_genomic_embedding_dim(self) -> int:
        """Get genomic embedding dimension (1536 = 4 × 384)."""
        return self.embedding_dim * 4


def load_tk_practices(json_path: str) -> List[Dict]:
    """Load TK practices from JSON file."""
    with open(json_path, 'r', encoding='utf-8') as f:
        data = json.load(f)
    return data if isinstance(data, list) else [data]


def load_genomic_targets(json_path: str) -> List[Dict]:
    """Load genomic targets from JSON file."""
    with open(json_path, 'r', encoding='utf-8') as f:
        data = json.load(f)
    return data if isinstance(data, list) else [data]


if __name__ == "__main__":
    print("=" * 70)
    print("GenomePath Tokenizer - Verification Test")
    print("=" * 70)
    
    # Initialize tokenizer
    tokenizer = GenomePathTokenizer()
    
    # Test TK encoding
    print("\n📚 Testing TK Practice Encoding...")
    tk_encoded = tokenizer.encode_tk_practice(
        practice_name="Cannabis for chronic pain relief",
        indications=["chronic pain", "inflammation", "neuropathic pain"],
        preparation_method="Infusion or tincture",
        cultural_context="Traditional use in pain management",
        source_community="traditional_western_herbal_medicine"
    )
    
    print(f"  TK embeddings shape: {tk_encoded['embeddings'].shape}")
    print(f"  TK attention mask: {tk_encoded['attention_mask'].shape}")
    print(f"  Number of components: {tk_encoded['metadata']['num_components']}")
    assert tk_encoded['embeddings'].shape[1] == 384, "TK embeddings should be 384-d"
    print("  ✅ TK encoding verified (384-d)")
    
    # Test Genomic encoding
    print("\n🧬 Testing Genomic Target Encoding...")
    genomic_encoded = tokenizer.encode_genomic_target(
        gene_id="CNR1",
        pathways=["endocannabinoid system", "GPCR signaling"],
        tissues=["brain", "CNS", "peripheral nervous system"],
        diseases=["chronic pain", "anxiety", "epilepsy"]
    )
    
    print(f"  Genomic embeddings shape: {genomic_encoded['embeddings'].shape}")
    print(f"  Genomic attention mask: {genomic_encoded['attention_mask'].shape}")
    assert genomic_encoded['embeddings'].shape[1] == 1536, "Genomic embeddings should be 1536-d"
    print("  ✅ Genomic encoding verified (1536-d = 4 × 384-d)")
    
    # Test batch encoding
    print("\n📦 Testing Batch Encoding...")
    tk_practices = [
        {
            'practice_name': 'Cannabis for pain',
            'indications': ['pain'],
            'preparation_method': 'Tea',
            'cultural_context': 'Traditional',
            'source_community': 'community1'
        },
        {
            'practice_name': 'Hemp for anxiety',
            'indications': ['anxiety', 'stress'],
            'preparation_method': 'Oil',
            'cultural_context': 'Modern',
            'source_community': 'community2'
        }
    ]
    
    batched_tk = tokenizer.batch_encode_tk_practices(tk_practices, max_length=128)
    print(f"  Batched TK shape: {batched_tk['embeddings'].shape}")
    print(f"  Batched mask shape: {batched_tk['attention_mask'].shape}")
    assert batched_tk['embeddings'].shape == (2, 128, 384), "Batch shape mismatch"
    print("  ✅ Batch encoding verified")
    
    print("\n" + "=" * 70)
    print("✅ All tokenizer tests passed!")
    print(f"✅ TK embedding dimension: {tokenizer.get_embedding_dim()} (384-d)")
    print(f"✅ Genomic embedding dimension: {tokenizer.get_genomic_embedding_dim()} (1536-d)")
    print("✅ Ready for dataset creation")
    print("=" * 70)
