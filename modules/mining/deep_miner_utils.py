import os
import torch
import re
from transformers import AutoTokenizer, AutoModel
from Bio.Seq import Seq


def _resolve_device():
    """Resolve compute device: CUDA (NVIDIA), ROCm (AMD), or CPU."""
    forced = os.getenv("FAMILY_DEVICE", "").strip().lower()
    if forced in ("cuda", "cpu"):
        return forced
    # PyTorch CUDA: works for NVIDIA; ROCm build on Linux also reports cuda.is_available()
    if torch.cuda.is_available():
        return "cuda"
    return "cpu"


class DeepEngine:
    def __init__(self, model_name="facebook/esm2_t12_35M_UR50D"):
        """
        Initializes the ESM-2 Model.
        Supports CUDA (NVIDIA), ROCm (AMD via PyTorch ROCm build), or CPU.
        Set FAMILY_DEVICE=cuda or cpu to force device.
        """
        self.device = _resolve_device()
        print(f"[*] Loading Deep Learning Model: {model_name} on {self.device.upper()}...")
        
        self.tokenizer = AutoTokenizer.from_pretrained(model_name)
        self.model = AutoModel.from_pretrained(model_name).to(self.device)
        self.model.eval() 
        
        # --- THE GOLDEN REFERENCE (Cas13d) ---
        self.ref_seq_segment = "RHYLDEIIEQISEFSKRVILADANLDKVLSAYNKHRDKPIREQAENIIHLFTLTNLGAPAAFKYFDTTIDRKRYTSTKEVLDATLIHQSITGLYETRIDLSQLGGD"
        
        # Calculate reference vector once on GPU
        self.ref_vector = self._get_embedding(self.ref_seq_segment)
        print(f"[+] Reference Vector Calculated on {self.device}.")

    def _get_embedding(self, sequence):
        """Converts Amino Acid sequence into Vector on GPU."""
        # Move inputs to the same device as model
        inputs = self.tokenizer(sequence, return_tensors="pt").to(self.device)
        
        with torch.no_grad():
            outputs = self.model(**inputs)
            
        # Mean pooling
        return outputs.last_hidden_state.mean(dim=1)

    def get_embeddings_batch(self, sequences, max_len=1000, batch_size=50):
        """
        Compute ESM-2 embeddings for a list of sequences.
        Processes in batches to avoid OOM. Returns tensor of shape (N, embed_dim).
        """
        if not sequences:
            return None, []
        processed = []
        valid_indices = []
        for i, seq in enumerate(sequences):
            s = str(seq).strip()
            if len(s) < 50:
                continue
            processed.append(s[:max_len])
            valid_indices.append(i)
        if not processed:
            return None, []
        batch_size = int(os.getenv("EMBED_BATCH_SIZE", batch_size))
        all_embeddings = []
        try:
            for start in range(0, len(processed), batch_size):
                batch = processed[start : start + batch_size]
                inputs = self.tokenizer(
                    batch,
                    return_tensors="pt",
                    padding=True,
                    truncation=True,
                    max_length=max_len,
                ).to(self.device)
                with torch.no_grad():
                    outputs = self.model(**inputs)
                mask = inputs["attention_mask"]
                hidden = outputs.last_hidden_state
                mask_expanded = mask.unsqueeze(-1).expand(hidden.size()).float()
                summed = torch.sum(hidden * mask_expanded, dim=1)
                sum_mask = mask_expanded.sum(dim=1).clamp(min=1e-9)
                emb = summed / sum_mask
                all_embeddings.append(emb)
            embeddings = torch.cat(all_embeddings, dim=0)
            return embeddings, valid_indices
        except Exception as e:
            print(f"[!] Embedding failed: {e}")
            return None, []

    def score_candidate(self, candidate_seq):
        """Returns Similarity Score (0.0 to 1.0)."""
        if len(candidate_seq) < 300: return 0.0
        process_seq = candidate_seq[:1000] # Truncate for speed
        
        try:
            cand_vector = self._get_embedding(process_seq)
            cosine_sim = torch.nn.functional.cosine_similarity(self.ref_vector, cand_vector)
            return cosine_sim.item()
        except Exception:
            return 0.0

class NeighborhoodWatch:
    """Finds CRISPR Arrays. Uses multiple chunk sizes (24-32bp) and accepts 2+ repeats."""

    def has_crispr_array(self, dna_sequence):
        seq_str = str(dna_sequence)
        length = len(seq_str)
        if length < 400:
            return False

        min_repeats = 3
        if length < 1500:
            min_repeats = 2

        for chunk_size in (24, 28, 32):
            seen = {}
            for i in range(0, length - chunk_size):
                chunk = seq_str[i : i + chunk_size]
                if chunk in seen:
                    seen[chunk] += 1
                    if seen[chunk] >= min_repeats:
                        return True
                else:
                    seen[chunk] = 1
        return False