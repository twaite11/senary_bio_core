import os
import torch
import re
from pathlib import Path
from transformers import AutoTokenizer, AutoModel
from Bio import SeqIO
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
    def __init__(self, model_name="facebook/esm2_t12_35M_UR50D", reference_fasta=None):
        """
        Initializes the ESM-2 Model.
        reference_fasta: optional path to FASTA with named refs (e.g. RfxCas13d, PspCas13a).
        If set, scores against each ref and returns max similarity (closest match).
        Set FAMILY_DEVICE=cuda or cpu to force device.
        """
        self.device = _resolve_device()
        print(f"[*] Loading Deep Learning Model: {model_name} on {self.device.upper()}...")
        
        self.tokenizer = AutoTokenizer.from_pretrained(model_name)
        self.model = AutoModel.from_pretrained(model_name).to(self.device)
        self.model.eval()

        self.ref_names = []
        self.ref_vectors = []
        path = reference_fasta or os.getenv("ESM_REFERENCE_FASTA", "").strip()
        if path and Path(path).exists():
            for rec in SeqIO.parse(path, "fasta"):
                seq = str(rec.seq).strip()
                if len(seq) >= 50:
                    self.ref_names.append(rec.id)
                    self.ref_vectors.append(self._get_embedding(seq[:1000]))
            if self.ref_vectors:
                print(f"[+] Loaded {len(self.ref_vectors)} reference(s) from {path}: {', '.join(self.ref_names)}")
        if not self.ref_vectors:
            # Default: single Cas13d-like segment (RfxCas13d-like)
            self.ref_seq_segment = "RHYLDEIIEQISEFSKRVILADANLDKVLSAYNKHRDKPIREQAENIIHLFTLTNLGAPAAFKYFDTTIDRKRYTSTKEVLDATLIHQSITGLYETRIDLSQLGGD"
            self.ref_names = ["Cas13d_ref"]
            self.ref_vectors = [self._get_embedding(self.ref_seq_segment)]
            print(f"[+] Using default single reference (Cas13d-like) on {self.device}.")

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
        """Returns max similarity (0.0 to 1.0) over all reference vectors (closest match)."""
        score, _ = self.score_candidate_with_ref(candidate_seq)
        return score

    def score_candidates_batch(self, sequences, batch_size=32):
        """
        Score multiple sequences in batches (faster on GPU/CPU than one-by-one).
        Returns list of float scores, same length as sequences; 0.0 for short/invalid.
        """
        if not sequences:
            return []
        # Filter short; we'll return 0.0 for them
        valid_seqs = []
        valid_indices = []
        for i, s in enumerate(sequences):
            s = str(s).strip()
            if len(s) < 300:
                continue
            valid_seqs.append(s[:1000])
            valid_indices.append(i)
        if not valid_seqs:
            return [0.0] * len(sequences)
        embeddings, _ = self.get_embeddings_batch(valid_seqs, max_len=1000, batch_size=batch_size)
        if embeddings is None:
            return [0.0] * len(sequences)
        scores_by_valid_idx = []
        for j in range(embeddings.size(0)):
            cand_vec = embeddings[j : j + 1]
            best = 0.0
            for ref_vec in self.ref_vectors:
                sim = torch.nn.functional.cosine_similarity(ref_vec, cand_vec).item()
                if sim > best:
                    best = sim
            scores_by_valid_idx.append(best)
        # Map back to original indices
        result = [0.0] * len(sequences)
        for k, idx in enumerate(valid_indices):
            result[idx] = scores_by_valid_idx[k]
        return result

    def score_candidate_with_ref(self, candidate_seq):
        """Returns (max_similarity, ref_name) for closest reference (RfxCas13d, PspCas13a, etc.)."""
        if len(candidate_seq) < 300:
            return 0.0, self.ref_names[0] if self.ref_names else ""
        process_seq = candidate_seq[:1000]
        try:
            cand_vector = self._get_embedding(process_seq)
            best = 0.0
            best_name = self.ref_names[0] if self.ref_names else ""
            for name, ref_vec in zip(self.ref_names, self.ref_vectors):
                sim = torch.nn.functional.cosine_similarity(ref_vec, cand_vector).item()
                if sim > best:
                    best = sim
                    best_name = name
            return best, best_name
        except Exception:
            return 0.0, self.ref_names[0] if self.ref_names else ""

    def passes_diversity_band(self, score):
        """
        When diversity mode is on (ESM_SIMILARITY_CEILING set), accept only scores
        within [floor, ceiling] to favor distant homologs over "perfect" matches.
        """
        floor_s = os.getenv("ESM_SIMILARITY_FLOOR", "")
        ceiling_s = os.getenv("ESM_SIMILARITY_CEILING", "")
        if not ceiling_s:
            return True  # no band: any score above threshold is ok
        try:
            floor = float(floor_s) if floor_s else 0.0
            ceiling = float(ceiling_s)
            return floor <= score <= ceiling
        except ValueError:
            return True

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

    def get_repeat_domains(self, dna_sequence):
        """
        Extract CRISPR repeat sequences from DNA (same logic as has_crispr_array).
        Returns list of repeat sequences (chunks that appear >= min_repeats) for synthesis metadata.
        """
        seq_str = str(dna_sequence)
        length = len(seq_str)
        if length < 400:
            return []

        min_repeats = 3
        if length < 1500:
            min_repeats = 2

        repeats = []
        for chunk_size in (24, 28, 32):
            seen = {}
            for i in range(0, length - chunk_size):
                chunk = seq_str[i : i + chunk_size]
                if chunk in seen:
                    seen[chunk] += 1
                    if seen[chunk] >= min_repeats and chunk not in repeats:
                        repeats.append(chunk)
                else:
                    seen[chunk] = 1
        return repeats

    def get_repeat_regions(self, dna_sequence):
        """
        Return nucleotide (start, end) for each CRISPR repeat array on the contig.
        Used to require CRISPR within N kb upstream/downstream of the ORF (e.g. 10 kb).
        """
        seq_str = str(dna_sequence)
        length = len(seq_str)
        if length < 400:
            return []

        min_repeats = 3
        if length < 1500:
            min_repeats = 2

        regions = []
        for chunk_size in (24, 28, 32):
            # seen[chunk] = list of start positions (so we can get span)
            seen = {}
            for i in range(0, length - chunk_size):
                chunk = seq_str[i : i + chunk_size]
                if chunk not in seen:
                    seen[chunk] = []
                seen[chunk].append(i)
            for chunk, positions in seen.items():
                if len(positions) >= min_repeats:
                    start_nt = min(positions)
                    end_nt = max(positions) + chunk_size
                    regions.append((start_nt, end_nt))
        return regions