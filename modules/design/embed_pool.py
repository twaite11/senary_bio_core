"""
Embed a FASTA pool into ESM-2 latent space. Saves embeddings and metadata for
interpolation and downstream design.
"""
import os
import sys
import argparse
import numpy as np
from pathlib import Path
from Bio import SeqIO

# Allow running from project root
sys.path.insert(0, os.getcwd())
from modules.mining.deep_miner_utils import DeepEngine


def embed_fasta(fasta_path: str, output_dir: str, max_len: int = 1000, batch_size: int = 50):
    """
    Load FASTA, embed all sequences with ESM-2, save embeddings and metadata.
    """
    fasta_path = Path(fasta_path)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    records = list(SeqIO.parse(fasta_path, "fasta"))
    if not records:
        print(f"[!] No sequences in {fasta_path}")
        return None, None

    ids = [r.id for r in records]
    sequences = [str(r.seq) for r in records]

    engine = DeepEngine()
    embeddings, valid_indices = engine.get_embeddings_batch(sequences, max_len=max_len, batch_size=batch_size)
    if embeddings is None:
        print("[!] Embedding failed.")
        return None, None

    # embeddings[j] corresponds to valid_indices[j]
    all_ids = [ids[i] for i in valid_indices]
    emb_matrix = embeddings.cpu().numpy()

    np.save(output_dir / "embeddings.npy", emb_matrix)
    with open(output_dir / "sequence_ids.txt", "w") as f:
        for sid in all_ids:
            f.write(sid + "\n")
    print(f"[+] Saved embeddings for {len(all_ids)} sequences to {output_dir}")
    return emb_matrix, all_ids


def main():
    parser = argparse.ArgumentParser(description="Embed FASTA pool into ESM-2 latent space.")
    parser.add_argument("--input", default="data/raw_sequences/deep_hits_latest.fasta",
                        help="Input FASTA (mined or merged pool)")
    parser.add_argument("--output-dir", default="data/design/embeddings",
                        help="Output directory for embeddings.npy and sequence_ids.txt")
    parser.add_argument("--max-len", type=int, default=1000, help="Max sequence length for ESM-2")
    parser.add_argument("--batch-size", type=int, default=50, help="Batch size for embedding")
    args = parser.parse_args()

    if not Path(args.input).exists():
        print(f"[!] Input not found: {args.input}. Run mining first or pass a FASTA.")
        return 1
    embed_fasta(args.input, args.output_dir, max_len=args.max_len, batch_size=args.batch_size)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
