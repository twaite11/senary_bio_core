"""
Latent space interpolation: interpolate between two sequence embeddings and
return nearest-neighbor sequence from the pool (ESM-2 has no decoder).
"""
import os
import sys
import argparse
import numpy as np
from pathlib import Path
from Bio import SeqIO

sys.path.insert(0, os.getcwd())
from modules.mining.deep_miner_utils import DeepEngine


def load_embeddings_and_ids(embed_dir: str):
    embed_dir = Path(embed_dir)
    emb = np.load(embed_dir / "embeddings.npy")
    with open(embed_dir / "sequence_ids.txt") as f:
        ids = [line.strip() for line in f if line.strip()]
    return emb, ids


def interpolate_and_nearest(embed_dir: str, fasta_path: str, id_a: str, id_b: str,
                            alphas=None, output_fasta: str = None):
    """
    Interpolate between embeddings of id_a and id_b at alphas; for each interpolated
    point find nearest sequence in pool; optionally write representative FASTA.
    """
    if alphas is None:
        alphas = [0.25, 0.5, 0.75]
    emb, ids = load_embeddings_and_ids(embed_dir)
    if id_a not in ids or id_b not in ids:
        raise ValueError(f"ids must be in {ids}")
    i_a, i_b = ids.index(id_a), ids.index(id_b)
    e_a = emb[i_a]
    e_b = emb[i_b]

    # Normalize for cosine-style nearest (optional; we use L2 here)
    results = []
    for alpha in alphas:
        e_interp = (1 - alpha) * e_a + alpha * e_b
        norms = np.linalg.norm(emb - e_interp, axis=1)
        nearest_idx = int(np.argmin(norms))
        nearest_id = ids[nearest_idx]
        results.append((alpha, nearest_id))
    if output_fasta and Path(fasta_path).exists():
        recs = {r.id: r for r in SeqIO.parse(fasta_path, "fasta")}
        with open(output_fasta, "w") as f:
            for alpha, nid in results:
                if nid in recs:
                    recs[nid].id = f"interp_alpha{alpha:.2f}_{nid}"
                    SeqIO.write(recs[nid], f, "fasta")
    return results


def main():
    parser = argparse.ArgumentParser(description="Interpolate between two sequences in ESM-2 space.")
    parser.add_argument("--embed-dir", default="data/design/embeddings", help="Embeddings from embed_pool")
    parser.add_argument("--fasta", default="data/raw_sequences/deep_hits_latest.fasta", help="Same FASTA as embed pool")
    parser.add_argument("--id-a", required=True, help="First sequence ID")
    parser.add_argument("--id-b", required=True, help="Second sequence ID")
    parser.add_argument("--alphas", type=float, nargs="+", default=[0.25, 0.5, 0.75])
    parser.add_argument("--output", default="data/design/interpolants.fasta", help="Output FASTA of nearest sequences")
    args = parser.parse_args()

    embed_dir = Path(args.embed_dir)
    if not (embed_dir / "embeddings.npy").exists():
        print("[!] Run embed_pool.py first.")
        return 1
    results = interpolate_and_nearest(
        args.embed_dir, args.fasta, args.id_a, args.id_b,
        alphas=args.alphas, output_fasta=args.output
    )
    for alpha, nid in results:
        print(f"  alpha={alpha:.2f} -> nearest: {nid}")
    print(f"[+] Wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
