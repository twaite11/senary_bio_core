"""
Propose point mutations, score with ESM-2 (stability), and keep variants
that stay below a max identity vs reference (drift goal).
"""
import os
import sys
import argparse
import random
from pathlib import Path
from Bio import SeqIO, Align
from Bio.Align import substitution_matrices

sys.path.insert(0, os.getcwd())
from modules.mining.deep_miner_utils import DeepEngine

# Standard amino acids
AA = list("ACDEFGHIKLMNPQRSTVWY")


def pairwise_identity(seq1: str, seq2: str) -> float:
    """Global alignment identity (fraction of aligned positions that match)."""
    aligner = Align.PairwiseAligner(mode="global", substitution_matrix=substitution_matrices.load("BLOSUM62"))
    alns = aligner.align(seq1, seq2)
    if not alns:
        return 0.0
    a = alns[0]
    matches = sum(1 for i, j in zip(a[0], a[1]) if i == j and i != "-")
    length = max(len(seq1), len(seq2))
    return matches / length if length else 0.0


def max_identity_to_refs(seq: str, ref_seqs: list) -> float:
    return max(pairwise_identity(seq, r) for r in ref_seqs) if ref_seqs else 0.0


def propose_mutations(sequence: str, num_mutants: int = 5, num_sites: int = 3):
    """Propose num_mutants variants, each with up to num_sites random AA substitutions."""
    seq = list(sequence)
    L = len(seq)
    mutants = []
    seen = set()
    for _ in range(num_mutants * 4):  # allow retries
        if len(mutants) >= num_mutants:
            break
        positions = random.sample(range(L), min(num_sites, L)) if L else []
        mut = seq[:]
        desc = []
        for pos in positions:
            wild = mut[pos]
            choices = [a for a in AA if a != wild]
            if not choices:
                continue
            mut[pos] = random.choice(choices)
            desc.append(f"{pos}{wild}{mut[pos]}")
        mut_str = "".join(mut)
        key = (tuple(positions), mut_str)
        if key not in seen:
            seen.add(key)
            mutants.append((mut_str, "_".join(desc)))
    return mutants


def main():
    parser = argparse.ArgumentParser(description="Mutate for drift: ESM-2 stability + identity < threshold.")
    parser.add_argument("--input", default="data/raw_sequences/deep_hits_latest.fasta", help="Input FASTA")
    parser.add_argument("--references", default="data/known_cas13.fasta",
                        help="FASTA of known Cas13 sequences (Lwa, Rfx, etc.) for identity check")
    parser.add_argument("--output", default="data/design/drift_variants.fasta", help="Output FASTA")
    parser.add_argument("--max-identity", type=float, default=0.85, help="Max identity to any reference (drift goal)")
    parser.add_argument("--min-esm-score", type=float, default=0.5, help="Min ESM-2 similarity to keep (stability proxy)")
    parser.add_argument("--num-mutants-per-seq", type=int, default=5, help="Max mutants to try per input sequence")
    parser.add_argument("--num-mutations", type=int, default=3, help="Number of positions to mutate per variant")
    args = parser.parse_args()

    input_path = Path(args.input)
    ref_path = Path(args.references)
    if not input_path.exists():
        print(f"[!] Input not found: {input_path}")
        return 1

    ref_seqs = []
    if ref_path.exists():
        ref_seqs = [str(r.seq) for r in SeqIO.parse(ref_path, "fasta")]
    else:
        print(f"[!] No reference FASTA at {ref_path}; identity check will be 0.")

    engine = DeepEngine()
    records = list(SeqIO.parse(input_path, "fasta"))
    kept = []
    for rec in records:
        seq = str(rec.seq)
        if len(seq) < 300:
            continue
        base_identity = max_identity_to_refs(seq, ref_seqs)
        mutants = propose_mutations(seq, num_mutants=args.num_mutants_per_seq, num_sites=args.num_mutations)
        added_any = False
        for mut_seq, desc in mutants:
            score = engine.score_candidate(mut_seq)
            ident = max_identity_to_refs(mut_seq, ref_seqs)
            if score >= args.min_esm_score and ident < args.max_identity:
                mut_id = f"{rec.id}_mut_{desc}_id{ident:.2f}"
                kept.append((mut_id, mut_seq, score, ident))
                added_any = True
        if not added_any and base_identity < args.max_identity:
            kept.append((rec.id, seq, engine.score_candidate(seq), base_identity))

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w") as f:
        for mid, mseq, score, ident in kept:
            f.write(f">{mid}_esm{score:.3f}\n{mseq}\n")
    print(f"[+] Wrote {len(kept)} drift variants to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
