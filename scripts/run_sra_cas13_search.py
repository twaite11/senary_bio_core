#!/usr/bin/env python3
"""
Run SRA Cas13 search: Diamond (bait) -> Hook (pairs) -> Megahit (assemble) -> Prodigal + filter.
Fetches Run accessions (SRR) by ESearch, downloads reads via fasterq-dump, runs Diamond blastx,
extracts hit + mate reads, assembles with Megahit, annotates with Prodigal, filters to full
Cas13-like genes (700-1400 aa, 2 HEPN, N/C intact, CRISPR on contig). Saves deep_hits_*.fasta
and metadata for the structure pipeline.

Requirements (on PATH):
  - SRA Toolkit (fasterq-dump)
  - Diamond, Megahit, Prodigal
  - pip: biopython

Usage:
  python scripts/run_sra_cas13_search.py
  python scripts/run_sra_cas13_search.py --max-runs 10 --max-reads 500000
  SRA_TERM="txid2[ORGN] AND metagenomic" python scripts/run_sra_cas13_search.py

Environment:
  SRA_TERM            - NCBI SRA search query (default: Bacteria|Archaea, WGS|METAGENOMIC)
  MAX_SRA_RUNS        - Max Run accessions to process (default 100000)
  REFERENCE_FASTA     - Protein Cas13 reference (default data/references/mining_refs.fasta)
  OUTPUT_DIR          - Where to write deep_hits_*.fasta (default data/raw_sequences)
  SRA_TEMP_DIR        - Temp dir for SRA dumps and per-run files (default system temp)
  DIAMOND_CMD         - diamond executable (default diamond)
  FASTERQ_DUMP_CMD    - fasterq-dump executable (default fasterq-dump)
  MEGAHIT_CMD         - megahit executable (default megahit)
  PRODIGAL_CMD        - prodigal executable (default prodigal)
  MAX_READS_PER_RUN   - Cap reads per run for bait step (default: no cap)
  DIAMOND_SENSITIVITY - sensitive|more-sensitive|very-sensitive (default sensitive)
  NUM_THREADS         - Threads per tool (default 4)
  SRA_WORKERS, NUM_THREADS - workers and threads per worker (default 4, 8).
"""
from __future__ import annotations

import os
import sys
import argparse

# Project root
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from modules.mining.sra_blast_miner import (
    DEFAULT_SRA_TERM,
    fetch_sra_run_accessions,
    mine_sra,
)


def main():
    parser = argparse.ArgumentParser(
        description="SRA Cas13 search: Diamond -> Hook -> Megahit -> Prodigal. Output: deep_hits_*.fasta for structure pipeline."
    )
    parser.add_argument(
        "--sra-term",
        default=os.environ.get("SRA_TERM", DEFAULT_SRA_TERM),
        help="NCBI SRA search query (organism, strategy, library)",
    )
    parser.add_argument(
        "--max-runs",
        type=int,
        default=int(os.environ.get("MAX_SRA_RUNS", "100000")),
        help="Max SRA Run accessions to process",
    )
    parser.add_argument(
        "--reference-fasta",
        default=os.environ.get("REFERENCE_FASTA", os.path.join(ROOT, "data", "references", "mining_refs.fasta")),
        help="Protein Cas13 reference for Diamond DB",
    )
    parser.add_argument(
        "--output-dir",
        default=os.environ.get("OUTPUT_DIR", os.path.join(ROOT, "data", "raw_sequences")),
        help="Output directory for deep_hits_*.fasta and metadata",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=int(os.environ.get("SRA_WORKERS", "4")),
        help="Parallel SRA runs (total CPU ≈ workers × num-threads)",
    )
    parser.add_argument(
        "--num-threads",
        type=int,
        default=int(os.environ.get("NUM_THREADS", "8")),
        help="Threads per worker (Diamond, Megahit)",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Only fetch and print SRA Run list; do not run pipeline",
    )
    args = parser.parse_args()

    print("[*] SRA Cas13 search (spot-check -> Diamond -> Hook -> Megahit -> Prodigal)")
    print(f"    SRA term: {args.sra_term[:80]}...")
    print(f"    Max runs: {args.max_runs}, workers: {args.workers}, threads/worker: {args.num_threads}")

    print("[*] Fetching SRA Run accessions (paginated)...")
    runs = fetch_sra_run_accessions(args.sra_term, max_records=args.max_runs, page_size=500)
    print(f"    Got {len(runs)} Run accessions (SRR)")

    if not runs:
        print("[!] No SRA runs found. Check SRA_TERM or NCBI availability.")
        return 1

    if args.dry_run:
        for i, r in enumerate(runs[:20]):
            print(f"      {r}")
        if len(runs) > 20:
            print(f"      ... and {len(runs) - 20} more")
        return 0

    if not os.path.exists(args.reference_fasta):
        print(f"[!] Reference FASTA not found: {args.reference_fasta}")
        return 1

    print("[*] Running miner (spot-check -> fetch -> Diamond -> Hook -> Megahit -> Prodigal)...")
    mine_sra(
        runs,
        args.reference_fasta,
        args.output_dir,
        workers=args.workers,
        threads_per_worker=args.num_threads,
    )
    print("[*] Done. Check output dir for deep_hits_*.fasta and deep_hits_*_metadata.csv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
