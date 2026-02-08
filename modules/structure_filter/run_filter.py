"""
Orchestrate structure prediction and bi-lobed/HEPN filter. Reads FASTA,
runs OmegaFold, runs TM-score + HEPN check, writes passed FASTA and failed log.
"""
import os
import json
import argparse
from pathlib import Path
from Bio import SeqIO

from .predict_structures import run_omegafold
from .bi_lobed_hepn_check import run_bi_lobed_hepn_filter


def run_full_filter(
    input_fasta: str,
    output_passed_fasta: str,
    output_failed_log: str,
    structures_dir: str,
    references_dir: str,
    tm_threshold: float = 0.4,
    omegafold_repo: str = None,
    device: str = None,
) -> int:
    """
    Run OmegaFold -> bi_lobed_hepn check -> write passed FASTA and failed log.
    Returns count of passed sequences.
    """
    input_path = Path(input_fasta)
    if not input_path.exists():
        print(f"[!] Input not found: {input_fasta}")
        return 0

    # 1. Predict structures
    run_omegafold(
        str(input_path),
        structures_dir,
        omegafold_repo=omegafold_repo,
        device=device,
    )

    # 2. Run filter
    results_path = Path(structures_dir).parent / "structure_filter_results.json"
    results = run_bi_lobed_hepn_filter(
        structures_dir,
        references_dir,
        str(input_path),
        tm_threshold=tm_threshold,
        output_json=str(results_path),
    )

    # 2b. Write homology_scores.json for dashboard (per-ref: cas13a, cas13b, cas13d)
    homology_path = Path(structures_dir).parent / "homology_scores.json"
    homology = {}
    for sid, r in results.items():
        h = {}
        for k in ("cas13a", "cas13b", "cas13d"):
            if k in r and r[k] is not None:
                h[k] = r[k]
        if h:
            homology[sid] = h
    homology_path.parent.mkdir(parents=True, exist_ok=True)
    with open(homology_path, "w") as f:
        json.dump(homology, f, indent=2)

    # 3. Write passed FASTA and failed log
    seqs = {r.id: r for r in SeqIO.parse(input_path, "fasta")}
    passed_ids = [sid for sid, r in results.items() if r["pass_overall"]]
    failed = [(sid, r) for sid, r in results.items() if not r["pass_overall"]]

    Path(output_passed_fasta).parent.mkdir(parents=True, exist_ok=True)
    with open(output_passed_fasta, "w") as f:
        for sid in passed_ids:
            if sid in seqs:
                SeqIO.write(seqs[sid], f, "fasta")
    with open(output_failed_log, "w") as f:
        for sid, r in failed:
            f.write(f"{sid}\ttm_score={r['tm_score']}\thepn={r['hepn_pass']}\n")

    print(f"[+] Structure filter: {len(passed_ids)} passed, {len(failed)} failed.")
    return len(passed_ids)


def main():
    parser = argparse.ArgumentParser(description="Run full structure filter (OmegaFold + bi-lobed/HEPN).")
    parser.add_argument("--input", default="data/design/drift_variants.fasta")
    parser.add_argument("--passed-fasta", default="data/structure_pipeline/passed_structures.fasta")
    parser.add_argument("--failed-log", default="data/structure_pipeline/failed_structures.log")
    parser.add_argument("--structures-dir", default="data/structure_pipeline/structures/omegafold")
    parser.add_argument("--references-dir", default="data/structure_pipeline/references")
    parser.add_argument("--tm-threshold", type=float, default=0.4)
    parser.add_argument("--omegafold-repo", default=os.environ.get("OMEGAFOLD_REPO"))
    parser.add_argument("--device", default=None)
    args = parser.parse_args()

    n = run_full_filter(
        args.input,
        args.passed_fasta,
        args.failed_log,
        args.structures_dir,
        args.references_dir,
        tm_threshold=args.tm_threshold,
        omegafold_repo=args.omegafold_repo,
        device=args.device,
    )
    return 0 if n >= 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
