"""
Bi-lobed Cas13 + HEPN position check using TM-score vs reference and
HEPN motif positions in the predicted structure.
"""
import os
import re
import json
import subprocess
from pathlib import Path
from typing import Optional, Dict, List, Tuple

HEPN_REGEX = re.compile(r"R.{4,6}H")


def get_hepn_positions(sequence: str) -> List[Tuple[int, int]]:
    """Return (start, end) 0-based for each HEPN motif in sequence."""
    return [(m.start(), m.end()) for m in HEPN_REGEX.finditer(sequence)]


def compute_tm_score(pdb_query: str, pdb_ref: str) -> Optional[float]:
    """TM-score via US-align or tmtools."""
    try:
        result = subprocess.run(
            ["USalign", pdb_query, pdb_ref],
            capture_output=True,
            text=True,
            timeout=120,
        )
        text = result.stdout or result.stderr or ""
        m = re.search(r"TM-score\s*=\s*([\d.]+)", text)
        if m:
            return float(m.group(1))
    except (subprocess.TimeoutExpired, FileNotFoundError):
        pass
    try:
        import tmtools
        res = tmtools.tm_score(pdb_query, pdb_ref)
        return getattr(res, "tm_norm_chain1", res.tm_score)
    except Exception:
        pass
    return None


def check_hepn_distance_in_pdb(pdb_path: str, seq_id: str, sequence: str) -> bool:
    """
    Heuristic: HEPN motifs should exist in sequence; we check that we have 2-3 HEPN.
    Full 3D distance check would require parsing PDB for C-alpha of motif positions;
    for now we only require 2-3 HEPN in the sequence (structure already predicted from it).
    """
    motifs = get_hepn_positions(sequence)
    return 2 <= len(motifs) <= 3


def run_bi_lobed_hepn_filter(
    structures_dir: str,
    references_dir: str,
    fasta_path: str,
    tm_threshold: float = 0.4,
    output_json: str = None,
) -> Dict[str, dict]:
    """
    For each PDB in structures_dir, compute TM-score vs Cas13 refs; check HEPN count from FASTA.
    Returns {seq_id: {tm_score, bi_lobed_pass, hepn_pass, pass_overall}}.
    """
    from Bio import SeqIO

    ref_dir = Path(references_dir)
    struct_dir = Path(structures_dir)
    fasta_path = Path(fasta_path)

    # Reference PDBs (Cas13a, Cas13b, RfxCas13d)
    ref_pdbs = {}
    for name, pdb_id in [("cas13a", "5W1H"), ("cas13b", "6DTD"), ("cas13d", "6IV9")]:
        p = ref_dir / f"{pdb_id}.pdb"
        if p.exists():
            ref_pdbs[name] = str(p)

    if not ref_pdbs:
        print(f"[!] No reference PDBs in {references_dir}. Run run_tmscore download or add 5W1H.pdb, 6DTD.pdb, 6IV9.pdb.")

    seqs = {r.id: str(r.seq) for r in SeqIO.parse(fasta_path, "fasta")} if fasta_path.exists() else {}

    # Map PDB stem -> seq_id (OmegaFold often names by first word of header)
    query_pdbs = {}
    for pdb in struct_dir.glob("*.pdb"):
        stem = pdb.stem
        # Match to FASTA id: try stem as-is, then first part before _ or space
        for sid in seqs:
            if sid == stem or stem.startswith(sid) or sid.startswith(stem.split("_")[0]):
                query_pdbs[sid] = str(pdb.resolve())
                break
        else:
            query_pdbs[stem] = str(pdb.resolve())

    results = {}
    for seq_id, pdb_path in query_pdbs.items():
        seq = seqs.get(seq_id, "")
        hepn_pass = check_hepn_distance_in_pdb(pdb_path, seq_id, seq)
        best_tm = 0.0
        per_ref = {}
        for ref_name, ref_path in ref_pdbs.items():
            tm = compute_tm_score(pdb_path, ref_path)
            if tm is not None:
                per_ref[ref_name] = round(tm, 4)
                best_tm = max(best_tm, tm)
        bi_lobed_pass = best_tm >= tm_threshold
        pass_overall = bi_lobed_pass and hepn_pass
        results[seq_id] = {
            "tm_score": round(best_tm, 4),
            "bi_lobed_pass": bi_lobed_pass,
            "hepn_pass": hepn_pass,
            "pass_overall": pass_overall,
            **per_ref,
        }

    if output_json:
        Path(output_json).parent.mkdir(parents=True, exist_ok=True)
        with open(output_json, "w") as f:
            json.dump(results, f, indent=2)
    return results


def main():
    import argparse
    parser = argparse.ArgumentParser(description="Bi-lobed + HEPN structure filter.")
    parser.add_argument("--structures-dir", default="data/structure_pipeline/structures/omegafold")
    parser.add_argument("--references-dir", default="data/structure_pipeline/references")
    parser.add_argument("--fasta", default="data/design/drift_variants.fasta")
    parser.add_argument("--tm-threshold", type=float, default=0.4)
    parser.add_argument("--output", default="data/structure_pipeline/structure_filter_results.json")
    args = parser.parse_args()

    results = run_bi_lobed_hepn_filter(
        args.structures_dir,
        args.references_dir,
        args.fasta,
        tm_threshold=args.tm_threshold,
        output_json=args.output,
    )
    passed = sum(1 for r in results.values() if r["pass_overall"])
    print(f"[+] {passed}/{len(results)} passed structure filter. Results: {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
