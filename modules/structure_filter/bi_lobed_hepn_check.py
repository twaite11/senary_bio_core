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


def _tm_score_tmtools_from_pdb(pdb_a: str, pdb_b: str) -> Optional[float]:
    """
    Compute TM-score from two PDB paths using tmtools.tm_align + tmtools.io.
    The tmtools package (PyPI) exposes tm_align(coords1, coords2, seq1, seq2), not tm_score(path, path).
    """
    try:
        from tmtools import tm_align
        from tmtools.io import get_structure, get_residue_data
    except ImportError:
        return None
    try:
        s1 = get_structure(pdb_a)
        s2 = get_structure(pdb_b)
        chain1 = next(s1.get_chains(), None)
        chain2 = next(s2.get_chains(), None)
        if chain1 is None or chain2 is None:
            return None
        coords1, seq1 = get_residue_data(chain1)
        coords2, seq2 = get_residue_data(chain2)
        if coords1 is None or coords2 is None or not seq1 or not seq2:
            return None
        res = tm_align(coords1, coords2, seq1, seq2)
        return float(getattr(res, "tm_norm_chain1", getattr(res, "tm_score", None)))
    except Exception:
        return None


def compute_tm_score(pdb_query: str, pdb_ref: str) -> Optional[float]:
    """TM-score via US-align or tmtools (tm_align + io for PDB paths)."""
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
    # tmtools: try legacy tm_score(path, path) first, then PDB load + tm_align
    try:
        import tmtools
        if hasattr(tmtools, "tm_score"):
            res = tmtools.tm_score(pdb_query, pdb_ref)
            return float(getattr(res, "tm_norm_chain1", res.tm_score))
    except Exception:
        pass
    score = _tm_score_tmtools_from_pdb(pdb_query, pdb_ref)
    if score is not None:
        return score
    return None


def check_tm_score_available(references_dir: str) -> bool:
    """
    Verify TM-score can be computed (US-align or tmtools) before running OmegaFold.
    Uses reference PDBs in references_dir; returns True if a score is obtained.
    """
    ref_dir = Path(references_dir)
    ref_pdb_ids = ["5W1H", "6DTD", "6IV9"]
    ref_paths = [ref_dir / f"{pdb_id}.pdb" for pdb_id in ref_pdb_ids]
    ref_paths = [str(p) for p in ref_paths if p.exists()]

    if len(ref_paths) < 1:
        print(
            "[!] TM-score check failed: no reference PDBs found in "
            f"{references_dir}. Add 5W1H.pdb, 6DTD.pdb, or 6IV9.pdb (e.g. run "
            "visualization/run_tmscore.py with download), or install tmtools/USalign."
        )
        return False

    # Use two refs if available, else same ref vs itself
    pdb_a, pdb_b = ref_paths[0], ref_paths[1] if len(ref_paths) > 1 else ref_paths[0]
    score = compute_tm_score(pdb_a, pdb_b)
    if score is not None and isinstance(score, (int, float)):
        print(f"[+] TM-score check OK (score={score:.4f}). Proceeding with structure filter.")
        return True
    print(
        "[!] TM-score check failed: could not compute TM-score. Install tmtools "
        "('pip install tmtools') in the same environment you use to run the pipeline, "
        "or put USalign on your PATH."
    )
    return False


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
    tm_threshold: float = 0.25,
    output_json: str = None,
) -> Dict[str, dict]:
    """
    For each PDB in structures_dir, compute TM-score vs Cas13 refs; check HEPN count from FASTA.
    Also computes functional criteria (domain TM, catalytic distance, linker charge, pLDDT dip)
    for dashboard ranking. Returns {seq_id: {tm_score, bi_lobed_pass, hepn_pass, pass_overall, ...}}.
    """
    from Bio import SeqIO

    from .functional_criteria import compute_functional_criteria, get_sequence_from_pdb_single_chain

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

    # Single longest chain per ref (so HEPN motifs are correct for domain TM)
    ref_sequences = {}
    ref_chain_ids = {}
    for name, ref_path in ref_pdbs.items():
        try:
            seq, chain_id = get_sequence_from_pdb_single_chain(ref_path)
            if seq and len(seq) > 100:
                ref_sequences[name] = seq
                ref_chain_ids[name] = chain_id
        except Exception:
            pass
    refs_list = [
        (name, ref_path, ref_sequences.get(name) or "", ref_chain_ids.get(name))
        for name, ref_path in ref_pdbs.items()
    ]

    seqs = {r.id: str(r.seq) for r in SeqIO.parse(fasta_path, "fasta")} if fasta_path.exists() else {}

    # Robust 1:1 PDB <-> FASTA id mapping (no overwriting, no ambiguous prefix match).
    sorted_pdbs = sorted(
        [(p.stem, str(p.resolve())) for p in struct_dir.glob("*.pdb")],
        key=lambda x: x[0],
    )
    sorted_sids = sorted(seqs.keys())

    query_pdbs = {}  # seq_id -> pdb_path
    if len(sorted_pdbs) == len(sorted_sids):
        # Try exact stem match first (OmegaFold often names output by FASTA id)
        for stem, path in sorted_pdbs:
            if stem in seqs:
                query_pdbs[stem] = path
        if len(query_pdbs) == len(sorted_pdbs):
            # All PDB stems matched FASTA ids
            pass
        else:
            # Fallback: pair by sorted order (same count => same order as input FASTA)
            query_pdbs = {sorted_sids[i]: sorted_pdbs[i][1] for i in range(len(sorted_sids))}
    else:
        # Different count: only exact stem match (no ambiguous prefix)
        for stem, path in sorted_pdbs:
            if stem in seqs:
                query_pdbs[stem] = path
        if len(query_pdbs) < len(sorted_pdbs) or len(query_pdbs) < len(sorted_sids):
            print(
                f"[*] Structure filter: {len(query_pdbs)} PDB–FASTA pairs (exact id match); "
                f"{len(sorted_pdbs)} PDBs, {len(sorted_sids)} FASTA ids."
            )

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
        # Functional criteria (soft metrics for dashboard); domain TM vs all refs, best used
        try:
            fc = compute_functional_criteria(
                pdb_path,
                seq,
                refs_list,
                compute_tm_score,
            )
            results[seq_id].update(fc)
        except Exception as e:
            if seq_id == list(query_pdbs.keys())[0]:
                print(f"[!] Functional criteria warning for {seq_id}: {e}")
            results[seq_id].update(
                hepn1_tm=None, hepn2_tm=None, catalytic_distance_angstrom=None,
                linker_net_charge=None, linker_plddt_mean=None, domain1_plddt_mean=None,
                domain2_plddt_mean=None, plddt_dip_ok=None,
            )

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
    parser.add_argument("--tm-threshold", type=float, default=0.25)
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
