"""
Functional criteria for Cas13d-like nucleases (soft metrics for dashboard ranking).
1. Domain TM-score: HEPN1 and HEPN2 local TM vs references (5W1H, 6DTD, 6IV9); single longest
   chain per ref; best hepn1_tm and best hepn2_tm across refs are used.
2. Catalytic distance: C-alpha distance between the two RxxxxH histidines (target < 25 Å).
3. Positive groove proxy: net charge (R+K-D-E) in linker between HEPN domains.
4. Linker pLDDT dip: high confidence in domains, lower in linker (flexibility signal).
"""
import re
import tempfile
from pathlib import Path
from typing import Dict, List, Optional, Tuple

HEPN_REGEX = re.compile(r"R.{4,6}H")

# 3-letter to 1-letter amino acid
RESIDUE_1L = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V", "UNK": "X",
}


def get_hepn_positions(sequence: str) -> List[Tuple[int, int]]:
    """Return (start, end) 0-based for each HEPN motif in sequence."""
    return [(m.start(), m.end()) for m in HEPN_REGEX.finditer(sequence)]


def _parse_pdb_residues(pdb_path: str) -> List[Tuple[str, int, float, float, float, float]]:
    """
    Parse PDB ATOM lines: list of (resname_3, resseq_1based, x, y, z, bfactor).
    Uses CA only; residues in order of appearance.
    """
    out = []
    seen = set()
    with open(pdb_path) as f:
        for line in f:
            if not line.startswith("ATOM "):
                continue
            if line[12:16].strip() != "CA":
                continue
            resname = line[17:20].strip()
            try:
                resseq = int(line[22:26].strip())
            except ValueError:
                continue
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            try:
                bfactor = float(line[60:66])
            except ValueError:
                bfactor = 0.0
            if resseq not in seen:
                seen.add(resseq)
                out.append((resname, resseq, x, y, z, bfactor))
    out.sort(key=lambda r: r[1])
    return out


def _parse_pdb_by_chains(pdb_path: str) -> Dict[str, List[Tuple[str, int, float, float, float, float]]]:
    """
    Parse PDB ATOM lines by chain: chain_id -> list of (resname, resseq, x, y, z, bfactor).
    CA only; each chain's list sorted by resseq. PDB chain ID at column 22 (index 21).
    """
    by_chain: Dict[str, List[Tuple]] = {}
    with open(pdb_path) as f:
        for line in f:
            if not line.startswith("ATOM "):
                continue
            if line[12:16].strip() != "CA":
                continue
            chain_id = line[21] if len(line) > 21 else "A"
            resname = line[17:20].strip()
            try:
                resseq = int(line[22:26].strip())
            except ValueError:
                continue
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            try:
                bfactor = float(line[60:66])
            except ValueError:
                bfactor = 0.0
            key = (chain_id, resseq)
            if chain_id not in by_chain:
                by_chain[chain_id] = []
            # avoid duplicate (chain, resseq) e.g. alt locs
            existing_resseqs = {r[1] for r in by_chain[chain_id]}
            if resseq not in existing_resseqs:
                by_chain[chain_id].append((resname, resseq, x, y, z, bfactor))
    for ch in by_chain:
        by_chain[ch].sort(key=lambda r: r[1])
    return by_chain


def get_sequence_from_pdb_single_chain(
    pdb_path: str,
) -> Tuple[Optional[str], Optional[str]]:
    """
    Extract sequence from the single longest chain (by residue count) in the PDB.
    Returns (sequence_str, chain_id) or (None, None) if no valid chain.
    """
    by_chain = _parse_pdb_by_chains(pdb_path)
    if not by_chain:
        return None, None
    best_chain = max(by_chain.keys(), key=lambda c: len(by_chain[c]))
    residues = by_chain[best_chain]
    seq = "".join(RESIDUE_1L.get(r[0], "X") for r in residues)
    return seq, best_chain


def get_sequence_from_pdb(pdb_path: str) -> str:
    """Extract one-letter sequence from PDB (CA atoms, ordered by resseq)."""
    residues = _parse_pdb_residues(pdb_path)
    return "".join(RESIDUE_1L.get(r[0], "X") for r in residues)


def get_residue_coords_bfactor(pdb_path: str) -> List[Tuple[float, float, float, float]]:
    """
    Return list of (x, y, z, bfactor) for each residue (1-based index = position in list).
    Length = number of residues. Index i (0-based) = residue i+1 in PDB order.
    """
    residues = _parse_pdb_residues(pdb_path)
    return [(r[2], r[3], r[4], r[5]) for r in residues]


def _write_pdb_slice(
    pdb_path: str,
    start_0: int,
    end_0: int,
    out_path: str,
    chain_id: Optional[str] = None,
) -> bool:
    """
    Write a PDB containing only residues [start_0, end_0] (0-based, inclusive).
    If chain_id is set, use only that chain's residues (for multi-chain refs).
    """
    if chain_id is not None:
        by_chain = _parse_pdb_by_chains(pdb_path)
        residues = by_chain.get(chain_id, [])
    else:
        residues = _parse_pdb_residues(pdb_path)
    if not residues or end_0 >= len(residues) or start_0 < 0:
        return False
    resseq_set = set()
    for i in range(start_0, min(end_0 + 1, len(residues))):
        resseq_set.add(residues[i][1])
    written = []
    with open(pdb_path) as f:
        for line in f:
            if not line.startswith("ATOM "):
                continue
            if chain_id is not None and (len(line) <= 21 or line[21] != chain_id):
                continue
            try:
                resseq = int(line[22:26].strip())
            except ValueError:
                continue
            if resseq in resseq_set:
                written.append(line)
    if not written:
        return False
    with open(out_path, "w") as f:
        f.writelines(written)
    return True


def _domain_spans(sequence: str, padding_before: int = 40, padding_after: int = 60) -> List[Tuple[int, int]]:
    """Return (start_0, end_0) inclusive for each HEPN domain (motif ± padding)."""
    motifs = get_hepn_positions(sequence)
    n = len(sequence)
    spans = []
    for start, end in motifs:
        s = max(0, start - padding_before)
        e = min(n - 1, end - 1 + padding_after)
        spans.append((s, e))
    return spans


def compute_domain_tm_scores(
    query_pdb: str,
    query_sequence: str,
    ref_pdb: str,
    ref_sequence: str,
    compute_tm_score_fn,
    ref_chain_id: Optional[str] = None,
) -> Tuple[Optional[float], Optional[float]]:
    """
    Compute TM-score of query HEPN1 vs ref HEPN1 and query HEPN2 vs ref HEPN2.
    Returns (hepn1_tm, hepn2_tm). If ref_chain_id is set, only that chain is used in ref.
    """
    q_spans = _domain_spans(query_sequence)
    r_spans = _domain_spans(ref_sequence)
    if len(q_spans) < 2 or len(r_spans) < 2:
        return None, None
    scores = []
    for i in (0, 1):
        qs, qe = q_spans[i]
        rs, re = r_spans[i]
        with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False) as fq:
            with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False) as fr:
                try:
                    _write_pdb_slice(query_pdb, qs, qe, fq.name)
                    _write_pdb_slice(ref_pdb, rs, re, fr.name, chain_id=ref_chain_id)
                    tm = compute_tm_score_fn(fq.name, fr.name)
                    scores.append(round(tm, 4) if tm is not None else None)
                finally:
                    Path(fq.name).unlink(missing_ok=True)
                    Path(fr.name).unlink(missing_ok=True)
    return (scores[0] if len(scores) > 0 else None), (scores[1] if len(scores) > 1 else None)


def catalytic_distance_angstrom(pdb_path: str, sequence: str) -> Optional[float]:
    """
    C-alpha distance between the histidine of the first RxxxxH and the histidine of the second RxxxxH.
    Returns distance in Å or None if not exactly two HEPN motifs.
    """
    motifs = get_hepn_positions(sequence)
    if len(motifs) < 2:
        return None
    # Histidine is the last residue of the motif (end-1 in 0-based)
    his1 = motifs[0][1] - 1  # 0-based index of His in first motif
    his2 = motifs[1][1] - 1  # 0-based index of His in second motif
    coords_bf = get_residue_coords_bfactor(pdb_path)
    if his1 >= len(coords_bf) or his2 >= len(coords_bf):
        return None
    x1, y1, z1, _ = coords_bf[his1]
    x2, y2, z2, _ = coords_bf[his2]
    return ((x1 - x2) ** 2 + (y1 - y2) ** 2 + (z1 - z2) ** 2) ** 0.5


def linker_net_charge(sequence: str) -> Optional[int]:
    """
    Net charge (R+K - D-E) for the linker region between first HEPN end and second HEPN start.
    Positive is favorable for RNA binding. Returns None if < 2 HEPN motifs.
    """
    motifs = get_hepn_positions(sequence)
    if len(motifs) < 2:
        return None
    start = motifs[0][1]  # first residue after first motif
    end = motifs[1][0]    # first residue of second motif
    if start >= end:
        return 0
    segment = sequence[start:end]
    pos = segment.count("R") + segment.count("K")
    neg = segment.count("D") + segment.count("E")
    return pos - neg


def plddt_metrics(pdb_path: str, sequence: str) -> Dict[str, Optional[float]]:
    """
    pLDDT from B-factor (OmegaFold/AlphaFold use 0-100 in B-factor).
    Returns linker_plddt_mean, domain1_plddt_mean, domain2_plddt_mean, plddt_dip_ok.
    dip_ok = linker mean < 60 and both domain means > 80 (flexible linker, confident domains).
    """
    coords_bf = get_residue_coords_bfactor(pdb_path)
    if len(coords_bf) != len(sequence):
        return {
            "linker_plddt_mean": None,
            "domain1_plddt_mean": None,
            "domain2_plddt_mean": None,
            "plddt_dip_ok": None,
        }
    bfactors = [c[3] for c in coords_bf]
    motifs = get_hepn_positions(sequence)
    if len(motifs) < 2:
        return {
            "linker_plddt_mean": None,
            "domain1_plddt_mean": None,
            "domain2_plddt_mean": None,
            "plddt_dip_ok": None,
        }
    pad = 45
    # Domain 1: around first motif
    d1_start = max(0, motifs[0][0] - pad)
    d1_end = min(len(bfactors), motifs[0][1] + pad)
    d1_vals = bfactors[d1_start:d1_end]
    # Domain 2: around second motif
    d2_start = max(0, motifs[1][0] - pad)
    d2_end = min(len(bfactors), motifs[1][1] + pad)
    d2_vals = bfactors[d2_start:d2_end]
    # Linker: between first motif end and second motif start
    link_start = motifs[0][1]
    link_end = motifs[1][0]
    if link_end <= link_start:
        link_vals = []
    else:
        link_vals = bfactors[link_start:link_end]

    def mean_opt(vals):
        return round(sum(vals) / len(vals), 2) if vals else None

    linker_mean = mean_opt(link_vals)
    dom1_mean = mean_opt(d1_vals)
    dom2_mean = mean_opt(d2_vals)
    dip_ok = (
        linker_mean is not None and dom1_mean is not None and dom2_mean is not None
        and linker_mean < 60 and dom1_mean > 80 and dom2_mean > 80
    )
    return {
        "linker_plddt_mean": linker_mean,
        "domain1_plddt_mean": dom1_mean,
        "domain2_plddt_mean": dom2_mean,
        "plddt_dip_ok": dip_ok,
    }


def compute_functional_criteria(
    pdb_path: str,
    sequence: str,
    refs_list: List[Tuple[str, str, str, Optional[str]]],
    compute_tm_score_fn,
) -> Dict[str, Optional[float]]:
    """
    Compute all four functional criteria for one structure.
    refs_list: [(ref_name, ref_path, ref_sequence, ref_chain_id), ...] for each reference.
    Domain TM is computed vs each ref (single chain); best hepn1_tm and best hepn2_tm are kept.
    """
    out = {
        "hepn1_tm": None,
        "hepn2_tm": None,
        "catalytic_distance_angstrom": None,
        "linker_net_charge": None,
        "linker_plddt_mean": None,
        "domain1_plddt_mean": None,
        "domain2_plddt_mean": None,
        "plddt_dip_ok": None,
    }
    out["catalytic_distance_angstrom"] = catalytic_distance_angstrom(pdb_path, sequence)
    out["linker_net_charge"] = linker_net_charge(sequence)
    plddt = plddt_metrics(pdb_path, sequence)
    out.update(plddt)

    best_hepn1, best_hepn2 = None, None
    for _name, ref_path, ref_sequence, ref_chain_id in refs_list:
        if not ref_path or not Path(ref_path).exists() or not ref_sequence:
            continue
        if len(_domain_spans(ref_sequence)) < 2:
            continue
        hepn1_tm, hepn2_tm = compute_domain_tm_scores(
            pdb_path,
            sequence,
            ref_path,
            ref_sequence,
            compute_tm_score_fn,
            ref_chain_id=ref_chain_id,
        )
        if hepn1_tm is not None and (best_hepn1 is None or hepn1_tm > best_hepn1):
            best_hepn1 = hepn1_tm
        if hepn2_tm is not None and (best_hepn2 is None or hepn2_tm > best_hepn2):
            best_hepn2 = hepn2_tm
    out["hepn1_tm"] = best_hepn1
    out["hepn2_tm"] = best_hepn2
    return out


def functional_score(criteria: Dict) -> Tuple[float, float]:
    """
    Composite score for ranking: (discrete_pass_count, continuous_tiebreaker).
    Higher is better. discrete = number of criteria passed (max 4); tiebreaker for same count.
    """
    pass_count = 0.0
    tiebreaker = 0.0
    # 1. HEPN1 TM >= 0.7
    if criteria.get("hepn1_tm") is not None:
        if criteria["hepn1_tm"] >= 0.7:
            pass_count += 1
        tiebreaker += criteria["hepn1_tm"] or 0
    # 2. HEPN2 TM >= 0.7
    if criteria.get("hepn2_tm") is not None:
        if criteria["hepn2_tm"] >= 0.7:
            pass_count += 1
        tiebreaker += criteria["hepn2_tm"] or 0
    # 3. Catalytic distance < 25 Å
    d = criteria.get("catalytic_distance_angstrom")
    if d is not None:
        if d < 25:
            pass_count += 1
        tiebreaker += max(0, 25 - d) / 25.0  # closer is better
    # 4. pLDDT dip (flexible linker)
    if criteria.get("plddt_dip_ok") is True:
        pass_count += 1
        tiebreaker += 1
    # Positive groove: use as tiebreaker only (linker_net_charge > 0)
    charge = criteria.get("linker_net_charge")
    if charge is not None and charge > 0:
        tiebreaker += 0.5
    return (pass_count, tiebreaker)
