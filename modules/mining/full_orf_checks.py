"""
Full-enzyme ORF checks for SRA mining: N-term start (M), C-term tail after HEPN,
and contig-boundary truncation. Used by sra_scout and autonomous_prospector.
"""
import os
import re
from typing import List

HEPN_REGEX = re.compile(r"R.{4,6}H")


def _env_bool(key: str, default: bool) -> bool:
    v = os.getenv(key, "").strip().lower()
    if not v:
        return default
    return v in ("1", "true", "yes")


def _env_int(key: str, default: int) -> int:
    v = os.getenv(key, "").strip()
    if not v:
        return default
    try:
        return int(v)
    except ValueError:
        return default


def passes_n_term(orf: str, require_m: bool = True) -> bool:
    """Require ORF to start with Met (M) for full-length (not internal fragment)."""
    if not require_m or not orf:
        return True
    return orf[0] == "M"


def passes_c_term(orf: str, min_tail: int = 15) -> bool:
    """Require min_tail aa after last HEPN motif (filters C-terminal fragments)."""
    matches = list(HEPN_REGEX.finditer(orf))
    if not matches:
        return False
    last_end = matches[-1].end()
    return len(orf) - last_end >= min_tail


def orf_nucleotide_bounds(
    frame_index: int,
    orf_index: int,
    orfs_list: List[str],
) -> tuple:
    """
    Return (start_nt, end_nt) in frame coordinates for the ORF at orf_index.
    Frame 0,1,2 = forward; 3,4,5 = reverse. Coordinates are 0-based.
    """
    frame_offset = frame_index % 3
    start_aa = sum(len(orfs_list[k]) + 1 for k in range(orf_index))
    orf_len_aa = len(orfs_list[orf_index])
    end_aa = start_aa + orf_len_aa
    start_nt = frame_offset + start_aa * 3
    end_nt = frame_offset + end_aa * 3
    return start_nt, end_nt


def is_truncated_at_boundary(
    start_nt: int,
    end_nt: int,
    contig_len: int,
    is_reverse: bool,
    margin: int = 30,
) -> bool:
    """
    True if ORF appears truncated at contig boundary (within margin nt).
    Forward: trunc 5' = start_nt < margin, trunc 3' = end_nt > contig_len - margin.
    Reverse: trunc 5' = end_nt > contig_len - 1 - margin, trunc 3' = start_nt < margin.
    """
    if is_reverse:
        trunc_5 = end_nt > contig_len - 1 - margin
        trunc_3 = start_nt < margin
    else:
        trunc_5 = start_nt < margin
        trunc_3 = end_nt > contig_len - margin
    return trunc_5 or trunc_3


def full_orf_passes(
    orf: str,
    contig_len: int,
    frame_index: int,
    orf_index: int,
    orfs_list: List[str],
    require_m: bool = True,
    min_tail: int = 15,
    boundary_margin: int = 30,
) -> bool:
    """
    Run all full-enzyme checks: N-term M, C-term tail after HEPN, not truncated at boundary.
    """
    if not passes_n_term(orf, require_m):
        return False
    if not passes_c_term(orf, min_tail):
        return False
    start_nt, end_nt = orf_nucleotide_bounds(frame_index, orf_index, orfs_list)
    is_reverse = frame_index >= 3
    if is_truncated_at_boundary(start_nt, end_nt, contig_len, is_reverse, boundary_margin):
        return False
    return True


def get_full_orf_config():
    """Read config from env for use in miners."""
    return {
        "require_m": _env_bool("REQUIRE_START_M", True),
        "min_tail": _env_int("MIN_CTERM_TAIL", 15),
        "boundary_margin": _env_int("CONTIG_BOUNDARY_MARGIN", 30),
    }
