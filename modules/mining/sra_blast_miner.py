"""
SRA Cas13 search using Magic-BLAST (blastn_vdb) against SRA runs without downloading.
Fetches Run accessions (SRR) by ESearch, runs magic-blast -sra, parses hits, applies
full-enzyme filters (700-1400 aa, exactly 2 HEPN, N-term M, C-term tail, mandatory CRISPR).
CRISPR domain sequences are extracted and saved in metadata with SRA accession for crRNA binding.
Saves deep_hits FASTA + metadata for the structure pipeline.

Pagination: SRA ESearch is paginated (retstart/retmax). We keep requesting pages until
we have max_total runs or NCBI returns no more. There is no fixed "page limit"—we stop
when the server returns fewer than page_size results or we hit our max_total cap.
"""
from __future__ import annotations

import os
import re
import subprocess
import tempfile
from pathlib import Path
from datetime import datetime
from typing import List, Tuple, Optional

from Bio import Entrez, SeqIO, Seq
from Bio.Seq import Seq

# Project root
_root = Path(__file__).resolve().parents[2]
if str(_root) not in __import__("sys").path:
    __import__("sys").path.insert(0, str(_root))

try:
    from modules.mining.full_orf_checks import passes_n_term, passes_c_term, get_full_orf_config
except ImportError:
    from full_orf_checks import passes_n_term, passes_c_term, get_full_orf_config

Entrez.email = os.environ.get("ENTREZ_EMAIL", "founder@senarybio.com")

HEPN_REGEX = re.compile(r"R.{4,6}H")

# Default SRA search: Bacteria OR Archaea, WGS or METAGENOMIC strategy, METAGENOMIC library source
DEFAULT_SRA_TERM = (
    '(txid2[ORGN] OR txid2157[ORGN]) '
    'AND (wgs[Strategy] OR metagenomic[Strategy]) '
    'AND metagenomic[LibrarySource]'
)

MIN_AA = 700
MAX_AA = 1400
HEPN_COUNT_EXACT = 2


# E. coli / bacterial preferred codons for back-translation (table 11)
_BACK_TABLE = {
    "M": "ATG", "F": "TTT", "L": "CTG", "S": "TCT", "Y": "TAT", "*": "TAA",
    "C": "TGC", "W": "TGG", "P": "CCG", "H": "CAT", "Q": "CAG", "R": "CGT",
    "I": "ATT", "T": "ACT", "N": "AAT", "K": "AAA", "V": "GTG", "A": "GCG",
    "D": "GAT", "E": "GAA", "G": "GGT",
}


def _back_translate(protein: str) -> str:
    """Back-translate protein to nucleotide using bacterial-preferred codons."""
    dna = []
    for aa in protein.upper():
        if aa == "*":
            break
        dna.append(_BACK_TABLE.get(aa, "ATG"))
    return "".join(dna)


def _parse_runs_from_xml(xml_handle) -> List[str]:
    """Parse RUN accessions from SRA EFetch XML."""
    import xml.etree.ElementTree as ET
    runs = []
    try:
        tree = ET.parse(xml_handle)
        root = tree.getroot()
        for elem in root.iter():
            if elem.tag.endswith("RUN") or (elem.tag == "RUN"):
                acc = elem.get("accession") or elem.get("acc")
                if acc and str(acc).startswith("SRR"):
                    runs.append(str(acc))
            if elem.tag.endswith("Run") and elem.text and elem.text.strip().startswith("SRR"):
                runs.append(elem.text.strip())
    except Exception:
        pass
    return runs


def fetch_sra_run_accessions(
    term: str = DEFAULT_SRA_TERM,
    max_records: int = 100_000,
    batch_size: int = 500,
    retstart: int = 0,
) -> List[str]:
    """
    ESearch SRA with term, then EFetch XML to resolve Run accessions (SRR).
    Returns list of SRR accessions. Paginate with retstart for large result sets.
    """
    run_ids = []
    try:
        h = Entrez.esearch(db="sra", term=term, retmax=min(batch_size, max_records), retstart=retstart)
        rec = Entrez.read(h)
        h.close()
        id_list = rec.get("IdList", [])
        if not id_list:
            return run_ids
        # EFetch XML to get RUN accessions (experiments contain one or more runs)
        h = Entrez.efetch(db="sra", id=",".join(id_list), retmode="xml")
        run_ids = _parse_runs_from_xml(h)
        h.close()
        if not run_ids:
            # Fallback: ESummary - some responses include Run list
            h = Entrez.esummary(db="sra", id=",".join(id_list))
            summary = Entrez.read(h)
            h.close()
            for item in (summary if isinstance(summary, list) else [summary]):
                if isinstance(item, dict):
                    for key in ("Runs", "Run", "runs", "run"):
                        if key in item:
                            val = item[key]
                            if isinstance(val, str) and val.startswith("SRR"):
                                run_ids.append(val)
                            elif isinstance(val, list):
                                for r in val:
                                    if isinstance(r, str) and r.startswith("SRR"):
                                        run_ids.append(r)
                                    elif isinstance(r, dict) and "acc" in r:
                                        run_ids.append(r["acc"])
                    acc = item.get("Accession") or item.get("Run") or item.get("acc")
                    if acc and str(acc).startswith("SRR"):
                        run_ids.append(str(acc))
    except Exception as e:
        print(f"[!] SRA fetch error: {e}")
    return list(dict.fromkeys(run_ids))


def get_all_sra_runs(
    term: str = DEFAULT_SRA_TERM,
    max_total: int = 100_000,
    page_size: int = 500,
) -> List[str]:
    """
    Paginate SRA search until we have up to max_total Run accessions.
    Keeps requesting pages (retstart) until NCBI returns no more runs or we hit max_total.
    There is no fixed page limit—stops only when the server returns fewer than page_size
    results or we have collected max_total.
    """
    all_runs = []
    retstart = 0
    while len(all_runs) < max_total:
        batch = fetch_sra_run_accessions(term, max_records=page_size, batch_size=page_size, retstart=retstart)
        if not batch:
            break
        for r in batch:
            if r not in all_runs:
                all_runs.append(r)
        if len(batch) < page_size:
            break
        retstart += page_size
        if retstart >= max_total:
            break
    return all_runs[:max_total]


def build_nucleotide_reference(protein_fasta: str, max_aa: int = 600, out_path: Optional[str] = None) -> str:
    """
    Back-translate first max_aa of first protein in FASTA to DNA; write to out_path or temp.
    Returns path to nucleotide FASTA for magic-blast -query.
    """
    path = Path(protein_fasta)
    if not path.exists():
        raise FileNotFoundError(f"Reference FASTA not found: {protein_fasta}")
    rec = next(SeqIO.parse(path, "fasta"), None)
    if not rec:
        raise ValueError(f"No sequence in {protein_fasta}")
    seq = str(rec.seq).strip()
    seq = seq[:max_aa].replace("*", "")
    dna = _back_translate(seq)
    out = out_path or tempfile.mktemp(suffix=".fasta", prefix="cas13_ref_")
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w") as f:
        f.write(f">{rec.id}_nt\n{dna}\n")
    return out


def run_magic_blast(
    sra_accessions: List[str],
    query_fasta: str,
    out_tsv: str,
    magicblast_cmd: str = "magicblast",
    num_threads: int = 4,
    evalue: float = 1e-5,
) -> bool:
    """
    Run magic-blast with -sra or -sra_batch against query FASTA. Writes tabular output
    with qseq and sseq so we get subject nucleotide sequence. Returns True if run succeeded.
    """
    Path(out_tsv).parent.mkdir(parents=True, exist_ok=True)
    batch_file = None
    try:
        if len(sra_accessions) == 1:
            sra_arg = ["-sra", sra_accessions[0]]
        else:
            batch_file = out_tsv + ".sra_batch.txt"
            with open(batch_file, "w") as f:
                f.write("\n".join(sra_accessions))
            sra_arg = ["-sra_batch", batch_file]
        cmd = [
            magicblast_cmd,
            "-query", str(query_fasta),
            "-out", out_tsv,
            "-outfmt", "6 std qseq sseq",
            "-num_threads", str(num_threads),
            "-evalue", str(evalue),
        ] + sra_arg
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
        if r.returncode != 0 and r.stderr:
            print(f"[!] magic-blast stderr: {r.stderr[:500]}")
        return r.returncode == 0
    finally:
        if batch_file and Path(batch_file).exists():
            try:
                Path(batch_file).unlink()
            except OSError:
                pass


def _translate_hit(sseq: str, sstart: int, send: int) -> str:
    """Translate subject nucleotide hit to protein. sstart/send 1-based; reverse if sstart > send."""
    s = sseq.replace("-", "").upper()
    if not s:
        return ""
    if sstart > send:
        s = str(Seq(s).reverse_complement())
    # Frame from start (1-based -> 0-based frame)
    frame = (abs(sstart) - 1) % 3
    s = s[frame:]
    return str(Seq(s).translate(to_stop=False))


def parse_magicblast_tsv(tsv_path: str) -> List[Tuple[str, str, str, str]]:
    """
    Parse magic-blast tabular output (format 6 std qseq sseq). Return list of
    (qseqid, sseqid, protein_sequence, nucleotide_sseq) for filtering.
    """
    # Columns: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq
    rows = []
    with open(tsv_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 14:
                continue
            qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore, qseq, sseq = parts[:14]
            try:
                sstart_i, send_i = int(sstart), int(send)
            except ValueError:
                continue
            protein = _translate_hit(sseq, sstart_i, send_i)
            if not protein or "*" in protein:
                protein = protein.split("*")[0]
            rows.append((qseqid, sseqid, protein, sseq))
    return rows


# CRISPR repeat: min length and min tandem count (required for crRNA binding)
CRISPR_MIN_LEN = 20
CRISPR_MIN_COUNT = 2


def _extract_crispr_repeat_sequences(nucleotide: str, min_len: int = 20, min_count: int = 2) -> List[str]:
    """
    Extract CRISPR repeat sequences from nucleotide hit: tandem repeats of length
    min_len repeated at least min_count times. Returns list of unique repeat sequences
    (for metadata; needed for crRNA binding / functional protein).
    """
    n = nucleotide.replace("-", "").upper()
    found = []
    seen_motifs = set()
    for L in range(min_len, min(55, len(n) // min_count + 1)):
        for i in range(len(n) - L * min_count + 1):
            motif = n[i : i + L]
            if motif in seen_motifs:
                continue
            count = 0
            j = i
            while j <= len(n) - L and n[j : j + L] == motif:
                count += 1
                j += L
            if count >= min_count:
                seen_motifs.add(motif)
                found.append(motif)
    return found


def _has_crispr_repeat(nucleotide: str, min_len: int = 20, min_count: int = 2) -> bool:
    """True if nucleotide has a CRISPR-like tandem repeat (required for functional crRNA binding)."""
    return len(_extract_crispr_repeat_sequences(nucleotide, min_len, min_count)) > 0


def passes_cas13_filters(
    protein: str,
    min_tail: int = 15,
    require_m: bool = True,
    nucleotide_hit: Optional[str] = None,
) -> bool:
    """
    Apply filters: length 700-1400 aa, exactly 2 HEPN, N-term M, C-term tail.
    CRISPR repeat is mandatory (nucleotide_hit must contain a tandem repeat) for crRNA binding.
    """
    if not protein or len(protein) < MIN_AA or len(protein) > MAX_AA:
        return False
    motifs = list(HEPN_REGEX.finditer(protein))
    if len(motifs) != HEPN_COUNT_EXACT:
        return False
    if not passes_n_term(protein, require_m):
        return False
    if not passes_c_term(protein, min_tail):
        return False
    if not nucleotide_hit or not _has_crispr_repeat(nucleotide_hit, CRISPR_MIN_LEN, CRISPR_MIN_COUNT):
        return False
    return True


# Delimiter for multiple CRISPR repeat sequences in metadata CSV (no commas)
CRISPR_REPEAT_DELIM = "|"


def mine_sra_with_magicblast(
    sra_runs: List[str],
    reference_fasta: str,
    output_dir: str = "data/raw_sequences",
    magicblast_cmd: str = "magicblast",
    run_batch_size: int = 50,
    num_threads: int = 4,
) -> Tuple[List[Tuple[str, str]], List[Tuple[str, str, str, str]]]:
    """
    For each batch of SRA runs, run magic-blast, parse hits, filter to full Cas13-like ORFs
    with mandatory CRISPR repeat. CRISPR domain sequences are saved in metadata with SRA accession.
    Returns (discoveries, metadata): discoveries = [(seq_id, seq)],
    metadata = [(seq_id, sra_accession, crispr_domain_sequences, score)].
    """
    cfg = get_full_orf_config()
    min_tail = cfg["min_tail"]
    require_m = cfg["require_m"]
    ref_nt_path = build_nucleotide_reference(reference_fasta, max_aa=600)
    discoveries = []
    metadata = []
    seen_seqs = set()

    for i in range(0, len(sra_runs), run_batch_size):
        batch = sra_runs[i : i + run_batch_size]
        tsv = str(Path(output_dir) / f"_magicblast_batch_{i}.tsv")
        ok = run_magic_blast(batch, ref_nt_path, tsv, magicblast_cmd=magicblast_cmd, num_threads=num_threads)
        if not ok or not Path(tsv).exists():
            continue
        for qseqid, sseqid, protein, sseq_nt in parse_magicblast_tsv(tsv):
            if not passes_cas13_filters(
                protein,
                min_tail=min_tail,
                require_m=require_m,
                nucleotide_hit=sseq_nt,
            ):
                continue
            crispr_seqs = _extract_crispr_repeat_sequences(sseq_nt, CRISPR_MIN_LEN, CRISPR_MIN_COUNT)
            crispr_str = CRISPR_REPEAT_DELIM.join(crispr_seqs) if crispr_seqs else ""
            # Dedupe by sequence
            key = protein[:200] + "|" + protein[-200:]
            if key in seen_seqs:
                continue
            seen_seqs.add(key)
            run_id = sseqid.split(".")[0] if "." in sseqid else sseqid.split("_")[0]
            seq_id = f"Cas13_{run_id}_{len(discoveries)}"
            discoveries.append((seq_id, protein))
            metadata.append((seq_id, run_id, crispr_str, "0"))
        try:
            Path(tsv).unlink()
        except OSError:
            pass

    if ref_nt_path.startswith(tempfile.gettempdir()):
        try:
            Path(ref_nt_path).unlink(missing_ok=True)
        except OSError:
            pass
    return discoveries, metadata


def save_discoveries(
    discoveries: List[Tuple[str, str]],
    metadata: List[Tuple[str, str, str, str]],
    output_dir: str = "data/raw_sequences",
) -> Tuple[str, str]:
    """Write deep_hits_YYYYMMDD_HHMMSS.fasta and matching _metadata.csv. Returns (fasta_path, csv_path)."""
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    fasta_path = str(Path(output_dir) / f"deep_hits_{ts}.fasta")
    csv_path = str(Path(output_dir) / f"deep_hits_{ts}_metadata.csv")
    with open(fasta_path, "w") as f:
        for seq_id, seq in discoveries:
            f.write(f">{seq_id}\n{seq}\n")
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        import csv as csv_module
        w = csv_module.writer(f)
        # repeat_domains = CRISPR domain sequences (pipe-separated if multiple) for crRNA binding
        w.writerow(["sequence_id", "sra_accession", "repeat_domains", "score"])
        for row in metadata:
            w.writerow(row)
    return fasta_path, csv_path
