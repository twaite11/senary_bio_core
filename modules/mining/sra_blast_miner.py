"""
SRA Cas13 Discovery Pipeline: Targeted Assembly Strategy
--------------------------------------------------------
Workflow:
1.  SPOT CHECK: Stream first N reads to Diamond. If hits < threshold, SKIP download.
2.  FETCH: Download full SRA reads (fasterq-dump) ONLY if Spot Check passes.
3.  BAIT:  Diamond blastx (Protein search) on full data.
4.  HOOK:  Extract hit reads + their mates (paired-end rescue) by base ID so both ends are kept.
5.  BUILD: Assemble ONLY the extracted reads using Megahit.
6.  HUNT:  Predict genes (Prodigal) and filter for valid Cas13 candidates.

Requirements:
    - Conda packages: sra-tools, diamond, megahit, prodigal, biopython
"""
from __future__ import annotations

import os
import re
import shutil
import subprocess
import tempfile
import argparse
from pathlib import Path
from datetime import datetime
from typing import List, Tuple, Optional, Dict, Set
from concurrent.futures import ProcessPoolExecutor, as_completed

from Bio import Entrez, SeqIO

# --- Configuration & Defaults ---
Entrez.email = os.environ.get("ENTREZ_EMAIL", "researcher@example.com")
HEPN_REGEX = re.compile(r"R.{4,6}H")

# Default Search: Metagenomes (Dark Matter) + WGS
DEFAULT_SRA_TERM = (
    '(txid2[ORGN] OR txid2157[ORGN] OR txid408169[ORGN] OR txid48479[ORGN]) '
    'AND (wgs[Strategy] OR metagenomic[Strategy]) '
    'AND metagenomic[LibrarySource]'
)

# Filters
MIN_AA = 700
MAX_AA = 1500
HEPN_COUNT_EXACT = 2
CRISPR_MIN_LEN = 20
CRISPR_MIN_COUNT = 2
CRISPR_FLANK_BP = 8000

# Optimization
SPOT_CHECK_READS = 3_000_000  # Number of reads to stream for initial check

# --- Helper Functions ---
def passes_n_term(protein: str, require_m: bool = True) -> bool:
    if not protein: return False
    return protein.startswith("M") if require_m else True

def passes_c_term(protein: str, min_tail: int = 15) -> bool:
    matches = list(HEPN_REGEX.finditer(protein))
    if not matches: return False
    last_hepn_end = matches[-1].end()
    return (len(protein) - last_hepn_end) >= min_tail

# --- SRA Interaction ---

def _parse_runs_from_xml(xml_handle) -> List[str]:
    import xml.etree.ElementTree as ET
    runs = []
    try:
        tree = ET.parse(xml_handle)
        for elem in tree.getroot().iter():
            if elem.tag == "RUN" or elem.tag.endswith("RUN"):
                acc = elem.get("accession") or elem.get("acc")
                if acc and str(acc).startswith("SRR"): runs.append(str(acc))
            if elem.tag.endswith("Run") and elem.text and elem.text.strip().startswith("SRR"):
                runs.append(elem.text.strip())
    except Exception: pass
    return runs

def fetch_sra_run_accessions(term: str, max_records: int = 1000, page_size: int = 500) -> List[str]:
    """Fetch SRA run accessions with pagination until we have max_records or no more results."""
    print(f"[*] Querying SRA for: {term}")
    run_ids: List[str] = []
    retstart = 0
    try:
        while len(run_ids) < max_records:
            fetch_max = min(page_size, max_records - len(run_ids))
            h = Entrez.esearch(db="sra", term=term, retmax=fetch_max, retstart=retstart)
            rec = Entrez.read(h)
            h.close()
            id_list = rec.get("IdList", [])
            if not id_list:
                break
            print(f"[*] Fetching details for {len(id_list)} experiments (retstart={retstart})...")
            h = Entrez.efetch(db="sra", id=",".join(id_list), retmode="xml")
            batch = _parse_runs_from_xml(h)
            h.close()
            for r in batch:
                if r not in run_ids:
                    run_ids.append(r)
            if len(run_ids) >= max_records:
                break
            retstart += len(id_list)
            if len(id_list) < fetch_max:
                break
    except Exception as e:
        print(f"[!] SRA Error: {e}")
    return list(dict.fromkeys(run_ids))[:max_records]

# --- Pipeline Components ---

def spot_check_sra(run_id: str, db_path: str) -> bool:
    """Stream first N reads to Diamond to check for ANY Cas13 signal before downloading."""
    # fastq-dump -X N -Z SRR... | diamond blastx ...
    # Note: fastq-dump is single-threaded but sufficient for streaming small chunks
    print(f"   [Spot Check] Peeking at first {SPOT_CHECK_READS} reads of {run_id}...")
    
    dump_cmd = ["fastq-dump", "-X", str(SPOT_CHECK_READS), "-Z", "--fasta", "0", run_id]
    diamond_cmd = [
        "diamond", "blastx", "-d", db_path, "-q", "-", "--quiet",
        "-f", "6", "qseqid", "--evalue", "1e-5", "--max-target-seqs", "1"
    ]
    
    try:
        # Pipe fastq-dump output directly to diamond stdin
        p1 = subprocess.Popen(dump_cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        p2 = subprocess.Popen(diamond_cmd, stdin=p1.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if p1.stdout: p1.stdout.close()
        
        output, _ = p2.communicate()
        
        # If output has content, we found a hit
        if output and len(output.strip()) > 0:
            print(f"   [Spot Check] HIT FOUND in {run_id}! Proceeding to full download.")
            return True
        else:
            return False
    except Exception:
        # If streaming fails, default to downloading (safer)
        return True

def dump_sra_run(run_id: str, out_dir: Path) -> Tuple[List[str], bool]:
    """Download SRA to FASTQ."""
    print(f"   [Download] Fetching full {run_id}...")
    cmd = ["fasterq-dump", run_id, "-O", str(out_dir), "-e", "8", "--split-files"]
    try:
        subprocess.run(cmd, check=True, timeout=7200, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        files = list(out_dir.glob(f"{run_id}*.fastq"))
        if len(files) == 2: return [str(f) for f in files], True
        if len(files) == 1: return [str(files[0])], False
        return [], False
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired):
        return [], False

def build_diamond_db(fasta_path: str, db_prefix: str) -> str:
    db_file = f"{db_prefix}.dmnd"
    if not os.path.exists(db_file):
        print(f"[*] Building Diamond DB from {fasta_path}...")
        subprocess.run(["diamond", "makedb", "--in", fasta_path, "-d", db_prefix, "--quiet"], check=True)
    return db_file

def run_diamond(fastq_paths: List[str], db_path: str, out_tsv: str, threads: int = 8) -> bool:
    # Diamond accepts one -q; concatenate R1+R2 so we search all reads
    query_path: str
    if len(fastq_paths) == 1:
        query_path = fastq_paths[0]
    else:
        out_dir = str(Path(out_tsv).parent)
        query_path = os.path.join(out_dir, "_diamond_query_combined.fastq")
        with open(query_path, "w") as out:
            for p in fastq_paths:
                with open(p) as f:
                    out.write(f.read())
    cmd = [
        "diamond", "blastx", "-d", db_path, "-q", query_path, "-o", out_tsv,
        "-p", str(threads), "-f", "6", "qseqid", "--evalue", "1e-5", "--sensitive", "--quiet"
    ]
    try:
        subprocess.run(cmd, check=True)
        ok = os.path.exists(out_tsv) and os.path.getsize(out_tsv) > 0
        if len(fastq_paths) > 1 and os.path.exists(query_path):
            try:
                os.unlink(query_path)
            except OSError:
                pass
        return ok
    except subprocess.CalledProcessError:
        if len(fastq_paths) > 1:
            try:
                os.unlink(query_path)
            except (OSError, NameError):
                pass
        return False

def _read_id_to_base(read_id: str) -> str:
    """Normalize read ID to base (strip .1/.2 or /1/2) for mate pairing."""
    base = read_id.split()[0].strip()
    for suffix in (".1", ".2", "/1", "/2"):
        if base.endswith(suffix):
            return base[: -len(suffix)]
    return base


def extract_reads(fastq_paths: List[str], hit_ids_file: str, out_r1: str, out_r2: Optional[str]):
    """Extract hit reads + their mates (paired-end rescue) by base ID; keeps both R1 and R2 for any hit."""
    # Build set of base IDs that had a hit (so we keep BOTH mates for any hit)
    bases_with_hits: Set[str] = set()
    with open(hit_ids_file) as f:
        for line in f:
            clean_id = line.strip().split()[0]
            if clean_id:
                bases_with_hits.add(_read_id_to_base(clean_id))

    # Write expanded ID list for seqtk (both mates per hit)
    if not bases_with_hits:
        return
    # seqtk subseq expects one ID per line; we need to pass concrete read IDs.
    # FASTQ headers vary (SRR123.1 vs SRR123/1), so we can't know exact IDs without scanning.
    # So we always use Python fallback for correct mate expansion: keep read if its base is in bases_with_hits.
    def filter_fq(in_fq: str, out_fq: str) -> None:
        Path(out_fq).parent.mkdir(parents=True, exist_ok=True)
        with open(in_fq) as fin, open(out_fq, "w") as fout:
            while True:
                head = fin.readline()
                if not head:
                    break
                seq = fin.readline()
                plus = fin.readline()
                qual = fin.readline()
                curr_id = head[1:].split()[0].split("/")[0].strip()
                base = _read_id_to_base(curr_id)
                if base in bases_with_hits:
                    fout.write(head + seq + plus + qual)

    filter_fq(fastq_paths[0], out_r1)
    if len(fastq_paths) > 1 and out_r2:
        filter_fq(fastq_paths[1], out_r2)

def run_megahit(r1: str, r2: Optional[str], out_dir: str, threads: int) -> Optional[str]:
    cmd = ["megahit", "-o", out_dir, "-t", str(threads), "--min-contig-len", "1000"]
    if r2: cmd.extend(["-1", r1, "-2", r2])
    else: cmd.extend(["-r", r1])
    
    if os.path.exists(out_dir): shutil.rmtree(out_dir)
    try:
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        contigs = Path(out_dir) / "final.contigs.fa"
        return str(contigs) if contigs.exists() else None
    except subprocess.CalledProcessError:
        return None

def run_prodigal_and_filter(contigs_path: str, out_base: Path) -> List[Tuple]:
    faa = out_base.with_suffix(".faa")
    gff = out_base.with_suffix(".gff")
    try:
        subprocess.run(["prodigal", "-i", contigs_path, "-a", str(faa), "-o", str(gff), "-p", "meta", "-q"], check=True)
    except subprocess.CalledProcessError:
        return []

    candidates = []
    contigs = {rec.id: str(rec.seq).upper() for rec in SeqIO.parse(contigs_path, "fasta")}

    for rec in SeqIO.parse(faa, "fasta"):
        prot_seq = str(rec.seq)
        if not (MIN_AA <= len(prot_seq) <= MAX_AA):
            print(f"   [Filter] {rec.id}: filtered out - length {len(prot_seq)} aa not in {MIN_AA}-{MAX_AA}")
            continue
        hepn_count = len(list(HEPN_REGEX.finditer(prot_seq)))
        if hepn_count != HEPN_COUNT_EXACT:
            print(f"   [Filter] {rec.id}: filtered out - HEPN count {hepn_count} (required {HEPN_COUNT_EXACT})")
            continue
        if not passes_n_term(prot_seq):
            print(f"   [Filter] {rec.id}: filtered out - missing N-term Met")
            continue
        if not passes_c_term(prot_seq):
            print(f"   [Filter] {rec.id}: filtered out - insufficient C-term tail after last HEPN")
            continue
        contig_id = "_".join(rec.id.rsplit("_", 1)[:-1])
        if contig_id not in contigs:
            print(f"   [Filter] {rec.id}: filtered out - contig {contig_id} not found")
            continue
        repeats = _find_repeats(contigs[contig_id])
        if not repeats:
            print(f"   [Filter] {rec.id}: filtered out - no CRISPR repeat on contig")
            continue
        candidates.append((rec.id, prot_seq, repeats))
    return candidates

def _find_repeats(dna: str) -> List[str]:
    """Find CRISPR-like tandem (non-overlapping) repeats: same motif repeated >= min_count times in a row."""
    found: List[str] = []
    seen_motifs: Set[str] = set()
    min_count = max(2, CRISPR_MIN_COUNT)
    for length in range(CRISPR_MIN_LEN, min(48, len(dna) // min_count + 1)):
        i = 0
        while i <= len(dna) - length * min_count:
            motif = dna[i : i + length]
            if motif in seen_motifs:
                i += 1
                continue
            count = 0
            j = i
            while j <= len(dna) - length and dna[j : j + length] == motif:
                count += 1
                j += length
            if count >= min_count:
                seen_motifs.add(motif)
                found.append(motif)
            i += 1
    return found

# --- Worker Function ---

def process_single_run(sra_id: str, db_path: str, threads: int) -> List[Tuple]:
    """Worker process for a single SRA run."""
    results = []
    temp_dir = Path(tempfile.gettempdir()) / f"cas13_{sra_id}"
    temp_dir.mkdir(parents=True, exist_ok=True)
    
    try:
        # 1. Spot Check
        if not spot_check_sra(sra_id, db_path):
            return []

        # 2. Dump
        fastq_files, is_paired = dump_sra_run(sra_id, temp_dir)
        if not fastq_files: return []

        # 3. Bait
        hits_tsv = str(temp_dir / "hits.tsv")
        if not run_diamond(fastq_files, db_path, hits_tsv, threads=threads):
            return []
        
        # 4. Hook
        r1_sub = str(temp_dir / "r1_sub.fq")
        r2_sub = str(temp_dir / "r2_sub.fq") if is_paired else None
        extract_reads(fastq_files, hits_tsv, r1_sub, r2_sub)
        
        # 5. Build
        asm_dir = temp_dir / "megahit_out"
        contigs_file = run_megahit(r1_sub, r2_sub, str(asm_dir), threads=threads)
        if not contigs_file: return []
        
        # 6. Hunt
        candidates = run_prodigal_and_filter(contigs_file, temp_dir / "genes")
        for cid, seq, crispr in candidates:
            results.append((sra_id, cid, seq, crispr))
            
    except Exception as e:
        print(f"   [Error] {sra_id}: {e}")
    finally:
        # Always remove temp (SRR dumps, Diamond output, Megahit dir) after bait/hook/build/hunt
        shutil.rmtree(temp_dir, ignore_errors=True)

    return results

# --- Main Driver ---

def mine_sra(sra_list: List[str], ref_fasta: str, out_dir: str, workers: int, threads_per_worker: int):
    root = Path(out_dir)
    root.mkdir(parents=True, exist_ok=True)
    db_path = build_diamond_db(ref_fasta, str(root / "cas13_ref"))
    
    import csv as csv_module
    total_hits = 0
    output_created = False
    out_fa_path = None
    out_csv_path = None

    print(f"[*] Starting mining with {workers} workers ({threads_per_worker} threads each)...")

    with ProcessPoolExecutor(max_workers=workers) as executor:
        futures = {executor.submit(process_single_run, sid, db_path, threads_per_worker): sid for sid in sra_list}

        for future in as_completed(futures):
            sid = futures[future]
            try:
                res = future.result()
                if res:
                    print(f"[*] {sid}: Found {len(res)} candidates! Saving incrementally...")
                    if not output_created:
                        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
                        out_fa_path = root / f"deep_hits_{timestamp}.fasta"
                        out_csv_path = root / f"deep_hits_{timestamp}_metadata.csv"
                        with open(out_csv_path, "w", newline="", encoding="utf-8") as fc:
                            w = csv_module.writer(fc)
                            w.writerow(["sequence_id", "sra_accession", "repeat_domains", "score"])
                        output_created = True

                    with open(out_fa_path, "a") as fa, open(out_csv_path, "a", newline="", encoding="utf-8") as fc:
                        w = csv_module.writer(fc)
                        for run, cid, seq, repeats in res:
                            seq_id = f"Cas13_{cid}_{total_hits}" if not cid.startswith("Cas13_") else cid
                            fa.write(f">{seq_id}\n{seq}\n")
                            repeat_domains = "|".join(repeats[:3]) if repeats else ""
                            w.writerow([seq_id, run, repeat_domains, "0"])
                            print(f"   [PASS] Candidate {seq_id} passed all requirements and written to {out_fa_path}")
                            total_hits += 1
                else:
                    print(f"[*] {sid}: No hits.")
            except Exception as e:
                print(f"[!] {sid} worker failed: {e}")

    if output_created:
        print(f"\n[DONE] Saved {total_hits} hits to {out_fa_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Cas13 SRA Miner")
    parser.add_argument("--ref", required=True, help="Path to known Cas13 proteins FASTA")
    parser.add_argument("--out", default="./cas13_results", help="Output directory")
    parser.add_argument("--limit", type=int, default=10, help="Max SRA runs")
    parser.add_argument("--workers", type=int, default=4, help="Parallel SRA downloads/process (Default: 4)")
    parser.add_argument("--threads", type=int, default=8, help="Threads per worker (Total CPU = workers*threads)")
    args = parser.parse_args()
    
    runs = fetch_sra_run_accessions(DEFAULT_SRA_TERM, max_records=args.limit)
    if runs:
        mine_sra(runs, args.ref, args.out, args.workers, args.threads)
    else:
        print("No SRA runs found matching query.")