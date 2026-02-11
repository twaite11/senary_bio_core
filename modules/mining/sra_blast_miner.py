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

Updates:
- Discovery mode: Relaxed filters (MIN_AA 500, MAX_AA 1800; MIN_HEPN_COUNT 1; accepts partials).
- N-term: 'Starts with M' check disabled by default (require_m=False); accepts partial genes.
- DEBUG: Prints raw vs filtered protein counts when all proteins are filtered out.
- hit_type: Outputs "System" (CRISPR on contig) vs "Orphan"; both accepted.
- seqtk: Uses seqtk for read extraction when available, else Python fallback.

Requirements:
    - Conda packages: sra-tools, diamond, megahit, prodigal, biopython, seqtk (optional)
"""
from __future__ import annotations

import os
import re
import shutil
import subprocess
import tempfile
import argparse
import csv
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

# Filters (relaxed for discovery: catch partials and smaller fragments)
MIN_AA = 500
MAX_AA = 1800
MIN_HEPN_COUNT = 1  # Keep at 1 to catch partials
CRISPR_MIN_LEN = 20
CRISPR_MIN_COUNT = 2
CRISPR_FLANK_BP = 8000

# Optimization
SPOT_CHECK_READS = 2_000_000  # Number of reads to stream for initial check

# --- Helper Functions ---
def passes_n_term(protein: str, require_m: bool = False) -> bool:
    """Discovery mode: default require_m=False to accept partial genes."""
    if not protein:
        return False
    return protein.startswith("M") if require_m else True

def passes_c_term(protein: str, min_tail: int = 10) -> bool:
    """Require at least one HEPN; for discovery mode tail check is lenient."""
    matches = list(HEPN_REGEX.finditer(protein))
    if not matches:
        return False
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
    seen_ids: Set[str] = set()
    retstart = 0
    try:
        while len(run_ids) < max_records:
            h = Entrez.esearch(db="sra", term=term, retmax=page_size, retstart=retstart)
            rec = Entrez.read(h)
            h.close()
            id_list = rec.get("IdList", [])
            if not id_list:
                break
            print(f"[*] Fetching details for {len(id_list)} experiments (retstart={retstart}, found={len(run_ids)}/{max_records})...")
            h = Entrez.efetch(db="sra", id=",".join(id_list), retmode="xml")
            batch = _parse_runs_from_xml(h)
            h.close()
            for r in batch:
                if r not in seen_ids:
                    seen_ids.add(r)
                    run_ids.append(r)
                    if len(run_ids) >= max_records:
                        break
            retstart += len(id_list)
            if len(id_list) < page_size:
                break
    except Exception as e:
        print(f"[!] SRA Error: {e}")
    return run_ids[:max_records]

# --- Pipeline Components ---

def spot_check_sra(run_id: str, db_path: str) -> bool:
    """Stream first N reads to Diamond to check for ANY Cas13 signal before downloading."""
    # fastq-dump -X N -Z SRR... | diamond blastx ...
    # Note: fastq-dump is single-threaded but sufficient for streaming small chunks
    print(f"   [Spot Check] Peeking at first {SPOT_CHECK_READS} reads of {run_id}...")
    
    dump_cmd = ["fastq-dump", "-X", str(SPOT_CHECK_READS), "-Z", "--fasta", "0", run_id]
    diamond_cmd = [
        "diamond", "blastx", "-d", db_path, "-q", "-", "--quiet",
        "-f", "6", "qseqid", "--evalue", "1e-5", "--max-target-seqs", "1",
        "--sensitive"
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
    """Download SRA to FASTQ. Uses -t/--temp in out_dir for cleanup."""
    print(f"   [Download] Fetching full {run_id}...")
    cmd = ["fasterq-dump", run_id, "-O", str(out_dir), "-t", str(out_dir), "-e", "8", "--split-files", "--force"]
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
    """Extract hit reads + their mates (paired-end rescue). Tries seqtk if available, else Python by base ID."""
    use_seqtk = shutil.which("seqtk") is not None
    if use_seqtk:
        try:
            Path(out_r1).parent.mkdir(parents=True, exist_ok=True)
            with open(out_r1, "w") as f:
                subprocess.run(["seqtk", "subseq", fastq_paths[0], hit_ids_file], stdout=f, check=True)
            if len(fastq_paths) > 1 and out_r2:
                with open(out_r2, "w") as f:
                    subprocess.run(["seqtk", "subseq", fastq_paths[1], hit_ids_file], stdout=f, check=True)
            return
        except subprocess.CalledProcessError:
            pass  # Fallback to Python

    bases_with_hits: Set[str] = set()
    with open(hit_ids_file) as f:
        for line in f:
            clean_id = line.strip().split()[0]
            if clean_id:
                bases_with_hits.add(_read_id_to_base(clean_id))
    if not bases_with_hits:
        return

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
    # --k-min 21 helps with lower coverage assemblies
    cmd = ["megahit", "-o", out_dir, "-t", str(threads), "--min-contig-len", "800", "--k-min", "21"]
    if r2: cmd.extend(["-1", r1, "-2", r2])
    else: cmd.extend(["-r", r1])
    
    if os.path.exists(out_dir): shutil.rmtree(out_dir)
    try:
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        contigs = Path(out_dir) / "final.contigs.fa"
        return str(contigs) if contigs.exists() else None
    except subprocess.CalledProcessError:
        return None

def run_prodigal_and_filter(contigs_path: str, out_base: Path) -> Tuple[List[Tuple], int]:
    """Returns (candidates, raw_protein_count). Candidates: (rec.id, prot_seq, repeats, hit_type)."""
    faa = out_base.with_suffix(".faa")
    gff = out_base.with_suffix(".gff")
    try:
        subprocess.run(["prodigal", "-i", contigs_path, "-a", str(faa), "-o", str(gff), "-p", "meta", "-q"], check=True)
    except subprocess.CalledProcessError:
        return [], 0

    candidates = []
    contigs = {rec.id: str(rec.seq).upper() for rec in SeqIO.parse(contigs_path, "fasta")}
    raw_proteins = list(SeqIO.parse(faa, "fasta"))
    raw_count = len(raw_proteins)

    for rec in raw_proteins:
        prot_seq = str(rec.seq)
        if not (MIN_AA <= len(prot_seq) <= MAX_AA):
            continue
        hepn_count = len(list(HEPN_REGEX.finditer(prot_seq)))
        if hepn_count < MIN_HEPN_COUNT:
            continue
        if not passes_n_term(prot_seq, require_m=False):
            continue
        if not passes_c_term(prot_seq, min_tail=10):
            continue
        repeats: List[str] = []
        hit_type = "Orphan"
        contig_id = "_".join(rec.id.rsplit("_", 1)[:-1])
        if contig_id in contigs:
            repeats = _find_repeats(contigs[contig_id])
            if repeats:
                hit_type = "System"
        candidates.append((rec.id, prot_seq, repeats, hit_type))
    return candidates, raw_count

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
        if not contigs_file:
            print(f"   [Warn] {sra_id}: Assembly failed (0 contigs).")
            return []

        # 6. Hunt
        candidates, raw_count = run_prodigal_and_filter(contigs_file, temp_dir / "genes")
        if raw_count > 0 and len(candidates) == 0:
            print(f"   [Debug] {sra_id}: Assembled {raw_count} proteins, but ALL were filtered out. (Check filters?)")
        for cid, seq, crispr, h_type in candidates:
            results.append((sra_id, cid, seq, crispr, h_type))

    except Exception as e:
        print(f"   [Error] {sra_id}: {e}")
    finally:
        shutil.rmtree(temp_dir, ignore_errors=True)
        try:
            for stray in Path.cwd().glob(f"{sra_id}*"):
                if stray.is_file() and stray.suffix in (".sra", ".fastq"):
                    stray.unlink()
        except Exception:
            pass

    return results

# --- Main Driver ---

def mine_sra(sra_list: List[str], ref_fasta: str, out_dir: str, workers: int, threads_per_worker: int):
    root = Path(out_dir)
    root.mkdir(parents=True, exist_ok=True)
    db_path = build_diamond_db(ref_fasta, str(root / "cas13_ref"))

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_fa_path = root / f"deep_hits_{timestamp}.fasta"
    out_csv_path = root / f"deep_hits_{timestamp}_metadata.csv"
    with open(out_csv_path, "w", newline="", encoding="utf-8") as fc:
        w = csv.writer(fc)
        w.writerow(["sequence_id", "sra_accession", "repeat_domains", "hit_type", "score"])
    open(out_fa_path, "a").close()

    total_hits = 0
    print(f"[*] Starting mining with {workers} workers ({threads_per_worker} threads each)...")
    print(f"[*] Results will be streamed to: {out_fa_path}")

    with ProcessPoolExecutor(max_workers=workers) as executor:
        futures = {executor.submit(process_single_run, sid, db_path, threads_per_worker): sid for sid in sra_list}

        for future in as_completed(futures):
            sid = futures[future]
            try:
                res = future.result()
                if res:
                    print(f"[*] {sid}: Found {len(res)} candidates! Saving immediately...")
                    with open(out_fa_path, "a") as fa, open(out_csv_path, "a", newline="", encoding="utf-8") as fc:
                        w = csv.writer(fc)
                        for idx, (run, cid, seq, repeats, h_type) in enumerate(res):
                            seq_id = f"Cas13_{cid}_{total_hits + idx}" if not cid.startswith("Cas13_") else cid
                            fa.write(f">{seq_id}\n{seq}\n")
                            repeat_domains = "|".join(repeats[:3]) if repeats else ""
                            w.writerow([seq_id, run, repeat_domains, h_type, "0"])
                    total_hits += len(res)
                else:
                    print(f"[*] {sid}: Finished. No candidates passing filters.")
            except Exception as e:
                print(f"[!] {sid} worker failed: {e}")

    print(f"\n[DONE] Processed all runs. Total unique hits saved: {total_hits}")

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