# SRA Cas13 miner on a VPS

## How the two entry points work together

- **`scripts/run_sra_cas13_search.py`** – Recommended CLI. Parses args (e.g. `--ref`, `--out`, `--max-runs`, `--workers`, `--num-threads`), fetches the SRA run list via **`fetch_sra_run_accessions()`**, then calls **`mine_sra()`** in `modules/mining/sra_blast_miner.py`. All pipeline logic lives in the miner module; the script is a thin wrapper.
- **`modules/mining/sra_blast_miner.py`** – The actual pipeline (spot-check → fetch → Diamond → Hook → Megahit → Prodigal → filter). You can run it **directly** with `python -m modules.mining.sra_blast_miner --ref ... --out ...`; it uses its own argparse and the same `fetch_sra_run_accessions` + `mine_sra`. Use the script when you want defaults and env-based config; use the module when you want one-liners from the repo root.

Both produce **`deep_hits_<timestamp>.fasta`** and **`deep_hits_<timestamp>_metadata.csv`** in the output directory.

---

## Dependencies to install on the VPS

### 1. System / Conda tools (on PATH)

- **SRA Toolkit** – for `fastq-dump` (spot-check) and `fasterq-dump` (full download)  
  - Conda: `conda install -c bioconda sra-tools`  
  - Or: [NCBI SRA Toolkit](https://github.com/ncbi/sra-tools)
- **Diamond** – protein search  
  - Conda: `conda install -c bioconda diamond`
- **Megahit** – assembler  
  - Conda: `conda install -c bioconda megahit`
- **Prodigal** – gene prediction  
  - Conda: `conda install -c bioconda prodigal`

One-liner (conda, with a bio env):

```bash
conda create -n cas13 python=3.10 -y
conda activate cas13
conda install -c bioconda sra-tools diamond megahit prodigal -y
```

### 2. Python

- From repo root: `pip install -r requirements.txt` (needs **biopython** for the miner).

---

## Commands to run on the VPS

Clone the repo and install Python deps:

```bash
cd /path/to/collateral_bio_core
pip install -r requirements.txt
```

**Option A – Use the script (recommended)**  
Reference FASTA is required (e.g. Cas13 proteins). Default ref: `data/references/mining_refs.fasta`.

```bash
# Small test (e.g. 5 runs, 4 workers, 8 threads each)
python scripts/run_sra_cas13_search.py \
  --ref data/references/mining_refs.fasta \
  --out data/raw_sequences \
  --max-runs 5 \
  --workers 4 \
  --num-threads 8

# Dry run: only fetch and print SRA run list
python scripts/run_sra_cas13_search.py --ref data/references/mining_refs.fasta --dry-run --max-runs 20

# Larger run (e.g. 50 runs)
python scripts/run_sra_cas13_search.py \
  --ref data/references/mining_refs.fasta \
  --out data/raw_sequences \
  --max-runs 50 \
  --workers 4 \
  --num-threads 8
```

**Option B – Run the miner module directly**

```bash
python -m modules.mining.sra_blast_miner \
  --ref data/references/mining_refs.fasta \
  --out data/raw_sequences \
  --limit 10 \
  --workers 4 \
  --threads 8
```

**Environment variables (optional)**  
- `SRA_TERM` – NCBI SRA query (default: metagenomes + WGS).  
- `REFERENCE_FASTA` / `OUTPUT_DIR` – override ref and output dir.  
- `SRA_WORKERS`, `NUM_THREADS` – workers and threads per worker when using the script.  
- `ENTREZ_EMAIL` – set to your email for NCBI (e.g. `export ENTREZ_EMAIL=you@example.com`).

**Output**  
- `data/raw_sequences/deep_hits_YYYYMMDD_HHMMSS.fasta` – candidate protein sequences.  
- `data/raw_sequences/deep_hits_YYYYMMDD_HHMMSS_metadata.csv` – columns: `sequence_id`, `sra_accession`, `repeat_domains` (CRISPR repeats), `score`.  

Use the FASTA (or a symlink like `deep_hits_latest.fasta`) as input to `run_full_pipeline.py` for design → structure → identity → matchmaker.
