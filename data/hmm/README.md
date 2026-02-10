# Pfam HMMs for Cas13 / Type VI mining

This directory holds **Pfam** profile HMMs used for HMMER-based screening of ORFs before ESM-2 scoring. Mining uses these only when **`USE_HMMER=1`** (off by default).

## Official Pfam accessions (exact denotations)

| Accession | ID   | Description |
|-----------|------|-------------|
| **PF05168** | HEPN | Higher Eukaryotes and Prokaryotes Nucleotide-binding domain; nuclease lobe in Cas13 and other Type VI / defense systems. |

## How to populate this directory

**Option 1: Run the fetch script (recommended)**

From the project root:

```bash
python scripts/fetch_pfam_cas13_hmms.py
```

This downloads the Pfam release, extracts the Cas13/HEPN HMMs listed above, and writes them into `data/hmm/`.

**Option 2: Manual fetch with HMMER**

If you have HMMER installed:

```bash
cd data/hmm
wget https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
gunzip Pfam-A.hmm.gz
hmmfetch Pfam-A.hmm PF05168 > PF05168.hmm
# Remove the full Pfam-A.hmm to save space: rm Pfam-A.hmm
```

## Env vars (mining)

- **`USE_HMMER`** – Set to `1` (or `true`/`yes`) to enable HMMER screening. **Default: off.** When off, all size+full_ORF passing ORFs go straight to ESM-2.
- **`HMM_DIR`** – Directory containing `.hmm` files (default: `data/hmm`).
- **`HMM_EVALUE`** – Max E-value for an ORF to count as a hit (default: `1e-5`).

When HMMER is on, only ORFs that hit at least one of these HMMs (e.g. PF05168) below the E-value threshold are sent to ESM-2 for embedding and scoring.
