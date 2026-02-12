<div align="center">

# SENARY BIO

### Computational Discovery Platform for Novel Type VI CRISPR Enzyme Therapeutics

<br>

```
  Mine  ──────>  Design  ──────>  Structure  ──────>  Identity  ──────>  Match
  NCBI WGS        ESM-2           OmegaFold          Drift <85%         Enzyme x
  6-frame ORF     Mutate          TM-score            vs Known          Fusion
  CRISPR loci     Drift           Bi-lobed HEPN       Cas13              Targets
```

<br>

![Python](https://img.shields.io/badge/Python-3.8+-3776AB?style=flat-square&logo=python&logoColor=white)
![PyTorch](https://img.shields.io/badge/PyTorch-2.0+-EE4C2C?style=flat-square&logo=pytorch&logoColor=white)
![ESM-2](https://img.shields.io/badge/ESM--2-35M-00897B?style=flat-square)
![OmegaFold](https://img.shields.io/badge/OmegaFold-Structure-6A1B9A?style=flat-square)
![Gemini](https://img.shields.io/badge/Gemini-AI_Agent-4285F4?style=flat-square&logo=google&logoColor=white)
![License](https://img.shields.io/badge/License-Proprietary-333333?style=flat-square)

</div>

---

<br>

## Table of Contents

| # | Section | Description |
|:-:|---------|-------------|
| 1 | [Mission](#mission) | What we're building and why |
| 2 | [Platform Overview](#platform-overview) | End-to-end discovery pipeline |
| 3 | [Structure Pipeline](#structure-pipeline--dashboard) | OmegaFold, TM-score, interactive dashboard |
| 4 | [Enzyme Mining](#enzyme-mining--the-armory) | SRA/WGS prospecting modules |
| 5 | [Target Discovery](#target-discovery--the-vault) | Fusion RNA target identification |
| 6 | [Matchmaker](#matchmaker) | Enzyme-target pairing engine |
| 7 | [Expert Agent](#expert-agent) | AI-driven lead validation |
| 8 | [Quick Start](#quick-start) | Setup, data, and workflow |
| 9 | [Project Structure](#project-structure) | Repository layout |
| 10 | [Configuration](#configuration) | Environment variables and tuning |

<br>

---

## Mission

**Senary Bio** is building precision cancer therapeutics through novel Type VI CRISPR enzyme discovery. The platform identifies and validates high-collateral Type VI CRISPR enzyme variants capable of targeting cancer-specific fusion RNAs -- creating a "suicide switch" mechanism that induces apoptosis specifically in tumor cells while preserving healthy tissue.

We combine computational biology, deep learning, and high-throughput structural screening to accelerate the discovery of next-generation RNA-guided therapeutics.

<br>

---

## Platform Overview

> See **[FILTERING_MAP.md](FILTERING_MAP.md)** for the complete filtering map from NCBI scraping through final candidate selection.

<br>

<table>
<tr><th>Stage</th><th>Process</th><th>Output</th></tr>
<tr>
  <td><b>1. Mine</b></td>
  <td>SRA/WGS in obscure environments; 6-frame translate; full-enzyme ORF validation (N-term M, C-term tail, not truncated); CRISPR repeat metadata</td>
  <td><code>deep_hits_*.fasta</code></td>
</tr>
<tr>
  <td><b>2. Design</b></td>
  <td>ESM-2 embedding; optional interpolation; mutate for drift (stability + &lt;85% identity to known Cas13)</td>
  <td><code>drift_variants.fasta</code></td>
</tr>
<tr>
  <td><b>3. Structure</b></td>
  <td>OmegaFold prediction; bi-lobed / HEPN structural check; TM-score vs 34 Cas13 family references</td>
  <td><code>passed_structures.fasta</code></td>
</tr>
<tr>
  <td><b>4. Identity</b></td>
  <td>Keep only sequences &lt;85% identical to known Cas13 (drift goal)</td>
  <td><code>passed.fasta</code></td>
</tr>
<tr>
  <td><b>5. Match</b></td>
  <td>Enzyme x fusion targets; PFS rule; ARCHS4 expression; Gemini AI verdict</td>
  <td><code>lead_candidates.csv</code></td>
</tr>
</table>

<br>

```
 THE ARMORY (Enzyme Mining)                         THE VAULT (Target Discovery)
 ┌─────────────────────────────────────┐            ┌─────────────────────────────────────┐
 │  sra_scout                          │            │  fusion_metadata                    │
 │  autonomous_prospector              │            │  specificity_filter                 │
 │  full_orf_checks                    │            │  archs4_loader                      │
 │  CRISPR repeat metadata             │            │  mutation_loader                    │
 └──────────────────┬──────────────────┘            └──────────────────┬──────────────────┘
                    │                                                  │
                    └──────────────────────┬───────────────────────────┘
                                           │
                                           v
 ┌─────────────────────────────────────────────────────────────────────────────────────────┐
 │  FULL PIPELINE                                                                         │
 │  embed_pool -> mutate_for_drift -> structure_filter -> identity_filter                  │
 │  run_full_pipeline.py  ->  data/identity_filtered/passed.fasta                         │
 └─────────────────────────────────────────┬───────────────────────────────────────────────┘
                                           │
                                           v
                              ┌────────────────────────────┐
                              │  MATCHMAKER -> EXPERT AGENT │
                              │  lead_candidates.csv        │
                              └────────────────────────────┘
```

<br>

---

## Structure Pipeline & Dashboard

For 2-3 HEPN filtered enzymes, the structure pipeline predicts 3D conformations with **OmegaFold**, computes TM-score homology against all known Cas13-family crystal structures (34 references across Cas13a/b/d/bt/X/Y), and builds an interactive dashboard with domain coloring, motif tables, and a 3D viewer.

Each sequence is predicted in its own subprocess to prevent CUDA memory fragmentation. The pipeline supports resume (skips existing PDBs) and early termination on persistent OOM.

<br>

### Domain Color Scheme

When InterPro annotations are available, the 3D viewer colors structures by Pfam domain:

| Color | Domain | Pfam | Role |
|:-----:|--------|:----:|------|
| **Cyan** | HEPN | PF05168 | Nucleotide-binding nuclease domain; cleaves RNA via RxxxxH motif. Cas13 has 2-3 HEPN domains. |
| **Orange** | HEL | PF01228 | Helicase-like domain; crRNA binding and target recognition. |
| **Slate** | WYL | PF18456 | WYL effector-associated domain; found in some Type VI systems. |
| *Spectrum* | -- | -- | Unmapped regions shown in rainbow gradient by chain position. |

<br>

### Dashboard Columns

| Column | Description |
|--------|-------------|
| **Protein ID / Length** | Identifier and total amino acids |
| **HEPN count** | Number of HEPN domains detected |
| **NUC / REC** | Predicted nuclease and recognition lobe residue spans |
| **HEPN1/2 TM** | Local TM-score for each HEPN domain vs reference |
| **Catalytic distance** | Distance between RxxxxH catalytic histidines (Angstroms) |
| **Best Reference** | Highest-scoring Cas13 family reference (dynamic across all 34) |
| **TM %** | TM-score percentage against best reference |
| **3D Viewer** | Interactive OmegaFold structure (3Dmol.js) with domain coloring |

<br>

### Running the Structure Pipeline

```bash
# 1. Filter to 2-3 HEPN sequences
python visualization/filter_23_hepn.py \
    --input data/raw_sequences/deep_hits_YYYYMMDD.fasta \
    --output data/structure_pipeline/input_2-3_hepn.fasta

# 2. Download all Cas13-family reference PDBs from RCSB (one-time)
python visualization/download_cas13_references.py

# 3. Run OmegaFold (one subprocess per sequence, auto-resume)
python visualization/run_omegafold.py \
    --omegafold-repo /path/to/OmegaFold \
    --max-residues 750          # optional: skip sequences too large for GPU

# 4. Compute TM-scores vs all references
python visualization/run_tmscore.py

# 5. Generate dashboard
python visualization/structure_dashboard.py

# 6. Serve and view
python -m http.server 8000
# Open: http://localhost:8000/visualization/structure_dashboard.html
```

> **OmegaFold on Python 3.11+:** Clone the repo and use a separate Python 3.10 venv.
> The pipeline auto-detects `venv_omegafold1` in the parent directory.
> See [docs/VPS_RUN_PLAN.md](docs/VPS_RUN_PLAN.md) for full GPU setup instructions.

<br>

---

## Enzyme Mining -- The Armory

| Module | Purpose | Details |
|--------|---------|---------|
| `sra_scout` | WGS metagenome search | Normalizes query; tries `wgs[Prop]` with fallback; 6-frame translate; HEPN + topology; full ORF checks |
| `autonomous_prospector` | AI-driven continuous mining | LLM strategy -> SRA search -> semantic filter -> ESM-2 scoring -> deep mine ORFs 600-1400 aa; optional CRISPR/HMMER screening; pagination + shotgun mode |
| `full_orf_checks` | Full-enzyme ORF validation | N-term M start; C-term tail after last HEPN; contig-boundary rejection |
| `deep_miner_utils` | Deep learning engine | ESM-2 35M scoring vs RfxCas13d / PspCas13a references; CRISPR array detection |
| `family_grouper` | Homology clustering | ESM-2 based clustering; SN01_001 naming scheme; family FASTA output |
| `hmmer_miner` | HMM pre-screening (optional) | Pfam Cas13/HEPN HMMs (PF05168) before ESM-2; requires `USE_HMMER=1` |

<br>

---

## Target Discovery -- The Vault

| Module | Purpose | Details |
|--------|---------|---------|
| `fusion_metadata` | Fusion-cancer mapping | Loads recurrence tables + novel matrix; builds fusion-to-TCGA mapping; organ keyword resolution |
| `specificity_filter` | High-specificity targets | Keeps fusions in <=3 tissue types; outputs `high_specificity_targets.csv` |
| `mutation_loader` | VCF mutation mining | Parses VCF for gene-specific mutations (e.g. KRAS G12C) |
| `archs4_loader` | Expression and safety | HDF5 human_matrix; organ-specific enrichment; fusion absent-in-normal / present-in-cancer |

<br>

---

## Matchmaker

Loads enzymes (FASTA) and targets (`high_specificity_targets.csv` or `known_fusions.csv`). Screens every enzyme-target pair against disease maps, applies PFS rules (no G at 3'), and outputs ranked `lead_candidates.csv`.

<br>

---

## Expert Agent

Loads matchmaker output, groups by (Target_Fusion, Associated_Disease) to minimize API calls, runs ARCHS4 organ-specific enrichment, and queries Gemini AI for a structured GO / NO-GO / HOLD verdict with screening strategy. Outputs `lead_candidates_filtered.csv`.

<br>

---

## Quick Start

### Prerequisites

```bash
python -m venv venv
source venv/bin/activate        # Linux/macOS
# venv\Scripts\activate         # Windows

pip install -r requirements.txt
```

<details>
<summary><b>Required data files</b></summary>

<br>

| File | Description |
|------|-------------|
| `data/targets/known_fusions.csv` | Validation fusion targets |
| `data/targets/novel_fusions.csv` | Discovery fusion targets |
| `data/matrices/disease_matrix_*.csv` | Fusion x cancer recurrence matrix |
| `data/matrices/KB_and_Pub_Recur_per_cancer.csv` | Knowledge-base recurrence |
| `data/expression_data/human_matrix.h5` | ARCHS4 expression ([download](https://maayanlab.cloud/archs4/)) |
| `data/references/known_cas13.fasta` | Known Cas13 sequences for identity/drift filter |
| `data/references/mining_refs.fasta` | RfxCas13d + PspCas13a for ESM-2 similarity scoring |

Mining outputs land in `data/raw_sequences/deep_hits_*.fasta` and `*_metadata.csv`.

</details>

<br>

### Workflow

```
Step 1   Mine enzymes
         python modules/mining/autonomous_prospector.py

Step 2   Family grouping
         python modules/mining/family_grouper.py

Step 3   Specificity filter
         python modules/targeting/specificity_filter.py

Step 4   Matchmaker
         python modules/matchmaker.py

Step 5   Expert agent
         python modules/analysis/expert_agent.py

Step 6   Structure pipeline
         See "Structure Pipeline & Dashboard" above
```

<br>

### Full Pipeline (single command)

After mining, run the integrated pipeline on any FASTA pool:

```bash
# With GPU (OmegaFold on this machine)
python run_full_pipeline.py \
    --input data/raw_sequences/deep_hits_latest.fasta \
    --run-matchmaker

# Without GPU (skip structure filter)
python run_full_pipeline.py --skip-structure --run-matchmaker

# Auto-detect latest deep_hits
python run_full_pipeline.py --run-matchmaker
```

<details>
<summary><b>Pipeline steps breakdown</b></summary>

<br>

| Step | What it does | Output |
|------|-------------|--------|
| **Embed** | ESM-2 embeddings for all sequences | `data/design/embeddings/` |
| **Mutate** | Score stability vs RfxCas13d/PspCas13a; optional Gemini trans-cleavage prompt; keep variants <85% identity | `data/design/drift_variants.fasta` |
| **Structure** | OmegaFold -> TM-score + HEPN check | `data/structure_pipeline/passed_structures.fasta` |
| **Identity** | Keep only <85% identity to known Cas13 | `data/identity_filtered/passed.fasta` |
| **Match** | Enzyme x fusion targets | `lead_candidates.csv` |

</details>

<details>
<summary><b>Exhaustive mining mode</b></summary>

<br>

Single fixed query (bacteria/archaea WGS + metagenomic), paginate through all NCBI Nucleotide contigs, mine every one without LLM/semantic filter:

```bash
EXHAUSTIVE_MODE=1 FILTER_23_HEPN=1 python modules/mining/autonomous_prospector.py
```

Output flows into the same `deep_hits_*.fasta` path, then filter_23_hepn -> structure pipeline -> dashboard.

</details>

<details>
<summary><b>SRA Diamond pipeline (assembly-based)</b></summary>

<br>

4-step workflow: Diamond blastx (protein-space bait) -> Hook (extract hit + mate reads) -> Megahit (assemble) -> Prodigal + filter (700-1400 aa, 2 HEPN, intact N/C, CRISPR).

```bash
python scripts/run_sra_cas13_search.py
```

Requires on PATH: `fasterq-dump`, `diamond`, `megahit`, `prodigal`.

</details>

<br>

---

## Project Structure

```
collateral_bio_core/
│
├── run_full_pipeline.py              # Orchestration: embed -> mutate -> structure -> identity -> match
├── requirements.txt
├── .env                              # Runtime configuration (from pipeline.env.example)
│
├── modules/
│   ├── mining/
│   │   ├── autonomous_prospector.py  # AI-driven continuous NCBI mining
│   │   ├── sra_scout.py             # SRA/WGS search and ORF extraction
│   │   ├── deep_miner_utils.py      # ESM-2 scoring engine + CRISPR detection
│   │   ├── full_orf_checks.py       # Full-enzyme ORF validation
│   │   ├── family_grouper.py        # ESM-2 homology clustering
│   │   └── hmmer_miner.py           # Optional HMM pre-screening
│   │
│   ├── design/
│   │   ├── embed_pool.py            # ESM-2 embedding
│   │   ├── interpolate.py           # Latent space interpolation
│   │   └── mutate_for_drift.py      # Stability-aware drift mutations
│   │
│   ├── structure_filter/
│   │   ├── predict_structures.py    # OmegaFold runner (1 subprocess per sequence)
│   │   ├── bi_lobed_hepn_check.py   # TM-score + HEPN structural validation
│   │   ├── functional_criteria.py   # Catalytic distance, linker, pLDDT checks
│   │   └── run_filter.py            # Orchestrates predict -> filter -> output
│   │
│   ├── targeting/
│   │   ├── archs4_loader.py         # ARCHS4 expression data
│   │   ├── fusion_metadata.py       # Fusion-cancer mapping
│   │   └── specificity_filter.py    # Tissue-type specificity
│   │
│   ├── analysis/
│   │   └── expert_agent.py          # Gemini AI verdict engine
│   │
│   ├── matchmaker.py                # Enzyme x target pairing
│   └── identity_filter.py           # Drift goal: max identity < 85%
│
├── visualization/
│   ├── download_cas13_references.py  # Fetch all Cas13 PDBs from RCSB
│   ├── filter_23_hepn.py            # 2-3 HEPN motif filter
│   ├── run_omegafold.py             # OmegaFold runner (standalone)
│   ├── run_tmscore.py               # TM-score computation (US-align / tmtools)
│   ├── structure_dashboard.py       # Interactive HTML dashboard generator
│   └── family_dashboard.py          # Family-level visualization
│
├── config/
│   └── pipeline.env.example          # Template for .env configuration
│
├── data/
│   ├── raw_sequences/                # deep_hits_*.fasta + metadata
│   ├── mined_sequences/              # family_grouped_*.fasta
│   ├── design/                       # embeddings/, drift_variants.fasta
│   ├── structure_pipeline/           # structures/omegafold/, references/, passed_structures.fasta
│   ├── identity_filtered/            # passed.fasta, identity_metadata.csv
│   ├── expression_data/              # human_matrix.h5
│   ├── targets/                      # known_fusions.csv, high_specificity_targets.csv
│   ├── matrices/                     # disease_matrix_*.csv
│   └── references/                   # known_cas13.fasta, mining_refs.fasta
│
├── docs/
│   └── VPS_RUN_PLAN.md              # Full VPS/RunPod deployment guide
│
├── scripts/                          # Utility scripts (SRA Diamond, HMM fetch, etc.)
└── utils/                            # split_excel, inspect_archs4_metadata
```

<br>

---

## Configuration

Copy `config/pipeline.env.example` to `.env` and configure:

<br>

<details>
<summary><b>Mining configuration</b></summary>

| Variable | Purpose | Default |
|----------|---------|---------|
| `ESM_THRESHOLD` | Minimum ESM-2 score to keep a hit | `0.75` |
| `ESM_SIMILARITY_CEILING` | Diversity mode: cap similarity | `0.85` |
| `ESM_REFERENCE_FASTA` | RfxCas13d + PspCas13a FASTA for closest-similarity mining | -- |
| `REQUIRE_CRISPR` | Require CRISPR array in contig | `1` |
| `CRISPR_FLANK_BP` | Max distance (bp) between ORF and CRISPR array | `10000` |
| `REQUIRE_START_M` | ORF must start with Met | `1` |
| `MIN_CTERM_TAIL` | Min residues after last HEPN motif | `15` |
| `CONTIG_BOUNDARY_MARGIN` | Reject ORFs within this many nt of contig end | `30` |
| `REQUIRE_FULL_STRUCTURE` | With CRISPR, also require repeat domains | `0` |
| `MIN_REPEAT_COUNT` | Min CRISPR repeat sequences per hit | `1` |
| `USE_HMMER` | Enable HMM pre-screening (requires `data/hmm/`) | `0` |
| `EXHAUSTIVE_MODE` | Single-query paginated mining (no LLM) | `0` |
| `PROSPECTOR_WORKERS` | Parallel mining workers | `1` |

</details>

<details>
<summary><b>Structure and design configuration</b></summary>

| Variable | Purpose | Default |
|----------|---------|---------|
| `OMEGAFOLD_REPO` | Path to cloned OmegaFold repository | -- |
| `FAMILY_DEVICE` | PyTorch device (`cuda` or `cpu`) | `cuda` |
| `EMBED_BATCH_SIZE` | Sequences per ESM-2 batch | `50` |
| `GEMINI_API_KEY` | For trans-cleavage mutation suggestions | -- |
| `USE_TRANS_CLEAVAGE_PROMPT` | Enable Gemini mutation suggestions | `0` |

</details>

<details>
<summary><b>LLM and API configuration</b></summary>

| Variable | Purpose | Default |
|----------|---------|---------|
| `LLM_PROVIDER` | Prospector LLM backend | `local` |
| `LLM_LOCAL_URL` | Ollama API endpoint | `http://localhost:11434/api/chat` |
| `LLM_MODEL` | Model name for strategy generation | `llama3-bio` |
| `GEMINI_API_KEY` | Google Gemini for expert agent | -- |

</details>

<br>

### Troubleshooting

| Issue | Solution |
|-------|----------|
| "Enzyme file not found" | Matchmaker falls back to mock enzymes. Run mining first. |
| "ARCHS4 file not found" | Download `human_matrix.h5` into `data/expression_data/` |
| Prospector import error | `pip install torch transformers requests` |
| OmegaFold on Python 3.12 | Clone repo + use `--omegafold-repo` or `OMEGAFOLD_REPO` |
| No hits from mining | Relax filters: `REQUIRE_CRISPR=0` or lower `ESM_THRESHOLD` |
| Identity filter keeps all | Add sequences to `data/references/known_cas13.fasta` |
| OmegaFold CUDA OOM | Use `--max-residues 750` or upgrade GPU; pipeline auto-resumes |
| TM-scores all null | Install US-align or verify tmtools is installed correctly |

<br>

---

<div align="center">

**Collateral Bio** -- Proprietary, 2026

*Built with precision for precision medicine.*

</div>
