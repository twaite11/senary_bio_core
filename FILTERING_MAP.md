# Collateral Bio: Full Filtering Map

**From NCBI Scraping → Final Type VI Cas13d Enzyme Candidates**

Wireframe map with filtering actions and criteria at each step.

---

## End-to-End Flow Overview

```
┌─────────────────────────────────────────────────────────────────────────────────────────────────┐
│  NCBI (Raw Data)  →  Mining Filters  →  Pipeline (Design+Structure+Identity)  →  Target Filters  │
│       →  Matchmaker  →  Expert Agent  →  FINAL                                                │
└─────────────────────────────────────────────────────────────────────────────────────────────────┘
```

---

## Phase 1: NCBI Scraping & Raw Collection

```
                                    ┌─────────────────────────────────────────┐
                                    │           NCBI DATABASES                 │
                                    └───────────────────┬─────────────────────┘
                                                        │
            ┌───────────────────────────────────────────┼───────────────────────────────────────────┐
            │                                           │                                           │
            ▼                                           ▼                                           ▼
┌───────────────────────────────┐     ┌───────────────────────────────┐     ┌───────────────────────────────┐
│       ncbi_miner              │     │       sra_scout               │     │   autonomous_prospector       │
├───────────────────────────────┤     ├───────────────────────────────┤     ├───────────────────────────────┤
│ FILTER: Search term match     │     │ FILTER: Environment keywords  │     │ FILTER: LLM picks strategy;   │
│ "Cas13d" in NCBI Protein DB   │     │ e.g. "hydrothermal vent       │     │ semantic filter picks top N   │
│ NCBI returns only annotated   │     │ metagenome". Try wgs[Prop]    │     │ datasets most likely to have   │
│ Cas13d proteins.              │     │ first; fallback to broader    │     │ uncultured CRISPR / viral      │
│                               │     │ search or BioProject elink.   │     │ defense / diverse microbes.    │
└───────────────┬───────────────┘     └───────────────┬───────────────┘     └───────────────┬───────────────┘
                │                                     │                                     │
                ▼                                     ▼                                     ▼
┌───────────────────────────────┐     ┌───────────────────────────────┐     ┌───────────────────────────────┐
│ search_*.fasta                │     │ undiscovered_cas13d_*.fasta   │     │ deep_hits_*.fasta             │
│ (annotated proteins)          │     │ (unannotated WGS hits)        │     │ (AI-validated deep hits)      │
│                               │     │ + full_orf_checks             │     │ + full_orf_checks + optional  │
│                               │     │                               │     │ CRISPR + ESM-2 similarity     │
└───────────────┬───────────────┘     └───────────────┬───────────────┘     └───────────────┬───────────────┘
                │                                     │                                     │
                └─────────────────────────────────────┼─────────────────────────────────────┘
                                                      │
                                                      ▼
                                    ┌─────────────────────────────────────────┐
                                    │   data/raw_sequences/                    │
                                    │   Raw FASTA files (input to Phase 2)    │
                                    └───────────────────┬─────────────────────┘
                                                        │
                                                        ▼
```

**Mining-time filters** (applied during collection; sra_scout & autonomous_prospector only for WGS paths):
- Size 600–1400 aa
- HEPN 2–3 motifs (R.{4,6}H)
- full_orf_checks: N-term M, C-term tail ≥15 aa after last HEPN, contig-boundary margin
- Prospector only (when REQUIRE_FULL_STRUCTURE): CRISPR array with ≥MIN_REPEAT_COUNT repeats; ESM-2 similarity (DeepEngine vs RfxCas13d/PspCas13a when ESM_REFERENCE_FASTA set)

---

## Phase 2: Enzyme Pipeline (run_full_pipeline.py)

The full pipeline takes raw FASTA (e.g. `deep_hits_*.fasta`, `family_grouped_*.fasta`) and runs: **Embed → Mutate for drift → Structure filter → Identity filter**. Matchmaker consumes the final output.

```
┌─────────────────────────────────────────────────────────────────────────────────────────────────────────┐
│                           ENZYME PIPELINE (run_full_pipeline.py)                                         │
└─────────────────────────────────────────────────────────────────────────────────────────────────────────┘

   RAW SEQUENCES (data/raw_sequences/deep_hits_*.fasta or family_grouped_*.fasta)
            │
            ▼
┌─────────────────────────────────────────────────────────────────────────────────────────────────────────┐
│ STEP 1: EMBED (optional; --skip-design skips)                                                            │
│ Action: ESM-2 embeddings for all sequences.                                                             │
│ Output: data/design/embeddings/                                                                         │
└─────────────────────────────────────────────────────────────────────────────────────────────────────────┘
            │
            ▼
┌─────────────────────────────────────────────────────────────────────────────────────────────────────────┐
│ STEP 2: MUTATE FOR DRIFT (optional; --skip-design skips)                                                 │
│ Keep: Variants scored for stability vs RfxCas13d/PspCas13a; keep <85% identity to known Cas13.           │
│ Throw out: Unstable or too similar to reference.                                                        │
│ Optional: Gemini trans-cleavage prompt for mutations that may increase activity.                        │
│ Output: data/design/drift_variants.fasta                                                                │
└─────────────────────────────────────────────────────────────────────────────────────────────────────────┘
            │
            ▼
┌─────────────────────────────────────────────────────────────────────────────────────────────────────────┐
│ STEP 3: STRUCTURE FILTER (optional; --skip-structure skips)                                              │
│ Keep: OmegaFold structure → TM-score vs Cas13 refs (5W1H, 6DTD, 6IV9) ≥ threshold (default 0.4);        │
│       2–3 HEPN motifs (R.{4,6}H) in sequence.                                                           │
│ Throw out: Low TM-score or wrong HEPN count.                                                            │
│ Output: data/structure_pipeline/passed_structures.fasta                                                 │
└─────────────────────────────────────────────────────────────────────────────────────────────────────────┘
            │
            ▼
┌─────────────────────────────────────────────────────────────────────────────────────────────────────────┐
│ STEP 4: IDENTITY FILTER (optional; --skip-identity skips)                                                │
│ Keep: Sequences with max identity to known_cas13.fasta &lt; 85% (drift goal).                              │
│ Throw out: Too similar to known Cas13 (Lwa, Rfx, etc.).                                                 │
│ Output: data/identity_filtered/passed.fasta, identity_metadata.csv                                      │
└─────────────────────────────────────────────────────────────────────────────────────────────────────────┘
            │
            ▼
┌─────────────────────────────────────────────────────────────────────────────────────────────────────────┐
│ FILTERED CAS13D ENZYMES → data/identity_filtered/passed.fasta                                           │
│ (Input to Matchmaker)                                                                                   │
└─────────────────────────────────────────────────────────────────────────────────────────────────────────┘
            │
            │  [Optional: Family Grouping on raw/mined sequences]
            ▼
┌─────────────────────────────────────────────────────────────────────────────────────────────────────────┐
│ FAMILY GROUPING (family_grouper.py – run separately)                                                    │
│ Partition by HEPN count → cluster by ESM-2 homology → name SN01_001, SN01_002, ...                      │
│ Output: data/mined_sequences/family_grouped_*.fasta, family_manifest_*.csv                              │
└─────────────────────────────────────────────────────────────────────────────────────────────────────────┘
```

---

## Phase 3: Target-Side Filtering

```
┌─────────────────────────────────────────────────────────────────────────────────────────────────────────┐
│                           TARGET FILTER PIPELINE                                                          │
└─────────────────────────────────────────────────────────────────────────────────────────────────────────┘

   ChimerDB / Excel:
   data/targets/known_fusions.csv, novel_fusions.csv
   data/matrices/disease_matrix_known.csv, disease_matrix_novel.csv
   KB_and_Pub_Recur_per_cancer.csv
            │
            ▼
┌─────────────────────────────────────────────────────────────────────────────────────────────────────────┐
│ FILTER 1: SPECIFICITY                                                                                    │
│ Keep: Fusions present in ≤3 cancer types (configurable).                                                │
│ Throw out: Fusions in many cancer types (promiscuous / artifact).                                       │
│ Why: Tissue-specific fusions are real drivers; pan-cancer = noise.                                      │
│ Output: data/targets/high_specificity_targets.csv                                                       │
└─────────────────────────────────────────────────────────────────────────────────────────────────────────┘
            │
            ▼
┌─────────────────────────────────────────────────────────────────────────────────────────────────────────┐
│ FILTER 2: FUSION METADATA (Mapping)                                                                      │
│ Action: Build fusion → [TCGA codes] and TCGA → organ keywords (e.g. PRAD→prostate, BRCA→breast).        │
│ Why: Needed for organ-specific ARCHS4 safety checks later.                                              │
└─────────────────────────────────────────────────────────────────────────────────────────────────────────┘
            │
            ▼
┌─────────────────────────────────────────────────────────────────────────────────────────────────────────┐
│ FILTERED FUSION TARGETS (cancer-specific, tissue-limited)                                               │
└─────────────────────────────────────────────────────────────────────────────────────────────────────────┘
```

---

## Phase 4: Matchmaker (Enzyme × Target Screening)

```
   FILTERED ENZYMES                         FILTERED TARGETS
   (data/identity_filtered/passed.fasta)    (high_specificity_targets.csv or known_fusions.csv)
   from run_full_pipeline --run-matchmaker  │
            │                                               │
            └───────────────────┬───────────────────────────┘
                                │
                                ▼
┌─────────────────────────────────────────────────────────────────────────────────────────────────────────┐
│ MATCHMAKER FILTER 1: VALID CUT SITES (PFS RULE)                                                         │
│ Keep: Pairs where fusion junction has valid cut sites (no G at position 23 of 22-base spacer).          │
│ Throw out: No valid cut sites.                                                                          │
│ Why: Cas13d won't cut well with G at 3′; more valid sites = more druggable.                             │
└─────────────────────────────────────────────────────────────────────────────────────────────────────────┘
                                │
                                ▼
┌─────────────────────────────────────────────────────────────────────────────────────────────────────────┐
│ MATCHMAKER FILTER 2: DISEASE ASSOCIATION                                                                │
│ Keep: Pairs where fusion has known Associated_Disease (from disease map).                               │
│ Throw out: Unknown or missing disease.                                                                  │
│ Why: We need to know which cancer we're targeting.                                                      │
└─────────────────────────────────────────────────────────────────────────────────────────────────────────┘
                                │
                                ▼
┌─────────────────────────────────────────────────────────────────────────────────────────────────────────┐
│ lead_candidates.csv (pre-safety, pre-AI)                                                                │
└─────────────────────────────────────────────────────────────────────────────────────────────────────────┘
```

---

## Phase 5: Expert Agent (Safety & AI Filter)

```
   lead_candidates.csv
            │
            ▼
┌─────────────────────────────────────────────────────────────────────────────────────────────────────────┐
│ EXPERT FILTER 1: ASSOCIATED DISEASE CHECK                                                               │
│ Keep: Candidates with valid Associated_Disease (cancer type).                                           │
│ Throw out: Blanks, "Unknown", "nan".                                                                    │
│ Why: Only pursue targets we can tie to a cancer.                                                        │
└─────────────────────────────────────────────────────────────────────────────────────────────────────────┘
            │
            ▼
┌─────────────────────────────────────────────────────────────────────────────────────────────────────────┐
│ EXPERT FILTER 2: ARCHS4 EXPRESSION (SAFETY)                                                             │
│ For each fusion's parent genes:                                                                         │
│ Keep: Low in normal tissue (max ≤ threshold), present in cancer (mean ≥ min).                           │
│       Organ-specific: cancer expression ≥ 2× normal (enrichment).                                       │
│ Throw out: High expression in healthy organs (risky) or absent in cancer.                               │
│ Why: Target should be "off" in normal, "on" in tumor.                                                   │
└─────────────────────────────────────────────────────────────────────────────────────────────────────────┘
            │
            ▼
┌─────────────────────────────────────────────────────────────────────────────────────────────────────────┐
│ EXPERT FILTER 3: GEMINI AI VERDICT                                                                      │
│ Action: AI returns GO / NO-GO / HOLD for each (fusion, disease) pair.                                   │
│ Considers: Target biology, patient count, cut sites, safety risks.                                      │
│ Why: First-pass prioritization before human review.                                                     │
└─────────────────────────────────────────────────────────────────────────────────────────────────────────┘
            │
            ▼
┌─────────────────────────────────────────────────────────────────────────────────────────────────────────┐
│ FINAL OUTPUT                                                                                            │
│ lead_candidates_filtered.csv  │  dashboard.html                                                         │
│ Type VI Cas13d enzyme–fusion pairs (cancer-specific, safety-validated)                                  │
└─────────────────────────────────────────────────────────────────────────────────────────────────────────┘
```

---

## Single-Page Wireframe Summary

```
╔══════════════════════════════════════════════════════════════════════════════════════════════════════════╗
║  COLLATERAL BIO: NCBI → TYPE VI CAS13D CANDIDATES                                                        ║
╚══════════════════════════════════════════════════════════════════════════════════════════════════════════╝

  ┌──────────────┐     ┌──────────────┐     ┌──────────────┐     ┌──────────────┐     ┌──────────────┐
  │    NCBI      │     │   PIPELINE   │     │   TARGET     │     │  MATCHMAKER  │     │   EXPERT     │
  │   SCRAPE     │────▶│   FILTERS    │────▶│   FILTERS    │────▶│   FILTERS    │────▶│   FILTERS    │
  └──────────────┘     └──────────────┘     └──────────────┘     └──────────────┘     └──────────────┘
       │                      │                     │                     │                     │
       ▼                      ▼                     ▼                     ▼                     ▼
  Search term          Embed → Mutate          Specificity            PFS cut sites        ARCHS4 safety
  Env keywords         Structure (OmegaFold   ≤3 cancer types        Disease known       Absent normal
  BioProject elink     + TM-score)             Fusion→organ           Rank by score       Present cancer
                       Identity (<85%)         mapping                                     Gemini GO/HOLD

  ═══════════════════════════════════════════════════════════════════════════════════════════════════════
  FINAL: lead_candidates_filtered.csv  │  Type VI Cas13d enzyme–fusion pairs (cancer-specific)
  ═══════════════════════════════════════════════════════════════════════════════════════════════════════
```

---

## Filter Summary Table (Quick Reference)

| Stage | Filter | Keep | Throw out |
|-------|--------|------|-----------|
| **NCBI** | Search | Matches query / env keywords | Rest of DB |
| **Mining** | Size | 600–1400 aa | Too short/long |
| **Mining** | HEPN | 2–3 motifs (R.{4,6}H) | Missing or 4+ domains |
| **Mining** | full_orf_checks | N-term M, C-term tail ≥15 aa, contig boundary | Fragments |
| **Mining** | CRISPR (Prospector) | ≥1 repeat when REQUIRE_FULL_STRUCTURE | No CRISPR context |
| **Mining** | ESM-2 (Prospector) | Similarity in range when ESM_REFERENCE_FASTA set | Outside range |
| **Pipeline** | Mutate for drift | Stable vs Rfx/Psp, <85% identity | Unstable, too similar |
| **Pipeline** | Structure | TM-score ≥0.4, 2–3 HEPN | Low TM-score, wrong HEPN |
| **Pipeline** | Identity | Max identity <85% to known_cas13.fasta | Too similar to known Cas13 |
| **Target** | Specificity | ≤3 cancer types | Promiscuous fusions |
| **Matchmaker** | PFS | Valid cut sites | No cut sites |
| **Matchmaker** | Disease | Known cancer | Unknown disease |
| **Expert** | ARCHS4 | Low normal, high cancer | Unsafe profile |
| **Expert** | Gemini | GO / HOLD | NO-GO |
