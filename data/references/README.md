# Reference sequences for mining and design

- **mining_refs.fasta** – Used when `ESM_REFERENCE_FASTA=data/references/mining_refs.fasta`. Should contain **RfxCas13d** and **PspCas13a** (or other Cas13 refs). Mining scores ESM similarity to each ref and keeps the **closest** match. For full-length scoring, paste full sequences from UniProt (RfxCas13d PDB 6IV9).
- Copy from `mining_refs.fasta.example` and replace RfxCas13d and PspCas13a with full sequences.
