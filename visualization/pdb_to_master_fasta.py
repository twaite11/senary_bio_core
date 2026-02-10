"""
Build a master FASTA from (1) sequences extracted from PDBs in structures/omegafold,
and (2) parent sequences from input_2-3_hepn.fasta. Use as --input for the structure
filter so every PDB gets an exact id match and parents are included.
"""
import argparse
from pathlib import Path
from typing import Optional

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def extract_sequence_from_pdb(pdb_path: str) -> Optional[str]:
    """Extract one-letter amino acid sequence from first chain of a PDB. Returns None on failure."""
    try:
        from tmtools.io import get_structure, get_residue_data
    except ImportError:
        return None
    try:
        s = get_structure(pdb_path)
        chain = next(s.get_chains(), None)
        if chain is None:
            return None
        coords, seq = get_residue_data(chain)
        return seq if seq else None
    except Exception:
        return None


def build_master_fasta(
    structures_dir: str,
    parents_fasta: str,
    output_fasta: str,
) -> int:
    """
    Write output_fasta with: all PDB-derived sequences (id=PDB stem) + all parent records from parents_fasta.
    Returns number of records written.
    """
    struct_dir = Path(structures_dir)
    parents_path = Path(parents_fasta)
    out_path = Path(output_fasta)

    records = []
    seen_ids = set()

    # 1. From PDBs: id = stem, seq = extracted
    if struct_dir.exists():
        for pdb in sorted(struct_dir.glob("*.pdb")):
            stem = pdb.stem
            seq = extract_sequence_from_pdb(str(pdb.resolve()))
            if seq:
                records.append(SeqRecord(Seq(seq), id=stem, description=""))
                seen_ids.add(stem)
            else:
                print(f"[!] Could not extract sequence from {pdb.name}, skipping.")
        print(f"[*] Extracted {len(records)} sequences from PDBs in {structures_dir}.")
    else:
        print(f"[!] Structures dir not found: {structures_dir}")

    # 2. Add parents from input_2-3_hepn (skip if id already from PDB)
    if parents_path.exists():
        n_added = 0
        for rec in SeqIO.parse(str(parents_path), "fasta"):
            if rec.id not in seen_ids:
                records.append(SeqRecord(rec.seq, id=rec.id, description=rec.description or ""))
                seen_ids.add(rec.id)
                n_added += 1
        print(f"[*] Added {n_added} parent sequences from {parents_fasta}.")
    else:
        print(f"[!] Parents FASTA not found: {parents_fasta}")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(records, str(out_path), "fasta")
    print(f"[+] Wrote {len(records)} records to {out_path}")
    return len(records)


def main():
    parser = argparse.ArgumentParser(
        description="Build master FASTA from PDBs + parent FASTA for structure filter."
    )
    parser.add_argument(
        "--structures-dir",
        default="data/structure_pipeline/structures/omegafold",
        help="Directory containing OmegaFold PDB files",
    )
    parser.add_argument(
        "--parents-fasta",
        default="data/structure_pipeline/input_2-3_hepn.fasta",
        help="FASTA of parent sequences (e.g. input_2-3_hepn.fasta)",
    )
    parser.add_argument(
        "--output",
        default="data/structure_pipeline/master_for_filter.fasta",
        help="Output master FASTA path",
    )
    args = parser.parse_args()

    n = build_master_fasta(args.structures_dir, args.parents_fasta, args.output)
    return 0 if n > 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
